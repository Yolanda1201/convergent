#!/usr/bin/env python3
"""
rerconverge_cli – batch launcher for RERconverge analyses
=========================================================

A modern command-line tool that wraps the RERconverge R package so that
large-scale, non-interactive analyses can be run with a single command.

Key features
------------
* Binary, continuous, and categorical trait analyses
* Automatic RER caching – expensive RER computation is done once and reused
  across all phenotype runs
* Multi-phenotype batch mode: pass a manifest TSV to test dozens of phenotypes
  in one invocation
* Optional parallel phenotype processing (phenotypes dispatched concurrently
  to Rscript workers)
* Rich progress reporting via tqdm (falls back gracefully if not installed)

Quick start
-----------
  # Single binary trait
  python rerconverge_cli.py run \\
      --trees  data/mammal_gene_trees.txt \\
      --phenotype data/marine_foreground.txt \\
      --trait-type binary \\
      --output results/marine.tsv

  # Batch mode (manifest TSV)
  python rerconverge_cli.py batch \\
      --trees data/mammal_gene_trees.txt \\
      --manifest data/phenotype_manifest.tsv \\
      --outdir results/ \\
      --rer-cache cache/mammal_RER.rds \\
      --jobs 4

Manifest format (tab-separated, no header required, but header is tolerated)
-----------------------------------------------------------------------------
  phenotype_name  phenotype_file      trait_type   [extra_flags …]
  marine          marine_fg.txt       binary
  body_size       body_size.tsv       continuous   --winsorize-rer 3
  diet            diet_categories.tsv categorical  --categorical-method kw
"""

import argparse
import os
import subprocess
import sys
import textwrap
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional

# Optional rich progress display
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

RSCRIPT = os.environ.get("RSCRIPT", "Rscript")
BATCH_R  = Path(__file__).parent / "scripts" / "rerconverge_batch.R"

_print_lock = threading.Lock()


def safe_print(*args, **kwargs):
    with _print_lock:
        print(*args, **kwargs)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def check_rscript():
    """Verify Rscript is available on PATH."""
    try:
        subprocess.run([RSCRIPT, "--version"], capture_output=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        sys.exit(
            f"[ERROR] '{RSCRIPT}' not found. Please install R or set the "
            "RSCRIPT environment variable to the full path of the Rscript executable."
        )


def check_batch_script():
    """Verify the companion R script exists."""
    if not BATCH_R.exists():
        sys.exit(
            f"[ERROR] Companion R script not found at: {BATCH_R}\n"
            "Make sure you are running rerconverge_cli.py from the repository root."
        )


def build_rscript_cmd(args_ns, extra_flags: Optional[List[str]] = None) -> List[str]:
    """
    Translate parsed argparse namespace (or raw kwargs dict) into the
    Rscript command list for rerconverge_batch.R.
    """
    a = vars(args_ns) if not isinstance(args_ns, dict) else args_ns

    cmd = [RSCRIPT, "--vanilla", str(BATCH_R)]

    def add(flag, key, store_true=False):
        val = a.get(key)
        if store_true:
            if val:
                cmd.append(flag)
        elif val is not None:
            cmd.extend([flag, str(val)])

    add("--trees",              "trees")
    add("--phenotype",          "phenotype")
    add("--output",             "output")
    add("--trait-type",         "trait_type")
    add("--clade",              "clade")
    add("--transition",         "transition")
    add("--weighted-pheno",     "weighted_pheno",     store_true=True)
    add("--categorical-method", "categorical_method")
    add("--transform",          "transform")
    add("--weighted-rer",       "weighted_rer",       store_true=True)
    add("--scale",              "scale",              store_true=True)
    add("--use-species",        "use_species")
    add("--max-read",           "max_read")
    add("--min-sp",             "min_sp")
    add("--min-pos",            "min_pos")
    add("--winsorize-rer",      "winsorize_rer")
    add("--winsorize-trait",    "winsorize_trait")
    add("--rer-cache",          "rer_cache")
    add("--trees-cache",        "trees_cache")
    add("--no-plots",           "no_plots",           store_true=True)
    add("--sort",               "sort",               store_true=True)
    add("--verbose",            "verbose",            store_true=True)

    if extra_flags:
        cmd.extend(extra_flags)

    return cmd


def run_rscript(cmd: List[str], label: str = "", verbose: bool = False) -> int:
    """Run Rscript command; stream stdout/stderr when verbose."""
    safe_print(f"[START] {label or ' '.join(cmd[:4])}")
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if verbose or proc.returncode != 0:
        with _print_lock:
            for line in proc.stdout.splitlines():
                print(f"  {label}: {line}")
    if proc.returncode != 0:
        safe_print(f"[FAIL]  {label} — exit code {proc.returncode}")
    else:
        safe_print(f"[OK]    {label}")
    return proc.returncode


# ---------------------------------------------------------------------------
# Sub-command: run
# ---------------------------------------------------------------------------

def cmd_run(args):
    check_rscript()
    check_batch_script()
    cmd = build_rscript_cmd(args)
    rc  = run_rscript(cmd, label=os.path.basename(args.output), verbose=args.verbose)
    sys.exit(rc)


# ---------------------------------------------------------------------------
# Sub-command: batch
# ---------------------------------------------------------------------------

def parse_manifest(manifest_path: str):
    """
    Parse a phenotype manifest TSV.

    Columns (tab-separated):
      1. phenotype_name  – used to name the output file
      2. phenotype_file  – path to phenotype data
      3. trait_type      – binary | continuous | categorical
      4. extra_flags     – (optional) additional flags forwarded verbatim to R script

    Lines starting with '#' and empty lines are ignored.
    A header line starting with 'phenotype_name' or 'name' is tolerated.
    """
    jobs = []
    with open(manifest_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if parts[0].lower() in ("phenotype_name", "name", "phenotype"):
                continue  # skip header
            if len(parts) < 3:
                sys.exit(
                    f"[ERROR] Manifest line has fewer than 3 columns:\n  {line}\n"
                    "Expected: name <TAB> phenotype_file <TAB> trait_type [<TAB> extra_flags]"
                )
            name  = parts[0].strip()
            pfile = parts[1].strip()
            ttype = parts[2].strip()
            extra = parts[3].strip().split() if len(parts) >= 4 else []
            if not os.path.exists(pfile):
                sys.exit(f"[ERROR] Phenotype file listed in manifest not found: {pfile}")
            jobs.append({"name": name, "phenotype": pfile, "trait_type": ttype, "extra": extra})
    return jobs


def cmd_batch(args):
    check_rscript()
    check_batch_script()

    jobs = parse_manifest(args.manifest)
    if not jobs:
        sys.exit("[ERROR] No valid entries found in manifest.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Build shared base args dict
    base = {
        "trees":              args.trees,
        "transform":          args.transform,
        "weighted_rer":       args.weighted_rer,
        "scale":              args.scale,
        "use_species":        getattr(args, "use_species", None),
        "max_read":           getattr(args, "max_read", None),
        "min_sp":             args.min_sp,
        "min_pos":            args.min_pos,
        "winsorize_rer":      getattr(args, "winsorize_rer", None),
        "winsorize_trait":    getattr(args, "winsorize_trait", None),
        "rer_cache":          args.rer_cache,
        "trees_cache":        args.trees_cache,
        "no_plots":           args.no_plots,
        "sort":               args.sort,
        "verbose":            args.verbose,
        "clade":              args.clade,
        "transition":         args.transition,
        "weighted_pheno":     args.weighted_pheno,
        "categorical_method": args.categorical_method,
    }

    def make_cmd(job):
        d = dict(base)
        d["phenotype"]  = job["phenotype"]
        d["trait_type"] = job["trait_type"]
        d["output"]     = str(outdir / f"{job['name']}.tsv")
        return build_rscript_cmd(d, extra_flags=job.get("extra", []))

    # Dispatch jobs
    n_jobs = args.jobs
    failed = []

    if n_jobs == 1:
        iterator = jobs if not HAS_TQDM else tqdm(jobs, desc="Phenotypes")
        for job in iterator:
            cmd = make_cmd(job)
            rc  = run_rscript(cmd, label=job["name"], verbose=args.verbose)
            if rc != 0:
                failed.append(job["name"])
    else:
        with ThreadPoolExecutor(max_workers=n_jobs) as pool:
            futs = {pool.submit(run_rscript, make_cmd(j), j["name"], args.verbose): j["name"]
                    for j in jobs}
            if HAS_TQDM:
                with tqdm(total=len(futs), desc="Phenotypes") as pbar:
                    for fut in as_completed(futs):
                        name = futs[fut]
                        rc   = fut.result()
                        if rc != 0:
                            failed.append(name)
                        pbar.update(1)
            else:
                for fut in as_completed(futs):
                    name = futs[fut]
                    rc   = fut.result()
                    if rc != 0:
                        failed.append(name)

    if failed:
        print(f"\n[WARN] {len(failed)} phenotype(s) failed: {', '.join(failed)}", file=sys.stderr)
        sys.exit(1)
    else:
        print(f"\n[DONE] All {len(jobs)} phenotype(s) completed successfully.")
        print(f"       Results written to: {outdir}/")


# ---------------------------------------------------------------------------
# Sub-command: install-deps
# ---------------------------------------------------------------------------

def cmd_install_deps(args):
    """Install R packages required by RERconverge."""
    r_code = textwrap.dedent(r"""
        pkgs <- c("optparse", "devtools", "ape", "phytools", "compiler",
                  "plotrix", "Rcpp", "RcppArmadillo", "phangorn", "weights",
                  "castor", "FSA", "Matrix", "data.table", "rsvd",
                  "progress", "TreeTools", "impute")

        missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
        if (length(missing)) {
          message("Installing missing CRAN packages: ", paste(missing, collapse = ", "))
          install.packages(missing, repos = "https://cloud.r-project.org", quiet = FALSE)
        } else {
          message("All CRAN dependencies already installed.")
        }

        if (!requireNamespace("RERconverge", quietly = TRUE)) {
          message("Installing RERconverge from GitHub …")
          devtools::install_github("nclark-lab/RERconverge")
        } else {
          message("RERconverge already installed.")
        }
    """)
    cmd = [RSCRIPT, "--vanilla", "-e", r_code]
    print("Installing R dependencies …")
    rc = subprocess.run(cmd).returncode
    if rc == 0:
        print("Dependencies installed successfully.")
    else:
        print("Some installations may have failed – check output above.", file=sys.stderr)
        sys.exit(rc)


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def common_rer_args(p):
    """Add RER computation arguments shared between run and batch sub-commands."""
    g = p.add_argument_group("RER computation")
    g.add_argument("--trees", required=True, metavar="FILE",
                   help="Gene trees file (tab-delimited: gene_name<TAB>newick_tree)")
    g.add_argument("--transform", default="sqrt", choices=["none", "sqrt", "log"],
                   help="Branch-length transform [default: sqrt]")
    g.add_argument("--weighted-rer", action="store_true", default=True,
                   help="Use weighted regression for RER estimation [default: True]")
    g.add_argument("--scale", action="store_true", default=True,
                   help="Scale gene-tree branches before RER [default: True]")
    g.add_argument("--use-species", metavar="SP1,SP2,…",
                   help="Comma-separated species to include (optional)")
    g.add_argument("--max-read", type=int, metavar="N",
                   help="Maximum number of gene trees to read (useful for testing)")
    g.add_argument("--rer-cache", metavar="FILE",
                   help="RDS file to save/load cached RER matrix")
    g.add_argument("--trees-cache", metavar="FILE",
                   help="RDS file to save/load cached treesObj")
    return p


def common_filter_args(p):
    g = p.add_argument_group("Filters")
    g.add_argument("--min-sp", type=int, default=10, metavar="N",
                   help="Min species per gene tree [default: 10]")
    g.add_argument("--min-pos", type=int, default=2, metavar="N",
                   help="Min foreground lineages per gene [default: 2]")
    g.add_argument("--winsorize-rer", type=int, metavar="N",
                   help="Winsorise N extreme RER values per gene (continuous: default 3)")
    g.add_argument("--winsorize-trait", type=int, metavar="N",
                   help="Winsorise N extreme trait values (continuous: default 3)")
    return p


def common_output_args(p):
    g = p.add_argument_group("Output")
    g.add_argument("--no-plots", action="store_true",
                   help="Suppress RERconverge diagnostic plots")
    g.add_argument("--sort", action="store_true",
                   help="Sort results by signed -log10(P)")
    g.add_argument("--verbose", action="store_true",
                   help="Print R output to terminal")
    return p


def common_phenotype_args(p):
    g = p.add_argument_group("Phenotype options")
    g.add_argument("--clade", default="ancestral",
                   choices=["ancestral", "terminal", "all"],
                   help="Binary: which clade lineages to mark foreground [default: ancestral]")
    g.add_argument("--transition", default="unidirectional",
                   choices=["unidirectional", "bidirectional"],
                   help="Binary: direction of trait transitions [default: unidirectional]")
    g.add_argument("--weighted-pheno", action="store_true",
                   help="Binary: distribute clade weight evenly across branches")
    g.add_argument("--categorical-method", default="kw", choices=["kw", "aov"],
                   help="Categorical: kw (Kruskal-Wallis) or aov (ANOVA) [default: kw]")
    return p


def build_parser():
    root = argparse.ArgumentParser(
        prog="rerconverge_cli.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = root.add_subparsers(dest="subcmd", metavar="COMMAND")
    sub.required = True

    # ---- run ----------------------------------------------------------------
    p_run = sub.add_parser("run", help="Run a single phenotype analysis",
                           formatter_class=argparse.RawDescriptionHelpFormatter)
    common_rer_args(p_run)
    common_phenotype_args(p_run)
    common_filter_args(p_run)
    common_output_args(p_run)
    p_run.add_argument("--phenotype", required=True, metavar="FILE",
                       help="Phenotype file (format depends on --trait-type)")
    p_run.add_argument("--trait-type", default="binary",
                       choices=["binary", "continuous", "categorical"],
                       help="Trait type [default: binary]")
    p_run.add_argument("--output", default="rerconverge_results.tsv", metavar="FILE",
                       help="Output TSV file [default: rerconverge_results.tsv]")
    p_run.set_defaults(func=cmd_run)

    # ---- batch --------------------------------------------------------------
    p_batch = sub.add_parser("batch", help="Run multiple phenotypes from a manifest",
                             formatter_class=argparse.RawDescriptionHelpFormatter,
                             description=textwrap.dedent("""
                               Run many phenotypes from a TSV manifest.

                               Manifest columns (tab-separated):
                                 1. name          – label for this phenotype (used for output filename)
                                 2. phenotype_file – path to phenotype data file
                                 3. trait_type    – binary | continuous | categorical
                                 4. extra_flags   – (optional) additional R-script flags
                             """))
    common_rer_args(p_batch)
    common_phenotype_args(p_batch)
    common_filter_args(p_batch)
    common_output_args(p_batch)
    p_batch.add_argument("--manifest", required=True, metavar="FILE",
                         help="Phenotype manifest TSV file")
    p_batch.add_argument("--outdir", default="rerconverge_results", metavar="DIR",
                         help="Output directory for per-phenotype TSV files [default: rerconverge_results]")
    p_batch.add_argument("--jobs", type=int, default=1, metavar="N",
                         help="Parallel phenotype jobs [default: 1]")
    p_batch.set_defaults(func=cmd_batch)

    # ---- install-deps -------------------------------------------------------
    p_inst = sub.add_parser("install-deps", help="Install R package dependencies")
    p_inst.set_defaults(func=cmd_install_deps)

    return root


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = build_parser()
    args   = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
