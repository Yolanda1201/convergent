# convergent — RERconverge Batch Pipeline

A **non-interactive, high-throughput batch pipeline** built on top of the
[RERconverge](https://github.com/nclark-lab/RERconverge) R package for
genome-wide association of evolutionary rate shifts with convergent phenotypes.

RERconverge itself is designed for interactive R sessions.  This repository
wraps its core computation functions so that you can:

* Run a **single phenotype** with one command (no R console needed)
* Run **hundreds of phenotypes** in batch from a manifest file
* **Cache** the expensive Relative Evolutionary Rate (RER) computation so it is
  done once and reused across all phenotype runs
* Dispatch multiple phenotype analyses **in parallel**

---

## Supported analyses

| Trait type | RERconverge function used | Input |
|------------|---------------------------|-------|
| **Binary** | `correlateWithBinaryPhenotype` | Text file with foreground species names *or* a Newick trait-tree |
| **Continuous** | `correlateWithContinuousPhenotype` | Two-column TSV: species, value |
| **Categorical** | `correlateWithCategoricalPhenotype` | Two-column TSV: species, category |

---

## Installation

### 1. Clone this repository

```bash
git clone <this-repo-url>
cd convergent
```

### 2. Set up the R environment

Either use the provided conda environment (recommended):

```bash
conda env create -f environment.yml
conda activate rerconverge-batch
```

Or install R packages manually inside R:

```r
install.packages(c(
  "optparse", "devtools", "ape", "phytools", "phangorn",
  "weights", "FSA", "data.table", "rsvd", "progress",
  "TreeTools", "impute", "Rcpp", "RcppArmadillo"
))
devtools::install_github("nclark-lab/RERconverge")
```

Or use the built-in helper:

```bash
python rerconverge_cli.py install-deps
```

### 3. Install Python dependencies (optional, for the CLI wrapper)

```bash
pip install -r requirements.txt
```

---

## Quick start

### Single binary trait run

```bash
python rerconverge_cli.py run \
  --trees  data/mammal_gene_trees.txt \
  --phenotype data/marine_foreground.txt \
  --trait-type binary \
  --clade ancestral \
  --output results/marine.tsv \
  --no-plots --verbose
```

### Single continuous trait run

```bash
python rerconverge_cli.py run \
  --trees  data/mammal_gene_trees.txt \
  --phenotype data/adult_weight.tsv \
  --trait-type continuous \
  --rer-cache cache/mammal_RER.rds \
  --output results/body_size.tsv
```

### Batch mode — many phenotypes at once

```bash
python rerconverge_cli.py batch \
  --trees  data/mammal_gene_trees.txt \
  --manifest data/phenotype_manifest.tsv \
  --outdir results/ \
  --rer-cache cache/mammal_RER.rds \
  --jobs 8          # run 8 phenotypes in parallel
```

### Use the R script directly (without the Python wrapper)

```bash
Rscript scripts/rerconverge_batch.R \
  --trees  data/mammal_gene_trees.txt \
  --phenotype data/marine_foreground.txt \
  --trait-type binary \
  --clade ancestral \
  --output results/marine.tsv \
  --no-plots
```

---

## Input formats

### Gene trees file (`--trees`)

Tab-delimited, one gene per line:

```
GENE1   (Species_A:0.1,Species_B:0.2,(Species_C:0.05,Species_D:0.07):0.15);
GENE2   (Species_A:0.08,...
```

This is the same format as `extdata/subsetMammalGeneTrees.txt` bundled with
RERconverge.

### Binary phenotype — foreground species list

Plain text, one species per line.  Lines starting with `#` are ignored:

```
# examples/data/marine_foreground.txt
Walrus
Seal
Killer_whale
Dolphin
Manatee
```

### Binary phenotype — Newick trait tree

A Newick tree file where foreground branches have length > 0 and background
branches have length 0.  The tool auto-detects which format you supplied based
on whether the first non-comment character is `(`.

### Continuous phenotype

Two-column TSV, no header (comments with `#` are tolerated):

```
Human       4.65
Chimpanzee  4.28
Dolphin     4.88
…
```

### Categorical phenotype

Two-column TSV, no header:

```
Dog        carnivore
Cow        herbivore
Human      omnivore
…
```

### Batch manifest (`--manifest`)

Tab-delimited, minimum three columns.  A header row is tolerated.  Lines
starting with `#` are ignored:

```tsv
# name <TAB> phenotype_file <TAB> trait_type [<TAB> extra_R_flags]
marine      examples/data/marine_foreground.txt      binary
body_size   examples/data/adult_weight_real.tsv    continuous  --winsorize-rer 3
diet        examples/data/diet_categories.tsv      categorical --categorical-method kw
```

---

## Output

Each run produces a **tab-separated results file** with one row per gene:

| Column | Description |
|--------|-------------|
| `gene` | Gene name (from the trees file) |
| `Rho`  | Correlation coefficient (Pearson / Kendall / Spearman) |
| `N`    | Number of branches used in the test |
| `P`    | Uncorrected p-value |
| `p.adj`| Benjamini–Hochberg adjusted p-value (FDR) |

For categorical analyses an additional `_pairwise.rds` file is saved alongside
the main TSV, containing per-category-pair Dunn test statistics.

---

## All options

### `rerconverge_cli.py run` / `Rscript scripts/rerconverge_batch.R`

| Flag | Default | Description |
|------|---------|-------------|
| `--trees FILE` | required | Gene trees file |
| `--phenotype FILE` | required | Phenotype data file |
| `--output FILE` | `rerconverge_results.tsv` | Output TSV |
| `--trait-type` | `binary` | `binary` \| `continuous` \| `categorical` |
| `--clade` | `ancestral` | Binary: `ancestral` \| `terminal` \| `all` |
| `--transition` | `unidirectional` | Binary: direction of transitions |
| `--weighted-pheno` | off | Binary: distribute clade weight evenly |
| `--categorical-method` | `kw` | Categorical: `kw` (Kruskal-Wallis) \| `aov` |
| `--transform` | `sqrt` | RER transform: `none` \| `sqrt` \| `log` |
| `--weighted-rer` | on | Weighted regression for RER |
| `--scale` | on | Scale gene-tree branches |
| `--use-species SP1,SP2,…` | all | Comma-separated species subset |
| `--max-read N` | all | Max gene trees to read (testing) |
| `--min-sp N` | `10` | Min species per gene |
| `--min-pos N` | `2` | Min foreground lineages per gene |
| `--winsorize-rer N` | `3` (continuous) | Winsorise N extreme RER values |
| `--winsorize-trait N` | `3` (continuous) | Winsorise N extreme trait values |
| `--rer-cache FILE` | none | Save/load RER matrix RDS |
| `--trees-cache FILE` | none | Save/load treesObj RDS |
| `--no-plots` | off | Suppress diagnostic plots |
| `--sort` | off | Sort by signed −log10(P) |
| `--verbose` / `-v` | off | Print extra progress |

### `rerconverge_cli.py batch` (additional)

| Flag | Default | Description |
|------|---------|-------------|
| `--manifest FILE` | required | Phenotype manifest TSV |
| `--outdir DIR` | `rerconverge_results` | Output directory |
| `--jobs N` | `1` | Parallel phenotype workers |

---

## Caching strategy

The two most expensive steps are:

1. **Reading gene trees** (`readTrees`) — O(n_genes × n_species)
2. **Computing RERs** (`getAllResiduals`) — O(n_genes × n_branches)

Both are cached with `--trees-cache` and `--rer-cache` respectively.  On a
subsequent run with a different phenotype, supply the same cache paths and both
steps will be skipped entirely.

```bash
# First phenotype: compute and cache everything
python rerconverge_cli.py run \
  --trees gene_trees.txt \
  --phenotype marine.txt --trait-type binary \
  --trees-cache cache/trees.rds \
  --rer-cache  cache/rer.rds \
  --output results/marine.tsv

# All subsequent phenotypes: load from cache (near-instant)
python rerconverge_cli.py run \
  --trees gene_trees.txt \
  --phenotype flight.txt --trait-type binary \
  --trees-cache cache/trees.rds \
  --rer-cache  cache/rer.rds \
  --output results/flight.tsv
```

---

## Lineage-specific acceleration (branch acceleration analysis)

RERconverge detects *convergent* rate shifts.  For a targeted
**lineage-specific accelerated evolution** test, use a binary trait tree where
only the lineage(s) of interest are set as foreground.  Example:

```text
# primate_ancestor_only.txt  — Newick trait tree
# foreground = single internal branch (primate ancestor)
((Human:0,(Chimpanzee:0,Gorilla:0):0):1,(Mouse:0,Rat:0):0);
```

Then:

```bash
python rerconverge_cli.py run \
  --trees gene_trees.txt \
  --phenotype primate_ancestor_only.txt \
  --trait-type binary \
  --output results/primate_acceleration.tsv
```

Genes with significant positive `Rho` and low `p.adj` are accelerated on that
lineage relative to genome-wide expectations.

---

## Citation

If you use this pipeline, please cite the original RERconverge publication:

> Kowalczyk A, Meyer WK, Partha R, Mao W, Clark NL, Chikina M.
> RERconverge: an R package for associating evolutionary rates with convergent traits.
> *Bioinformatics* 2019. <https://doi.org/10.1093/bioinformatics/btz468>

---

## License

GPL-3.0 (same as RERconverge)
