#!/usr/bin/env Rscript
# =============================================================================
# RERconverge Batch Runner
# =============================================================================
# Non-interactive, command-line driver for the RERconverge R package.
# Supports binary, continuous and categorical trait analyses; caches the
# expensive RER computation step so multiple phenotypes can be tested against
# the same gene-tree set without re-computing from scratch.
#
# Usage (see --help):
#   Rscript rerconverge_batch.R [options]
#
# Required packages: RERconverge, optparse, parallel, data.table
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE))
    install.packages("optparse", repos = "https://cloud.r-project.org", quiet = TRUE)
  library(optparse)
})

# ---------------------------------------------------------------------------
# CLI argument definitions
# ---------------------------------------------------------------------------
option_list <- list(

  # ---- Input ----------------------------------------------------------------
  make_option(c("-t", "--trees"),
    type = "character", default = NULL,
    help = "Path to gene trees file (tab-delimited: gene_name <TAB> newick_tree). REQUIRED."),

  make_option(c("-p", "--phenotype"),
    type = "character", default = NULL,
    help = paste0("Path to phenotype file. Format depends on --trait-type:\n",
      "  binary     : one foreground species name per line  OR  a Newick trait tree file\n",
      "  continuous : two-column TSV (species <TAB> value), no header\n",
      "  categorical: two-column TSV (species <TAB> category_label), no header\n",
      "REQUIRED.")),

  make_option(c("-o", "--output"),
    type = "character", default = "rerconverge_results.tsv",
    help = "Output TSV file for correlation results [default: %default]"),

  # ---- Analysis type --------------------------------------------------------
  make_option(c("--trait-type"),
    type = "character", default = "binary",
    help = "Trait type: binary | continuous | categorical [default: %default]"),

  make_option(c("--clade"),
    type = "character", default = "ancestral",
    help = "For binary: which lineages mark as foreground: ancestral | terminal | all [default: %default]"),

  make_option(c("--transition"),
    type = "character", default = "unidirectional",
    help = "For binary: unidirectional | bidirectional [default: %default]"),

  make_option(c("--weighted-pheno"),
    action = "store_true", default = FALSE,
    help = "For binary/all: distribute clade weight evenly across branches (foreground2Tree weighted=TRUE)"),

  make_option(c("--categorical-method"),
    type = "character", default = "kw",
    help = "For categorical: kw (Kruskal-Wallis/Dunn) | aov (ANOVA) [default: %default]"),

  # ---- RER computation ------------------------------------------------------
  make_option(c("--transform"),
    type = "character", default = "sqrt",
    help = "Branch-length transform: none | sqrt | log [default: %default]"),

  make_option(c("--weighted-rer"),
    action = "store_true", default = TRUE,
    help = "Use weighted regression for RER estimation [default: TRUE]"),

  make_option(c("--scale"),
    action = "store_true", default = TRUE,
    help = "Scale gene-tree branches before RER computation [default: TRUE]"),

  make_option(c("--use-species"),
    type = "character", default = NULL,
    help = "Comma-separated list of species to include (others excluded)"),

  make_option(c("--max-read"),
    type = "integer", default = NA,
    help = "Maximum number of gene trees to read (useful for testing) [default: all]"),

  # ---- Filters --------------------------------------------------------------
  make_option(c("--min-sp"),
    type = "integer", default = 10,
    help = "Minimum number of species per gene tree [default: %default]"),

  make_option(c("--min-pos"),
    type = "integer", default = 2,
    help = "Minimum foreground lineages per gene (binary/categorical) [default: %default]"),

  make_option(c("--winsorize-rer"),
    type = "integer", default = NULL,
    help = "Winsorise this many extreme RER values at each end per gene (continuous) [default: 3 if continuous]"),

  make_option(c("--winsorize-trait"),
    type = "integer", default = NULL,
    help = "Winsorise this many extreme trait values (continuous) [default: 3 if continuous]"),

  # ---- Caching --------------------------------------------------------------
  make_option(c("--rer-cache"),
    type = "character", default = NULL,
    help = "Path to save/load RER matrix RDS cache. If the file exists it is loaded instead of recomputed."),

  make_option(c("--trees-cache"),
    type = "character", default = NULL,
    help = "Path to save/load treesObj RDS cache. If the file exists it is loaded instead of re-read."),

  # ---- Parallel -------------------------------------------------------------
  make_option(c("--cores"),
    type = "integer", default = 1,
    help = "Number of CPU cores to use where RERconverge supports parallelism [default: %default]"),

  # ---- Misc -----------------------------------------------------------------
  make_option(c("--no-plots"),
    action = "store_true", default = FALSE,
    help = "Suppress all diagnostic plots generated by getAllResiduals"),

  make_option(c("--sort"),
    action = "store_true", default = FALSE,
    help = "Sort output by signed -log10(P)"),

  make_option(c("-v", "--verbose"),
    action = "store_true", default = FALSE,
    help = "Print extra progress messages")
)

parser <- OptionParser(
  usage       = "Rscript %prog [options]",
  option_list = option_list,
  description = paste0(
    "\nNon-interactive batch runner for RERconverge.\n",
    "Computes relative evolutionary rates (RER) genome-wide and correlates\n",
    "them with a user-supplied phenotype without any interactive prompts.\n\n",
    "Examples:\n",
    "  # Binary trait (foreground species list)\n",
    "  Rscript rerconverge_batch.R \\\n",
    "    --trees mammal_gene_trees.txt \\\n",
    "    --phenotype marine_foreground.txt \\\n",
    "    --trait-type binary --clade ancestral \\\n",
    "    --output marine_results.tsv\n\n",
    "  # Continuous trait with pre-computed RER cache\n",
    "  Rscript rerconverge_batch.R \\\n",
    "    --trees mammal_gene_trees.txt \\\n",
    "    --phenotype adult_weight.tsv \\\n",
    "    --trait-type continuous \\\n",
    "    --rer-cache mammal_RER.rds \\\n",
    "    --output weight_results.tsv\n"
  )
)

opt <- parse_args(parser)

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
if (is.null(opt$trees))     stop("--trees is required. Run with --help for usage.")
if (is.null(opt$phenotype)) stop("--phenotype is required. Run with --help for usage.")
if (!file.exists(opt$trees))     stop(paste("Trees file not found:", opt$trees))
if (!file.exists(opt$phenotype)) stop(paste("Phenotype file not found:", opt$phenotype))

trait_type <- tolower(opt[["trait-type"]])
if (!trait_type %in% c("binary", "continuous", "categorical"))
  stop("--trait-type must be one of: binary, continuous, categorical")

verbose <- opt$verbose
msg <- function(...) if (verbose) message(paste0("[RERconverge] ", ...))

# ---------------------------------------------------------------------------
# Load RERconverge
# ---------------------------------------------------------------------------
msg("Loading RERconverge package …")
if (!requireNamespace("RERconverge", quietly = TRUE)) {
  message("RERconverge not found – attempting installation from GitHub …")
  if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools", repos = "https://cloud.r-project.org", quiet = TRUE)
  devtools::install_github("nclark-lab/RERconverge", quiet = TRUE)
}
suppressPackageStartupMessages(library(RERconverge))

# Suppress diagnostic plots when requested
if (opt[["no-plots"]]) {
  # redirect graphics to a null device for the duration of the script
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
}

# ---------------------------------------------------------------------------
# Species subset
# ---------------------------------------------------------------------------
use_species <- NULL
if (!is.null(opt[["use-species"]])) {
  use_species <- trimws(strsplit(opt[["use-species"]], ",")[[1]])
  msg("Species subset: ", paste(use_species, collapse = ", "))
}

# ---------------------------------------------------------------------------
# Step 1 – Read gene trees (with optional cache)
# ---------------------------------------------------------------------------
trees_cache_path <- opt[["trees-cache"]]

if (!is.null(trees_cache_path) && file.exists(trees_cache_path)) {
  message("Loading treesObj from cache: ", trees_cache_path)
  treesObj <- readRDS(trees_cache_path)
} else {
  message("Reading gene trees from: ", opt$trees, " …")
  max_read <- if (is.na(opt[["max-read"]])) NA else opt[["max-read"]]
  treesObj <- readTrees(
    opt$trees,
    max.read   = max_read,
    useSpecies = use_species
  )
  if (!is.null(trees_cache_path)) {
    message("Saving treesObj cache to: ", trees_cache_path)
    saveRDS(treesObj, file = trees_cache_path)
  }
}
message("Gene trees loaded: ", treesObj$numTrees, " trees, ",
        treesObj$maxSp, " species")

# ---------------------------------------------------------------------------
# Step 2 – Compute RERs (with optional cache)
# ---------------------------------------------------------------------------
rer_cache_path <- opt[["rer-cache"]]

if (!is.null(rer_cache_path) && file.exists(rer_cache_path)) {
  message("Loading RER matrix from cache: ", rer_cache_path)
  mamRERw <- readRDS(rer_cache_path)
} else {
  message("Computing relative evolutionary rates (RERs) …")
  # getAllResiduals signature varies between RERconverge versions.
  # Inspect at runtime and only pass arguments that are accepted.
  rer_formals <- names(formals(getAllResiduals))

  rer_args <- list(treesObj = treesObj)
  if (!is.null(use_species) && "useSpecies" %in% rer_formals)
    rer_args$useSpecies <- use_species
  if ("transform" %in% rer_formals)
    rer_args$transform <- opt$transform
  # Older API: weighted / scale / plot
  if ("weighted" %in% rer_formals)
    rer_args$weighted <- opt[["weighted-rer"]]
  if ("scale" %in% rer_formals)
    rer_args$scale <- opt$scale
  if ("plot" %in% rer_formals)
    rer_args$plot <- !opt[["no-plots"]]
  # Newer API (v0.3+): use.weights instead of weighted; norm instead of scale
  if ("use.weights" %in% rer_formals)
    rer_args[["use.weights"]] <- opt[["weighted-rer"]]
  if ("norm" %in% rer_formals) {
    rer_args$norm <- if (opt$scale) "scale" else "none"
  }

  mamRERw <- do.call(getAllResiduals, rer_args)
  if (!is.null(rer_cache_path)) {
    message("Saving RER cache to: ", rer_cache_path)
    saveRDS(mamRERw, file = rer_cache_path)
  }
}
message("RER matrix dimensions: ", nrow(mamRERw), " genes × ", ncol(mamRERw), " branches")

# ---------------------------------------------------------------------------
# Step 3 – Build phenotype paths
# ---------------------------------------------------------------------------
message("Building phenotype paths for trait type: ", trait_type, " …")

build_binary_paths <- function() {
  # Detect whether the phenotype file is a Newick tree or a species list
  first_line <- readLines(opt$phenotype, n = 1, warn = FALSE)
  is_newick  <- grepl("^\\(", trimws(first_line))

  if (is_newick) {
    msg("Phenotype file detected as Newick trait tree")
    trait_tree <- read.tree(opt$phenotype)
    phenv      <- tree2Paths(trait_tree, treesObj, binaryCol = TRUE)
  } else {
    msg("Phenotype file detected as foreground species list")
    fg_species <- trimws(readLines(opt$phenotype, warn = FALSE))
    fg_species <- fg_species[nchar(fg_species) > 0]
    msg("Foreground species: ", paste(fg_species, collapse = ", "))
    fp_formals <- names(formals(foreground2Paths))
    fp_args <- list(
      foreground = fg_species,
      treesObj   = treesObj,
      clade      = opt$clade,
      plotTree   = !opt[["no-plots"]]
    )
    if (!is.null(use_species) && "useSpecies" %in% fp_formals)
      fp_args$useSpecies <- use_species
    if ("transition" %in% fp_formals)
      fp_args$transition <- opt$transition
    if ("weighted" %in% fp_formals)
      fp_args$weighted <- opt[["weighted-pheno"]]
    phenv <- do.call(foreground2Paths, fp_args)
  }
  phenv
}

build_continuous_paths <- function() {
  df <- read.table(opt$phenotype, header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE, col.names = c("species", "value"))
  trait_vec        <- as.numeric(df$value)
  names(trait_vec) <- df$species
  msg("Continuous trait: ", length(trait_vec), " species loaded")
  char2Paths(trait_vec, treesObj)
}

build_categorical_paths <- function() {
  df <- read.table(opt$phenotype, header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE, col.names = c("species", "category"))
  cat_vec        <- df$category
  names(cat_vec) <- df$species
  msg("Categorical trait: ", length(unique(cat_vec)), " categories, ",
      length(cat_vec), " species loaded")
  char2PathsCategorical(cat_vec, treesObj)
}

phenv <- switch(trait_type,
  binary     = build_binary_paths(),
  continuous = build_continuous_paths(),
  categorical = build_categorical_paths()
)

# ---------------------------------------------------------------------------
# Step 4 – Correlation / association test
# ---------------------------------------------------------------------------
message("Running association test …")

results <- switch(trait_type,

  binary = correlateWithBinaryPhenotype(
    mamRERw, phenv,
    min.sp         = opt[["min-sp"]],
    min.pos        = opt[["min-pos"]],
    weighted       = if (opt[["weighted-pheno"]]) TRUE else "auto",
    winsorizeRER   = opt[["winsorize-rer"]],
    winsorizetrait = opt[["winsorize-trait"]],
    sort           = opt$sort
  ),

  continuous = {
    wRER   <- if (is.null(opt[["winsorize-rer"]]))   3L else opt[["winsorize-rer"]]
    wTrait <- if (is.null(opt[["winsorize-trait"]])) 3L else opt[["winsorize-trait"]]
    correlateWithContinuousPhenotype(
      mamRERw, phenv,
      min.sp         = opt[["min-sp"]],
      winsorizeRER   = wRER,
      winsorizetrait = wTrait,
      sort           = opt$sort
    )
  },

  categorical = {
    res_list <- correlateWithCategoricalPhenotype(
      mamRERw, phenv,
      min.sp         = opt[["min-sp"]],
      min.pos        = opt[["min-pos"]],
      method         = opt[["categorical-method"]],
      winsorizeRER   = opt[["winsorize-rer"]],
      winsorizetrait = opt[["winsorize-trait"]],
      sort           = opt$sort
    )
    # res_list[[1]] is the main corout table, res_list[[2]] is pairwise tables
    # save pairwise tables alongside main output
    pw_path <- sub("(\\.tsv|\\.txt)?$", "_pairwise.rds", opt$output, ignore.case = TRUE)
    message("Saving pairwise categorical tables to: ", pw_path)
    saveRDS(res_list[[2]], file = pw_path)
    res_list[[1]]
  }
)

# ---------------------------------------------------------------------------
# Step 5 – Write results
# ---------------------------------------------------------------------------
message("Writing results to: ", opt$output)

# Ensure gene names are present as a column
out_df <- as.data.frame(results)
out_df <- cbind(gene = rownames(out_df), out_df)
rownames(out_df) <- NULL

write.table(out_df, file = opt$output, sep = "\t", row.names = FALSE, quote = FALSE)

message("Done. ", nrow(out_df), " genes written to ", opt$output)
