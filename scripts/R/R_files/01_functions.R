# 01_functions.R
# All custom functions used in the ambient RNA benchmarking pipeline.
# Source this file at the top of each analysis script:
#   source(here::here("scripts/R/R_files/01_functions.R"))

# Add project-local R library to search path
.libPaths(c(here::here("dependencies/r_packages"), .libPaths()))


# Preprocessing and quality control ---------------------------------------

#' Load a STARsolo / CellRanger output folder
#'
#' Reads barcodes.tsv, features.tsv, and matrix.mtx from a directory and
#' returns a sparse count matrix.
#'
#' @param folder Character. Path to the folder containing the three files.
#'
#' @return A sparse CsparseMatrix (genes x barcodes).
load_folder <- function(folder) {
  barcodes <- read.delim(file.path(folder, "barcodes.tsv"), header = FALSE)
  features <- read.delim(file.path(folder, "features.tsv"), header = FALSE)
  counts   <- Matrix::readMM(file.path(folder, "matrix.mtx"))
  colnames(counts) <- barcodes$V1
  rownames(counts) <- features$V2
  counts <- as(counts, "CsparseMatrix")
  return(counts)
}


#' Combine two single-cell RNA-seq count matrices
#'
#' Combines two scRNA-seq count matrices, handling cases where gene lists
#' differ between runs. Suffixes are appended to cell names to ensure
#' uniqueness in the combined matrix.
#'
#' @param counts1 First count matrix (genes x cells).
#' @param counts2 Second count matrix (genes x cells).
#' @param suffix1 Character. Suffix for cell names in counts1 (default "-run1").
#' @param suffix2 Character. Suffix for cell names in counts2 (default "-run2").
#' @param sparse Logical. Return a sparse matrix (default TRUE).
#'
#' @return A combined count matrix with all cells from both inputs.
combine_count_matrices <- function(counts1, counts2,
                                   suffix1 = "-run1",
                                   suffix2 = "-run2",
                                   sparse  = TRUE) {
  if (!is.matrix(counts1) && !inherits(counts1, "dgCMatrix"))
    stop("counts1 must be a matrix or sparse matrix")
  if (!is.matrix(counts2) && !inherits(counts2, "dgCMatrix"))
    stop("counts2 must be a matrix or sparse matrix")

  colnames(counts1) <- paste0(colnames(counts1), suffix1)
  colnames(counts2) <- paste0(colnames(counts2), suffix2)

  if (identical(rownames(counts1), rownames(counts2))) {
    message("Gene lists are identical — performing direct cbind.")
    return(cbind(counts1, counts2))
  }

  message("Gene lists differ — performing union operation.")
  all_genes <- union(rownames(counts1), rownames(counts2))

  if (sparse && (inherits(counts1, "dgCMatrix") || inherits(counts2, "dgCMatrix"))) {
    counts1 <- as(counts1, "dgCMatrix")
    counts2 <- as(counts2, "dgCMatrix")
    counts1_full <- Matrix::Matrix(0, nrow = length(all_genes), ncol = ncol(counts1),
                                   sparse = TRUE, dimnames = list(all_genes, colnames(counts1)))
    counts2_full <- Matrix::Matrix(0, nrow = length(all_genes), ncol = ncol(counts2),
                                   sparse = TRUE, dimnames = list(all_genes, colnames(counts2)))
  } else {
    counts1_full <- matrix(0, nrow = length(all_genes), ncol = ncol(counts1),
                           dimnames = list(all_genes, colnames(counts1)))
    counts2_full <- matrix(0, nrow = length(all_genes), ncol = ncol(counts2),
                           dimnames = list(all_genes, colnames(counts2)))
  }

  counts1_full[rownames(counts1), ] <- counts1
  counts2_full[rownames(counts2), ] <- counts2
  cbind(counts1_full, counts2_full)
}


# Generate synthetic datasets ---------------------------------------------

#' Generate a synthetic single-cell RNA-seq dataset
#'
#' Simulates a dataset with real cells, ambient RNA contamination, and empty
#' droplets, based on negative-binomial / Poisson models fit to real cell-line
#' data. Ambient fraction per cell is drawn from a clipped normal distribution.
#'
#' @param param_list Named list of model parameter sets (one per cell line),
#'   as returned by \code{fit_models()}.
#' @param sizes Numeric vector of library sizes from real data, as returned
#'   by \code{get_real_sizes()}.
#' @param ambient_mean Numeric. Mean ambient fraction (default 0.2).
#' @param ambient_sd Numeric. SD of ambient fraction (default 0.05).
#' @param mix_fraction Numeric vector. Relative proportions of cell types
#'   (default equal proportions). Set to NULL for equal mixing.
#' @param threshold Integer. Min UMIs for a droplet to be considered a real
#'   cell (default 1000).
#' @param group_size Integer. Cells simulated per group (default 20).
#' @param n_cpu Integer. CPU cores for parallel simulation (default 60).
#' @param min_count Integer. Min UMIs for an empty droplet (default 10).
#' @param ambient_max Numeric. Maximum ambient fraction cap (default 0.8).
#' @param ambient_min Numeric. Minimum ambient fraction floor (default 0.01).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{empty}{Sparse matrix of empty droplet counts.}
#'     \item{real_counts}{Sparse matrix of endogenous counts per cell.}
#'     \item{real_ambient}{Sparse matrix of ambient counts per cell.}
#'     \item{metadata}{Data frame with cell type, class, and ambient fraction.}
#'   }
generate_dataset <- function(param_list,
                              sizes,
                              ambient_mean = 0.2,
                              ambient_sd   = 0.05,
                              mix_fraction = NULL,
                              threshold    = 1000,
                              group_size   = 20,
                              n_cpu        = 60,
                              min_count    = 10,
                              ambient_max  = 0.8,
                              ambient_min  = 0.01) {

  real  <- sizes[sizes >= threshold]
  empty <- sizes[sizes < threshold & sizes > min_count]

  # Ambient fraction per cell (clipped normal)
  ambient_fraction <- rnorm(length(real), mean = ambient_mean, sd = ambient_sd)
  ambient_fraction <- pmax(ambient_min, pmin(ambient_max, ambient_fraction))

  # Mixing proportions
  if (is.null(mix_fraction)) mix_fraction <- rep(1, length(param_list))
  mix_fraction <- mix_fraction / sum(mix_fraction)

  # Sample cell type labels
  celltypes <- sample(names(param_list), size = length(real),
                      replace = TRUE, prob = mix_fraction)

  # Superset of genes across all cell lines
  genes <- unique(unlist(lapply(param_list, function(p) rownames(p$params))))

  # Simulate real cells
  metadata      <- data.frame()
  ambient_list  <- list()
  real_list     <- list()

  for (celltype in unique(celltypes)) {
    ct_idx    <- which(celltypes == celltype)
    n_groups  <- ceiling(length(ct_idx) / group_size)
    ct_model  <- which(names(param_list) == celltype)

    adjusted_sizes <- floor(
      real[ct_idx] * ((1 - ambient_fraction[ct_idx]) +
                        (ambient_fraction[ct_idx] * mix_fraction[ct_model]))
    )
    ct_simulated <- simulate_cells(param_list[[ct_model]], n_groups, adjusted_sizes)
    ct_simulated <- insert_missing_genes(ct_simulated, genes)
    ct_simulated[is.na(ct_simulated)] <- 0

    # Simulate ambient from the other cell lines
    amb_simulated <- Matrix::Matrix(0, nrow = nrow(ct_simulated),
                                    ncol = ncol(ct_simulated), sparse = TRUE)
    for (model in which(names(param_list) != celltype)) {
      adj <- floor(real[ct_idx] * (ambient_fraction[ct_idx] * mix_fraction[model]))
      tmp <- simulate_cells(param_list[[model]], n_groups, adj)
      tmp <- insert_missing_genes(tmp, genes)
      tmp[is.na(tmp)] <- 0
      amb_simulated <- amb_simulated + tmp
    }

    metadata <- rbind(metadata, data.frame(
      celltype      = rep(celltype, length(ct_idx)),
      class         = "real",
      ambient_counts = colSums(amb_simulated),
      real_counts   = colSums(ct_simulated)
    ))
    real_list[[length(real_list)   + 1]] <- ct_simulated
    ambient_list[[length(ambient_list) + 1]] <- amb_simulated
  }

  real_counts        <- do.call("cbind", real_list)
  real_ambient_counts <- do.call("cbind", ambient_list)
  metadata$ambient_fraction <- metadata$ambient_counts /
    (metadata$ambient_counts + metadata$real_counts)

  # Simulate empty droplets
  empty_simulated <- Matrix::Matrix(0, nrow = length(genes),
                                    ncol = length(empty), sparse = TRUE)
  n_groups <- ceiling(length(empty) / group_size)

  for (model in seq_along(param_list)) {
    adj <- floor(empty * mix_fraction[model])
    adj[adj == 0] <- 1
    tmp <- simulate_cells(param_list[[model]], n_groups, adj)
    tmp <- insert_missing_genes(tmp, genes)
    empty_simulated <- empty_simulated + tmp
  }

  list(
    empty       = empty_simulated,
    real_counts = real_counts,
    real_ambient = real_ambient_counts,
    metadata    = metadata
  )
}


#' Pad a count matrix with zero rows for missing genes
#'
#' @param count_mat Sparse count matrix to pad.
#' @param genes Character vector of all expected gene names.
#'
#' @return Padded and reordered sparse count matrix.
insert_missing_genes <- function(count_mat, genes) {
  missing_genes  <- genes[!(genes %in% rownames(count_mat))]
  missing_counts <- Matrix::Matrix(0, nrow = length(missing_genes),
                                   ncol = ncol(count_mat), sparse = TRUE,
                                   dimnames = list(missing_genes, colnames(count_mat)))
  count_mat <- rbind(count_mat, missing_counts)
  count_mat[match(genes, rownames(count_mat)), ]
}


#' Get per-cell library sizes from a raw STARsolo output folder
#'
#' @param folder Character. Path to a Solo.out/Gene/raw folder.
#'
#' @return Named numeric vector of total UMI counts per barcode.
get_real_sizes <- function(folder) {
  counts   <- Matrix::readMM(file.path(folder, "matrix.mtx"))
  barcodes <- read.delim(file.path(folder, "barcodes.tsv"), header = FALSE)
  features <- read.delim(file.path(folder, "features.tsv"), header = FALSE)

  keep              <- Matrix::colSums(counts) >= 1
  counts            <- counts[, keep]
  barcodes          <- barcodes[keep, 1]
  rownames(counts)  <- features[, 1]
  colnames(counts)  <- barcodes

  Matrix::colSums(counts)
}


#' Fit negative-binomial / Poisson models to a count matrix
#'
#' Genes with variance > mean are fit with a negative-binomial model; all
#' other genes use a Poisson model.
#'
#' @param counts Count matrix (genes x cells) from a filtered cell-line dataset.
#' @param min_mean Numeric. Genes with mean expression below this are modelled
#'   as Poisson regardless of variance (default 0.01).
#'
#' @return A named list with elements \code{params} (data frame of fitted
#'   parameters) and \code{sf_mean} (geometric mean of size factors).
fit_models <- function(counts, min_mean = 0.01) {
  sf      <- colSums(counts)
  sf_mean <- exp(mean(log(sf)))

  gene_mean <- rowMeans(counts)
  gene_var  <- apply(counts, 1, var)
  nb_genes  <- names(which(gene_var > gene_mean & gene_mean >= min_mean))
  po_genes  <- rownames(counts)[!(rownames(counts) %in% nb_genes)]

  nb_model  <- glmGamPoi::glm_gp(counts[nb_genes, ], design = ~ 1,
                                   size_factors = sf / sf_mean)
  po_model  <- glmGamPoi::glm_gp(counts[po_genes, ], design = ~ 1,
                                   overdispersion = FALSE, size_factors = sf / sf_mean)

  model_params <- rbind(
    data.frame(beta = nb_model$Beta[, 1], theta = nb_model$overdispersions, model = "NB"),
    data.frame(beta = po_model$Beta[, 1], theta = 0,                        model = "Poisson")
  )
  model_params <- model_params[match(rownames(counts), rownames(model_params)), ]
  list(params = model_params, sf_mean = sf_mean)
}


#' Filter a cell-line dataset for simulation use
#'
#' @param folder Character. Path to a Solo.out/Gene/raw folder.
#' @param min_counts Integer. Minimum total UMI count per cell.
#' @param min_cells Integer. Minimum number of cells a gene must appear in.
#' @param max_mito Numeric. Maximum mitochondrial fraction (0–1).
#'
#' @return A dense count matrix of filtered cells and genes.
filter_dataset <- function(folder, min_counts, min_cells, max_mito) {
  counts   <- Matrix::readMM(file.path(folder, "matrix.mtx"))
  barcodes <- read.delim(file.path(folder, "barcodes.tsv"), header = FALSE)
  features <- read.delim(file.path(folder, "features.tsv"), header = FALSE)

  keep     <- Matrix::colSums(counts) >= min_counts
  counts   <- counts[, keep]
  barcodes <- barcodes[keep, 1]

  keep     <- Matrix::rowSums(counts) >= min_cells
  counts   <- counts[keep, ]
  features <- features[keep, ]

  rownames(counts) <- features[, 1]
  colnames(counts) <- barcodes

  # Mitochondrial filtering
  mt_idx   <- grep("^MT-", features$V2)
  mt_frac  <- Matrix::colSums(counts[features[mt_idx, 1], ]) / Matrix::colSums(counts)
  counts   <- counts[, mt_frac <= max_mito]

  as.matrix(counts)
}


#' Simulate cells from a fitted model
#'
#' @param params Named list with elements \code{params} (model parameters) and
#'   \code{sf_mean} (size factor mean), as returned by \code{fit_models()}.
#' @param n_groups Integer. Number of groups to split simulation into.
#' @param counts Numeric vector of target library sizes per cell.
#' @param n_cpu Integer. CPU cores for parallel execution (default 16).
#'
#' @return A sparse count matrix (genes x cells).
simulate_cells <- function(params, n_groups, counts, n_cpu = 16) {
  model_params <- params$params
  sf_mean      <- params$sf_mean

  counts     <- counts / sf_mean
  group_list <- split(counts, rep(seq_len(n_groups),
                                  each = length(counts) %/% n_groups,
                                  length.out = length(counts)))

  mat_list <- BiocParallel::bplapply(
    group_list,
    BPPARAM = BiocParallel::MulticoreParam(workers = n_cpu),
    FUN = function(x) {
      mat        <- matrix(ncol = length(x), nrow = nrow(model_params))
      is_negbinom <- which(model_params[, "theta"] > 0)
      is_poisson  <- which(model_params[, "theta"] == 0)

      for (i in seq_along(x)) {
        offset <- x[i]
        mu     <- model_params[is_negbinom, "beta"]
        theta  <- 1 / (model_params[is_negbinom, "theta"] *
                         (exp(mu) / exp(mu + log(offset))))
        mu     <- exp(mu + log(offset))
        mat[is_negbinom, i] <- mapply(stats::rnbinom, mu = mu, size = theta,
                                      MoreArgs = list(n = 1))
        mu_p  <- exp(model_params[is_poisson, "beta"] + log(offset))
        mat[is_poisson, i]  <- mapply(stats::rpois, lambda = mu_p,
                                      MoreArgs = list(n = 1))
      }
      mat
    }
  )

  mat           <- do.call("cbind", mat_list)
  rownames(mat) <- rownames(model_params)
  colnames(mat) <- names(counts)
  Matrix::Matrix(mat, sparse = TRUE)
}


#' Build raw count matrices from simulated samples (with ambient RNA)
#'
#' @param simulated_samples Named list of simulation objects from
#'   \code{generate_dataset()}.
#'
#' @return Named list of sparse count matrices (real cells + empty droplets).
create_count_matrices <- function(simulated_samples) {
  lapply(simulated_samples, function(s) {
    cbind(s$real_counts + s$real_ambient, s$empty)
  })
}


#' Build raw count matrices from simulated samples (without ambient RNA)
#'
#' @param simulated_samples Named list as returned by \code{generate_dataset()}.
#'
#' @return Named list of sparse count matrices (real counts only + empty).
create_count_matrices_wo_ambient <- function(simulated_samples) {
  lapply(simulated_samples, function(s) cbind(s$real_counts, s$empty))
}


#' Build filtered count matrices (with ambient RNA)
#'
#' @param simulated_samples Named list as returned by \code{generate_dataset()}.
#'
#' @return Named list of filtered sparse count matrices.
create_filtered_count_matrices <- function(simulated_samples) {
  lapply(simulated_samples, function(s) s$real_counts + s$real_ambient)
}


#' Build filtered count matrices (without ambient RNA)
#'
#' @param simulated_samples Named list as returned by \code{generate_dataset()}.
#'
#' @return Named list of filtered sparse count matrices.
create_filtered_count_matrices_wo_ambient <- function(simulated_samples) {
  lapply(simulated_samples, function(s) s$real_counts)
}


#' Extract filtered barcodes from filtered count matrices
#'
#' @param filtered_count_matrices Named list of filtered count matrices.
#'
#' @return Named list of character vectors of cell barcodes.
get_filtered_barcodes <- function(filtered_count_matrices) {
  lapply(filtered_count_matrices, colnames)
}


#' Combine and return metadata from all simulated samples
#'
#' @param simulated_samples Named list as returned by \code{generate_dataset()}.
#'
#' @return A data frame of cell metadata with sample-tagged rownames.
get_metadata <- function(simulated_samples) {
  metadata_list <- lapply(names(simulated_samples), function(sample) {
    md            <- simulated_samples[[sample]]$metadata
    rownames(md)  <- paste0(rownames(md), "_", sample)
    md
  })
  dplyr::bind_rows(metadata_list)
}


# Running decontamination methods -----------------------------------------

## Memory and timing helpers -----------------------------------------------

#' Get current R heap usage in MB (via gc)
get_gc_mb <- function() gc(reset = FALSE)["Vcells", "used"]

#' Get current process RSS in MB (via ps)
get_rss_mb <- function() {
  as.numeric(system("ps -o rss= -p $$", intern = TRUE)) / 1024
}

#' Log memory use at a labelled checkpoint
#'
#' @param label Character. Label for this checkpoint.
#' @param file Character. Path to the log file (appended to).
#' @param before_gc Numeric. gc MB at the previous checkpoint (for delta).
#' @param before_rss Numeric. RSS MB at the previous checkpoint (for delta).
#'
#' @return Named list with elements \code{gc} and \code{rss} for chaining.
log_mem <- function(label, file, before_gc = NULL, before_rss = NULL) {
  gc_now  <- get_gc_mb()
  rss_now <- get_rss_mb()
  delta_gc  <- if (!is.null(before_gc))  gc_now  - before_gc  else NA
  delta_rss <- if (!is.null(before_rss)) rss_now - before_rss else NA

  cat(
    "=== ", label, " ===\n",
    "Time:         ", format(Sys.time()), "\n",
    "GC_MB:        ", gc_now,   "\n",
    "RSS_MB:       ", rss_now,  "\n",
    "Delta_GC_MB:  ", delta_gc,  "\n",
    "Delta_RSS_MB: ", delta_rss, "\n\n",
    file = file, append = TRUE
  )
  list(gc = gc_now, rss = rss_now)
}


## R-based methods ---------------------------------------------------------

#' Run all R-based ambient RNA decontamination methods
#'
#' Clusters cells interactively (prompts for PCA dimensions and clustering
#' resolution after displaying elbow plots and clustree/UMAP panels), then
#' runs FastCAR, SoupX (full and reduced), scCDC, and DecontX (full and
#' reduced) on each sample. Running times and memory use are logged to file.
#'
#' @param counts Named list of raw count matrices (genes x barcodes), one
#'   per sample.
#' @param barcodes Named list of character vectors giving called cell
#'   barcodes, one per sample.
#' @param resolutions Numeric vector of resolutions to test in clustree
#'   (default 0.1 – 1.0 in steps of 0.1).
#' @param nfeat Integer. Variable features for PCA (default 2000).
#' @param empty_cutoff_upper Numeric. Upper UMI boundary for empty droplets
#'   (passed to FastCAR and SoupX). NULL = use defaults.
#' @param empty_cutoff_lower Numeric. Lower UMI boundary (SoupX soupRange).
#'   NULL = use defaults.
#' @param dataset_name Character. Label used for output file names.
#'
#' @return A Seurat object containing the RNA assay plus one assay per
#'   method: FastCAR, scCDC, SoupXfull, SoupXreduced, DecontXfull,
#'   DecontXreduced.
run_decontamination_methods <- function(counts,
                                        barcodes,
                                        resolutions        = seq(0.1, 1.0, by = 0.1),
                                        nfeat              = 2000,
                                        empty_cutoff_upper = NULL,
                                        empty_cutoff_lower = NULL,
                                        dataset_name       = NULL) {

  message("Setting up data.")

  # ------------------------------------------------------------------
  # 1. Subset and build per-sample Seurat objects
  # ------------------------------------------------------------------
  seu_list    <- list()
  subset_list <- list()

  for (dataset in seq_along(counts)) {
    colnames(counts[[dataset]]) <- paste(
      colnames(counts[[dataset]]), names(counts)[dataset], sep = "_"
    )
    barcodes[[dataset]] <- paste(barcodes[[dataset]], names(counts)[dataset], sep = "_")

    subset_list[[dataset]] <- counts[[dataset]][
      , colnames(counts[[dataset]]) %in% barcodes[[dataset]]
    ]
    subset_list[[dataset]] <- subset_list[[dataset]][
      Matrix::rowSums(subset_list[[dataset]]) >= 1, 
    ]
    names(subset_list)[dataset]          <- names(counts)[dataset]
    rownames(subset_list[[dataset]])     <- make.unique(rownames(subset_list[[dataset]]))

    seu_list[[dataset]] <- suppressWarnings(
      Seurat::CreateSeuratObject(
        counts       = Seurat::CreateAssayObject(subset_list[[dataset]]),
        project      = names(counts)[dataset],
        min.cells    = 0,
        min.features = 0
      )
    )
    seu_list[[dataset]]@meta.data$barcode <- colnames(subset_list[[dataset]])
    seu_list[[dataset]] <- subset(
      seu_list[[dataset]],
      cells = rownames(seu_list[[dataset]]@meta.data)[
        seu_list[[dataset]]$nFeature_RNA >= 200
      ]
    )
    seu_list[[dataset]] <- subset(
      seu_list[[dataset]],
      features = names(which(
        rowSums(seu_list[[dataset]]@assays$RNA$counts > 0) >= 3
      ))
    )
  }

  # Merge
  if (length(counts) > 1) {
    seu <- suppressWarnings(merge(seu_list[[1]], seu_list[2:length(seu_list)]))
  } else {
    seu <- seu_list[[1]]
  }

  # Join layers after merging (required for Seurat v5 multi-sample objects)
  seu <- JoinLayers(seu)
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, verbose = FALSE, nfeatures = nfeat)
  seu <- Seurat::ScaleData(seu, features = rownames(seu), verbose = FALSE)
  seu <- Seurat::RunPCA(seu, verbose = FALSE)

  # ------------------------------------------------------------------
  # 2. Interactive dimension selection
  # ------------------------------------------------------------------
  message("Visualising PCs.")
  print(Seurat::DimHeatmap(seu, dims = 1:20, cells = 500, balanced = TRUE))
  print(Seurat::ElbowPlot(seu))

  dims <- as.numeric(readline("Select number of PCs based on elbow plot: "))
  seu@misc$dims_used <- dims

  # Save diagnostic PDFs
  dir.create(here::here("results/dim_tests"), recursive = TRUE, showWarnings = FALSE)
  pdf(here::here(paste0("results/dim_tests/", dataset_name, "_dim_test.pdf")),
      width = 14, height = 10)
  print(Seurat::ElbowPlot(seu))
  print(Seurat::DimHeatmap(seu, dims = 1:20, cells = 500, balanced = TRUE))
  dev.off()
  message("Saved elbow plot PDF.")

  if (length(counts) > 1) {
    seu <- suppressWarnings(harmony::RunHarmony(seu, group.by.vars = "orig.ident",
                                                verbose = FALSE))
    seu <- Seurat::FindNeighbors(seu, dims = 1:dims, reduction = "harmony",
                                 verbose = FALSE)
  } else {
    seu <- Seurat::FindNeighbors(seu, dims = 1:dims, verbose = FALSE)
  }

  # ------------------------------------------------------------------
  # 3. Interactive resolution selection
  # ------------------------------------------------------------------
  message("Testing clustering resolutions.")
  for (res in resolutions) {
    seu <- Seurat::FindClusters(seu, resolution = res, verbose = FALSE)
  }
  if (!("umap" %in% names(seu@reductions))) {
    seu <- Seurat::RunUMAP(seu, dims = 1:dims, verbose = FALSE)
  }

  umap_plots <- lapply(resolutions, function(res) {
    rname <- paste0("RNA_snn_res.", res)
    suppressMessages(
      Seurat::DimPlot(seu, reduction = "umap", group.by = rname, label = TRUE) +
        ggplot2::ggtitle(paste("Resolution:", res))
    )
  })
  ct <- clustree::clustree(seu, prefix = "RNA_snn_res.")
  print(patchwork::wrap_plots(umap_plots, ncol = 3))
  print(ct)

  dir.create(here::here("results/res_tests"), recursive = TRUE, showWarnings = FALSE)
  pdf(here::here(paste0("results/res_tests/", dataset_name, "_res_test.pdf")),
      width = 14, height = 10)
  print(patchwork::wrap_plots(umap_plots, ncol = 3))
  print(ct)
  dev.off()
  message("Saved resolution test PDF.")

  chosen_resolution <- as.numeric(readline("Enter resolution to use: "))
  seu@misc$resolution_used <- chosen_resolution
  seu <- Seurat::FindClusters(seu, resolution = chosen_resolution, verbose = FALSE)

  # Species fraction metadata
  seu[["input.percent.mm10"]]   <- Seurat::PercentageFeatureSet(seu, pattern = "^mm10")
  seu[["input.percent.GRCh38"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^GRCh38")

  # Re-split into per-sample objects
  for (dataset in seq_along(counts)) {
    seu_list[[dataset]] <- subset(
      seu,
      cells = rownames(seu@meta.data)[seu$orig.ident == names(counts)[dataset]]
    )
  }

  # ------------------------------------------------------------------
  # 4. Run decontamination methods
  # ------------------------------------------------------------------
  fastcar         <- list()
  cdc             <- list()
  soupx_full      <- list()
  soupx_reduced   <- list()
  decontx_full    <- list()
  decontx_reduced <- list()

  running_times <- data.frame(tool_and_dataset = character(),
                               time_seconds     = numeric(),
                               stringsAsFactors = FALSE)

  dir.create(here::here("results/memory_use"), recursive = TRUE, showWarnings = FALSE)
  log_file <- here::here(paste0(
    "results/memory_use/memory_log_", dataset_name, "_",
    format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"
  ))

  for (dataset in seq_along(counts)) {
    filt_matrix <- Seurat::GetAssayData(seu_list[[dataset]], assay = "RNA", layer = "counts")
    colnames(filt_matrix) <- seu_list[[dataset]]$barcode

    raw_matrix            <- counts[[dataset]]
    rownames(raw_matrix)  <- gsub("_", "-", rownames(raw_matrix))
    raw_matrix            <- raw_matrix[!duplicated(rownames(raw_matrix)), ]
    rownames(filt_matrix) <- gsub("_", "-", rownames(filt_matrix))
    filt_matrix           <- filt_matrix[!duplicated(rownames(filt_matrix)), ]
    filt_matrix           <- filt_matrix[rownames(filt_matrix) %in% rownames(raw_matrix), ]
    raw_matrix            <- raw_matrix[rownames(raw_matrix) %in% rownames(filt_matrix), ]
    filt_matrix           <- filt_matrix[match(rownames(raw_matrix), rownames(filt_matrix)), ]

    clusters <- data.frame(
      barcode = seu_list[[dataset]]$barcode,
      cluster = seu_list[[dataset]]$seurat_clusters
    )
    clusters <- clusters[match(colnames(filt_matrix), clusters$barcode), ]

    # Save inputs for memory profiling (low complexity only)
    if (!is.null(dataset_name) && dataset_name == "low_complexity") {
      dir.create(here::here("data/processed/memory_use"), recursive = TRUE,
                 showWarnings = FALSE)
      saveRDS(raw_matrix,   here::here("data/processed/memory_use/lc_raw_matrix.Rds"))
      saveRDS(filt_matrix,  here::here("data/processed/memory_use/lc_filt_matrix.Rds"))
      saveRDS(clusters,     here::here("data/processed/memory_use/lc_clusters.Rds"))
      saveRDS(seu_list,     here::here("data/processed/memory_use/lc_seu_list.Rds"))
    }

    message("\nRunning decontamination for: ", names(counts)[dataset])

    # --- FastCAR ---
    message("\tRunning FastCAR.")
    bef <- log_mem(paste("Before FastCAR", names(counts)[dataset]), log_file)
    tictoc::tic()
    cutoff <- if (!is.null(empty_cutoff_upper)) empty_cutoff_upper else 100
    ambientProfile <- FastCAR::determine.background.to.remove(
      raw_matrix, filt_matrix,
      emptyDropletCutoff        = cutoff,
      contaminationChanceCutoff = 0.05
    )
    fastcar[[dataset]] <- FastCAR::remove.background(filt_matrix, ambientProfile)
    t_fc <- tictoc::toc(quiet = TRUE)
    log_mem(paste("After FastCAR", names(counts)[dataset]), log_file, bef$gc, bef$rss)
    running_times <- rbind(running_times, data.frame(
      tool_and_dataset = paste0("Dataset", dataset, "_FastCAR"),
      time_seconds     = t_fc$toc - t_fc$tic
    ))

    # --- SoupX full ---
    message("\tRunning SoupX.")
    bef <- log_mem(paste("Before SoupX", names(counts)[dataset]), log_file)
    tictoc::tic()
    if (is.null(empty_cutoff_upper)) {
      sc <- SoupX::SoupChannel(raw_matrix, filt_matrix)
    } else {
      sc <- SoupX::SoupChannel(raw_matrix, filt_matrix, calcSoupProfile = FALSE)
      sc <- SoupX::estimateSoup(sc, soupRange = c(empty_cutoff_lower, empty_cutoff_upper))
    }
    sc <- SoupX::setClusters(sc, setNames(clusters$cluster, clusters$barcode))
    sc <- SoupX::autoEstCont(sc, forceAccept = TRUE, doPlot = FALSE, verbose = FALSE)
    soupx_full[[dataset]] <- suppressWarnings(SoupX::adjustCounts(sc, verbose = FALSE))
    t_sx <- tictoc::toc(quiet = TRUE)
    log_mem(paste("After SoupX", names(counts)[dataset]), log_file, bef$gc, bef$rss)
    running_times <- rbind(running_times, data.frame(
      tool_and_dataset = paste0("Dataset", dataset, "_SoupX_full"),
      time_seconds     = t_sx$toc - t_sx$tic
    ))

    # --- SoupX reduced ---
    message("\tRunning SoupX (filtered only).")
    bef <- log_mem(paste("Before SoupX_filt", names(counts)[dataset]), log_file)
    tictoc::tic()
    sc_r      <- SoupX::SoupChannel(filt_matrix, filt_matrix, calcSoupProfile = FALSE)
    sc_r      <- SoupX::setClusters(sc_r, setNames(clusters$cluster, clusters$barcode))
    soup_prof <- data.frame(
      row.names = rownames(filt_matrix),
      est       = Matrix::rowSums(filt_matrix) / sum(filt_matrix),
      counts    = Matrix::rowSums(filt_matrix)
    )
    sc_r      <- SoupX::setSoupProfile(sc_r, soup_prof)
    sc_r      <- SoupX::autoEstCont(sc_r, forceAccept = TRUE, doPlot = FALSE, verbose = FALSE)
    soupx_reduced[[dataset]] <- suppressWarnings(SoupX::adjustCounts(sc_r, verbose = FALSE))
    t_sxr <- tictoc::toc(quiet = TRUE)
    log_mem(paste("After SoupX_filt", names(counts)[dataset]), log_file, bef$gc, bef$rss)
    running_times <- rbind(running_times, data.frame(
      tool_and_dataset = paste0("Dataset", dataset, "_SoupX_reduced"),
      time_seconds     = t_sxr$toc - t_sxr$tic
    ))

    # --- scCDC ---
    message("\tRunning scCDC.")
    bef <- log_mem(paste("Before scCDC", names(counts)[dataset]), log_file)
    tictoc::tic()
    if (!is.null(dataset_name) && dataset_name == "smart-seq2") {
      out <- capture.output(GCGs <- suppressMessages(suppressWarnings(
        scCDC::ContaminationDetection(seu_list[[dataset]], min.cell = 10)
      )))
      out <- capture.output(seuratobj_corrected <- suppressMessages(suppressWarnings(
        scCDC::ContaminationCorrection(seu_list[[dataset]], rownames(GCGs), min.cell = 10)
      )))
    } else {
      out <- capture.output(GCGs <- suppressMessages(suppressWarnings(
        scCDC::ContaminationDetection(seu_list[[dataset]])
      )))
      out <- capture.output(seuratobj_corrected <- suppressMessages(suppressWarnings(
        scCDC::ContaminationCorrection(seu_list[[dataset]], rownames(GCGs))
      )))
    }
    Seurat::DefaultAssay(seuratobj_corrected) <- "Corrected"
    cdc[[dataset]] <- Seurat::GetAssayData(seuratobj_corrected, "Corrected", layer = "counts")
    t_cd <- tictoc::toc(quiet = TRUE)
    log_mem(paste("After scCDC", names(counts)[dataset]), log_file, bef$gc, bef$rss)
    running_times <- rbind(running_times, data.frame(
      tool_and_dataset = paste0("Dataset", dataset, "_scCDC"),
      time_seconds     = t_cd$toc - t_cd$tic
    ))

    # --- DecontX full ---
    message("\tRunning DecontX.")
    bef <- log_mem(paste("Before DecontX", names(counts)[dataset]), log_file)
    tictoc::tic()
    out <- capture.output(decontx_full[[dataset]] <- celda::decontX(
      filt_matrix, background = raw_matrix, verbose = FALSE
    ))
    decontx_full[[dataset]] <- decontx_full[[dataset]]$decontXcounts
    t_dx <- tictoc::toc(quiet = TRUE)
    log_mem(paste("After DecontX", names(counts)[dataset]), log_file, bef$gc, bef$rss)
    running_times <- rbind(running_times, data.frame(
      tool_and_dataset = paste0("Dataset", dataset, "_DecontX_full"),
      time_seconds     = t_dx$toc - t_dx$tic
    ))

    # --- DecontX reduced ---
    message("\tRunning DecontX (filtered only).")
    bef <- log_mem(paste("Before DecontX_filt", names(counts)[dataset]), log_file)
    tictoc::tic()
    out <- capture.output(decontx_reduced[[dataset]] <- celda::decontX(
      filt_matrix, verbose = FALSE
    ))
    decontx_reduced[[dataset]] <- decontx_reduced[[dataset]]$decontXcounts
    t_dxr <- tictoc::toc(quiet = TRUE)
    log_mem(paste("After DecontX_filt", names(counts)[dataset]), log_file, bef$gc, bef$rss)
    running_times <- rbind(running_times, data.frame(
      tool_and_dataset = paste0("Dataset", dataset, "_DecontX_reduced"),
      time_seconds     = t_dxr$toc - t_dxr$tic
    ))
  }

  # Save running times
  dir.create(here::here("results/running_times"), recursive = TRUE, showWarnings = FALSE)
  write.table(
    running_times,
    file      = here::here(paste0("results/running_times/", dataset_name, "_running_times.txt")),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )

  # ------------------------------------------------------------------
  # 5. Add corrected counts as assays
  # ------------------------------------------------------------------
  message("\nAdding assays.")
  seu <- add_assay(seu, fastcar,         "FastCAR")
  seu <- add_assay(seu, cdc,             "scCDC")
  seu <- add_assay(seu, soupx_full,      "SoupXfull")
  seu <- add_assay(seu, soupx_reduced,   "SoupXreduced")
  seu <- add_assay(seu, decontx_full,    "DecontXfull")
  seu <- add_assay(seu, decontx_reduced, "DecontXreduced")

  return(seu)
}


## Python-based methods ---------------------------------------------------

#' Read a corrected count matrix from an H5 or H5AD file
#'
#' Parses sparse matrix indices via rhdf5. Supports CellBender (H5),
#' scAR (H5AD), and CellClear (H5AD) output formats.
#'
#' @param path Character. Path to the H5 or H5AD file.
#' @param type Character. One of "CellBender", "scAR", or "CellClear"
#'   (default "CellBender").
#' @param name Character. Sample name appended to barcodes for uniqueness.
#'
#' @return A sparse dgCMatrix (genes x barcodes).
read_H5AD <- function(path, type = "CellBender", name) {
  root <- rhdf5::H5Fopen(path)

  if (type %in% c("scAR", "CellClear")) {
    i <- as.vector(unlist(rhdf5::h5read(root, "/X/indices")))
    p <- as.vector(unlist(rhdf5::h5read(root, "/X/indptr")))
    x <- as.vector(unlist(rhdf5::h5read(root, "/X/data")))
    o <- as.vector(rhdf5::h5read(root, "/obs")$"_index")
    v <- as.vector(rhdf5::h5read(root, "/var")$"_index")
  } else if (type == "CellBender") {
    i <- as.vector(unlist(rhdf5::h5read(root, "/matrix/indices")))
    p <- as.vector(unlist(rhdf5::h5read(root, "/matrix/indptr")))
    x <- as.vector(unlist(rhdf5::h5read(root, "/matrix/data")))
    o <- as.vector(rhdf5::h5read(root, "/matrix/barcodes"))
    v <- as.vector(rhdf5::h5read(root, "/matrix/features/name"))
  }

  rhdf5::H5Fclose(root)
  dims   <- c(length(v), length(o))
  counts <- Matrix::sparseMatrix(i = i, p = p, x = x, index1 = FALSE, dims = dims)
  rownames(counts) <- v
  colnames(counts) <- paste(o, name, sep = "_")
  return(counts)
}


#' Add a list of corrected count matrices as a new Seurat assay
#'
#' Aligns each matrix to the features and barcodes in `object`, padding
#' missing genes with zeros, then inserts the combined matrix as a new
#' assay and adds species fraction metadata columns.
#'
#' @param object A Seurat object.
#' @param mat Named list of corrected count matrices (one per sample).
#' @param name Character. Name for the new assay (e.g. "CellBender").
#'
#' @return The Seurat object with the new assay added.
add_assay <- function(object, mat, name) {
  for (dataset in seq_along(mat)) {
    rownames(mat[[dataset]]) <- gsub("_", "-", rownames(mat[[dataset]]))
    mat[[dataset]] <- mat[[dataset]][
      rownames(mat[[dataset]]) %in% rownames(object[["RNA"]]$counts), 
    ]
    mat[[dataset]] <- mat[[dataset]][!duplicated(rownames(mat[[dataset]])), ]

    if (nrow(object) > nrow(mat[[dataset]])) {
      missing_genes  <- rownames(object)[!(rownames(object) %in% rownames(mat[[dataset]]))]
      missing_counts <- matrix(0, ncol = ncol(mat[[dataset]]),
                                nrow = length(missing_genes),
                                dimnames = list(missing_genes, colnames(mat[[dataset]])))
      missing_counts <- as(missing_counts, "dgCMatrix")
      mat[[dataset]] <- rbind(missing_counts, mat[[dataset]])
    }
    mat[[dataset]] <- mat[[dataset]][
      match(rownames(object[["RNA"]]$counts), rownames(mat[[dataset]])), 
    ]
  }

  mat <- do.call("cbind", mat)

  barcodes.keep <- intersect(object$barcode, colnames(mat))
  object        <- subset(object, cells = barcodes.keep)
  mat           <- mat[, match(object$barcode, colnames(mat))]
  colnames(mat) <- colnames(object)

  object[[name]] <- Seurat::CreateAssayObject(mat, min.cells = 0, min.features = 0)
  object[[paste0(name, ".percent.mm10")]]   <- Seurat::PercentageFeatureSet(
    object, pattern = "^mm10",   assay = name
  )[, 1]
  object[[paste0(name, ".percent.GRCh38")]] <- Seurat::PercentageFeatureSet(
    object, pattern = "^GRCh38", assay = name
  )[, 1]
  return(object)
}


# Post-decontamination processing -----------------------------------------

#' Fix triple-dash artefacts in feature names across all assays
#'
#' Replaces "---" with "-" in rownames of all assays and their
#' data / scale.data layers.
#'
#' @param seu A Seurat object.
#'
#' @return The Seurat object with corrected feature names.
correct_rownames <- function(seu) {
  var_feats <- gsub("---", "-", Seurat::VariableFeatures(seu, assay = "RNA"))
  for (assay in names(seu@assays)) {
    cnts <- Seurat::GetAssayData(seu, assay = assay, layer = "counts")
    rownames(cnts) <- gsub("---", "-", rownames(cnts))
    dat  <- Seurat::GetAssayData(seu, assay = assay, layer = "data")
    rownames(dat)  <- gsub("---", "-", rownames(dat))
    new_assay              <- Seurat::CreateAssayObject(cnts)
    new_assay@data         <- dat
    scdat <- Seurat::GetAssayData(seu, assay = assay, layer = "scale.data")
    if (length(scdat) > 0) {
      rownames(scdat)         <- gsub("---", "-", rownames(scdat))
      new_assay@scale.data    <- scdat
    }
    seu[[assay]] <- new_assay
    Seurat::VariableFeatures(seu, assay = assay) <- var_feats
  }
  return(seu)
}


#' Convert Ensembl IDs to HGNC gene symbols (local, offline)
#'
#' Uses the `org.Hs.eg.db` Bioconductor annotation package to map Ensembl
#' IDs to gene symbols. Version numbers (e.g. ".1") are stripped before
#' mapping. Falls back to the original versioned Ensembl ID where no
#' symbol is found, and deduplicates with \code{make.unique()}.
#'
#' Works on both Seurat objects and plain sparse matrices.
#'
#' @param mat A Seurat object or sparse matrix with Ensembl IDs as rownames.
#'
#' @return The input object with rownames replaced by gene symbols.
convert_ensembl_ids_to_gene_names <- function(mat) {
  ensembl_ids_with_version <- rownames(mat)
  ensembl_ids              <- sub("\\..*", "", ensembl_ids_with_version)

  gene_names <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys      = ensembl_ids,
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )

  new_rownames <- gene_names[ensembl_ids]
  new_rownames[is.na(new_rownames)] <- ensembl_ids_with_version[is.na(new_rownames)]
  rownames(mat) <- make.unique(unname(new_rownames))
  return(mat)
}


#' Remove cells with high mitochondrial content
#'
#' Calculates the fraction of mitochondrial UMIs per cell and removes cells
#' above the specified threshold. The fraction is stored as `fraction.mt`
#' in metadata.
#'
#' @param seu A Seurat object.
#' @param mito_pattern Character vector. One or more regex patterns matching
#'   mitochondrial gene names (e.g. \code{c("^mm10-mt-", "^GRCh38-MT-")}).
#' @param mt_threshold Numeric. Maximum mitochondrial fraction (0–1).
#'
#' @return The filtered Seurat object.
remove_mito_genes <- function(seu, mito_pattern, mt_threshold) {
  mito_genes <- unlist(lapply(mito_pattern, function(p) {
    rownames(seu)[grepl(paste0("^", p), rownames(seu))]
  }))
  seu[["fraction.mt"]] <- Matrix::colSums(
    seu@assays$RNA$counts[rownames(seu@assays$RNA$counts) %in% mito_genes, ]
  ) / Matrix::colSums(seu@assays$RNA$counts)
  subset(seu, subset = fraction.mt <= mt_threshold)
}


#' Assign human / mouse species-of-origin to cells
#'
#' Uses the mm10 UMI fraction to assign cells to human (≤ 25% mm10) or
#' mouse (≥ 75% mm10). Cells within the 25–75% grey zone are removed.
#'
#' @param seu A Seurat object with `input.percent.mm10` metadata.
#'
#' @return The Seurat object with an `origin` metadata column.
determine_species_origin <- function(seu) {
  seu <- subset(seu, cells = names(which(abs(seu$input.percent.mm10 - 50) >= 25)))
  seu$origin <- "NA"
  seu@meta.data[seu@meta.data$input.percent.mm10 <= 25, "origin"] <- "human"
  seu@meta.data[seu@meta.data$input.percent.mm10 >= 75, "origin"] <- "mouse"
  subset(seu, origin != "NA")
}


# Ambient RNA analysis ----------------------------------------------------

#' Extract raw counts for cells of a given species or cell type
#'
#' @param seu A Seurat object with `origin` or `celltype` metadata.
#' @param origin Character. One of "human", "mouse", "MCF12A", "HCC1500",
#'   "HS578T".
#'
#' @return A sparse count matrix for the selected cells.
extract_raw_counts <- function(seu, origin = c("human", "mouse",
                                                "MCF12A", "HCC1500", "HS578T")) {
  origin <- match.arg(origin)
  if (origin %in% c("human", "mouse")) {
    seu@assays$RNA$counts[, seu$origin   == origin, drop = FALSE]
  } else {
    seu@assays$RNA$counts[, seu$celltype == origin, drop = FALSE]
  }
}


#' Subset a count matrix to ambient genes (by species prefix)
#'
#' @param raw_counts Sparse count matrix.
#' @param ambient_pattern Character. Regex matching ambient gene names.
#'
#' @return Subsetted count matrix.
extract_ambient_counts <- function(raw_counts,
                                   ambient_pattern = c("^mm10-", "^GRCh38-")) {
  ambient_pattern <- match.arg(ambient_pattern)
  raw_counts[grepl(ambient_pattern, rownames(raw_counts)), , drop = FALSE]
}


#' Subset a count matrix to endogenous genes (by species prefix)
#'
#' @param raw_counts Sparse count matrix.
#' @param endogenous_pattern Character. Regex matching endogenous gene names.
#'
#' @return Subsetted count matrix.
extract_endogenous_counts <- function(raw_counts,
                                      endogenous_pattern = c("^GRCh38-", "^mm10-")) {
  endogenous_pattern <- match.arg(endogenous_pattern)
  raw_counts[grepl(endogenous_pattern, rownames(raw_counts)), , drop = FALSE]
}


#' Calculate per-cell ambient RNA load percentages
#'
#' Computes the ambient fraction (as a percentage) per cell. When `samples`
#' and `seu` are provided, results are grouped by sample.
#'
#' @param ambient_list Named list of ambient count matrices.
#' @param raw_list Named list of raw count matrices (same names as ambient_list).
#' @param samples Optional character vector of sample names.
#' @param seu Optional Seurat object (required when `samples` is not NULL).
#'
#' @return A data frame with columns: group, barcode, ambient_umis,
#'   total_umis, ambient_percent.
calc_ambient_load_percentages <- function(ambient_list, raw_list,
                                           samples = NULL, seu = NULL) {
  if (is.null(names(ambient_list)) || any(names(ambient_list) == ""))
    stop("All elements of ambient_list must be named.")
  if (is.null(names(raw_list)) || any(names(raw_list) == ""))
    stop("All elements of raw_list must be named.")
  if (!setequal(names(ambient_list), names(raw_list)))
    stop("ambient_list and raw_list must have the same names.")

  raw_list <- raw_list[names(ambient_list)]

  if (is.null(samples)) {
    out <- lapply(names(ambient_list), function(name) {
      barcodes        <- colnames(raw_list[[name]])
      ambient_umis    <- Matrix::colSums(ambient_list[[name]], na.rm = TRUE)
      total_umis      <- Matrix::colSums(raw_list[[name]],    na.rm = TRUE)
      data.frame(
        group           = "S1",
        barcode         = barcodes,
        ambient_umis    = ambient_umis,
        total_umis      = total_umis,
        ambient_percent = (ambient_umis / total_umis) * 100,
        stringsAsFactors = FALSE
      )
    })
    df_out <- do.call(rbind, out)
    return(df_out[!is.na(df_out$ambient_percent), ])
  }

  out <- lapply(seq_along(samples), function(i) {
    sample_barcodes <- rownames(seu@meta.data[seu$orig.ident == samples[i], ])
    rows <- lapply(names(ambient_list), function(name) {
      keep         <- colnames(raw_list[[name]]) %in% sample_barcodes
      if (!any(keep)) return(NULL)
      barcodes     <- colnames(raw_list[[name]])[keep]
      ambient_umis <- Matrix::colSums(ambient_list[[name]][, keep, drop = FALSE], na.rm = TRUE)
      total_umis   <- Matrix::colSums(raw_list[[name]][,    keep, drop = FALSE], na.rm = TRUE)
      data.frame(
        group           = paste0("S", i),
        barcode         = barcodes,
        ambient_umis    = ambient_umis,
        total_umis      = total_umis,
        ambient_percent = (ambient_umis / total_umis) * 100,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, Filter(Negate(is.null), rows))
  })

  df_out <- do.call(rbind, out)
  df_out[!is.na(df_out$ambient_percent), ]
}


#' Calculate average ambient load percentage per sample group
#'
#' @param ambient_load_percentages Data frame as returned by
#'   \code{calc_ambient_load_percentages()}.
#'
#' @return A data frame with columns `group` and `value` (mean ambient %).
calc_avg_ambient_load_percentage <- function(ambient_load_percentages) {
  ambient_load_percentages |>
    dplyr::group_by(group) |>
    dplyr::summarise(value = signif(mean(ambient_percent, na.rm = TRUE), 3),
                     .groups = "drop") |>
    dplyr::arrange(as.numeric(sub("S", "", group)))
}


#' Calculate per-cell ambient fraction from raw matrices
#'
#' @param raw_ambient_1 Ambient count matrix (first species/direction).
#' @param raw_1 Raw count matrix matching raw_ambient_1.
#' @param raw_ambient_2 Ambient count matrix (second species/direction).
#'   NULL for single-species.
#' @param raw_2 Raw count matrix matching raw_ambient_2. NULL for single-species.
#'
#' @return Named numeric vector of per-cell ambient fractions.
calc_ambient_load_fractions <- function(raw_ambient_1 = NULL, raw_1 = NULL,
                                         raw_ambient_2 = NULL, raw_2 = NULL) {
  f1 <- if (!is.null(raw_ambient_1) && !is.null(raw_1))
    Matrix::colSums(raw_ambient_1, na.rm = TRUE) / Matrix::colSums(raw_1, na.rm = TRUE)
  f2 <- if (!is.null(raw_ambient_2) && !is.null(raw_2))
    Matrix::colSums(raw_ambient_2, na.rm = TRUE) / Matrix::colSums(raw_2, na.rm = TRUE)
  c(f1, f2)
}


#' Calculate per-cell and per-gene exon fractions
#'
#' Compares Gene (exon-only) and GeneFull_Ex50pAS (full gene body) STARsolo
#' output folders to derive the fraction of UMIs from exonic regions.
#'
#' @param path Character. Path to the STARsolo Solo.out directory.
#' @param seu A Seurat object used to match barcodes.
#' @param name Character. Sample name used in barcode construction.
#' @param paste1 Logical. Format barcodes as "<bc>-1_<name>" (TRUE, default)
#'   or "<bc>_<name>" (FALSE).
#'
#' @return A named list with elements \code{cells} and \code{features}.
calc_exon_fraction <- function(path, seu, name, paste1 = TRUE) {
  total <- load_folder(file.path(path, "GeneFull_Ex50pAS", "raw"))
  exon  <- load_folder(file.path(path, "Gene", "raw"))

  total <- total[rownames(total) %in% rownames(exon), ]
  exon  <- exon[rownames(total), ]

  if (paste1) {
    colnames(total) <- paste0(colnames(total), "-1_", name)
    colnames(exon)  <- paste0(colnames(exon),  "-1_", name)
  } else {
    colnames(total) <- paste0(colnames(total), "_", name)
    colnames(exon)  <- paste0(colnames(exon),  "_", name)
  }

  total <- total[, colnames(total) %in% colnames(seu), drop = FALSE]
  exon  <- exon[,  colnames(exon)  %in% colnames(seu), drop = FALSE]

  features_exon <- rowSums(exon) / rowSums(total)
  # Normalise gene name formatting
  names(features_exon) <- gsub("---", "-", gsub("_", "-", names(features_exon)))

  list(
    cells    = colSums(exon) / colSums(total),
    features = features_exon
  )
}


#' Match ambient load with exon fractions (barnyard / high-complexity datasets)
#'
#' Builds a data frame of per-gene exon fractions and log-ambient expression
#' for both species.
#'
#' @param exon_fraction_per_feature Named numeric vector of per-gene exon
#'   fractions (used for both species when a single vector is provided).
#' @param raw_human_ambient Sparse matrix of ambient counts in human cells.
#' @param raw_mouse_ambient Sparse matrix of ambient counts in mouse cells.
#'
#' @return A data frame with columns exon_fraction, log_ambient_expression,
#'   species, features.
match_ambient_load_with_exon_fractions <- function(exon_fraction_per_feature = NULL,
                                                    raw_human_ambient,
                                                    raw_mouse_ambient) {
  ef_human <- exon_fraction_per_feature
  ef_mouse <- exon_fraction_per_feature

  make_df <- function(ef, raw_ambient, species_label) {
    raw_ambient <- raw_ambient[rownames(raw_ambient) %in% names(ef), , drop = FALSE]
    ef          <- ef[names(ef) %in% rownames(raw_ambient)]
    raw_ambient <- raw_ambient[match(names(ef), rownames(raw_ambient)), , drop = FALSE]
    data.frame(
      exon_fraction          = ef,
      log_ambient_expression = log1p(rowSums(raw_ambient)),
      species                = species_label,
      features               = rownames(raw_ambient)
    )
  }

  rbind(
    make_df(ef_human, raw_human_ambient, "human"),
    make_df(ef_mouse, raw_mouse_ambient, "mouse")
  )
}


#' Match ambient load with exon fractions (medium complexity, cross-species)
#'
#' Species-aware version for the medium complexity dataset, where human and
#' mouse exon fractions come from different sample sets.
#'
#' @param exon_fraction_per_feature_1 Per-gene exon fractions for human samples.
#' @param exon_fraction_per_feature_2 Per-gene exon fractions for mouse samples.
#' @param raw_human_ambient Ambient matrix from human samples.
#' @param raw_mouse_ambient Ambient matrix from mouse samples.
#'
#' @return A data frame with columns exon_fraction, log_ambient_expression,
#'   species.
match_ambient_load_with_exon_fractions_mc <- function(exon_fraction_per_feature_1,
                                                       exon_fraction_per_feature_2,
                                                       raw_human_ambient,
                                                       raw_mouse_ambient) {
  make_df <- function(ef, prefix, raw_ambient, species_label) {
    ef          <- ef[grepl(paste0("^", prefix), names(ef))]
    ef          <- ef[names(ef) %in% rownames(raw_ambient)]
    raw_ambient <- raw_ambient[rownames(raw_ambient) %in% names(ef), , drop = FALSE]
    raw_ambient <- raw_ambient[match(names(ef), rownames(raw_ambient)), , drop = FALSE]
    data.frame(
      exon_fraction          = ef,
      log_ambient_expression = log1p(rowSums(raw_ambient)),
      species                = species_label
    )
  }

  rbind(
    make_df(exon_fraction_per_feature_1, "GRCh38", raw_human_ambient, "human"),
    make_df(exon_fraction_per_feature_2, "mm10",   raw_mouse_ambient, "mouse")
  )
}


#' Get gene biotype information from a STARsolo geneInfo.tab file
#'
#' @param path Character. Path to the geneInfo.tab file.
#' @param seu A Seurat object used to filter to present genes.
#'
#' @return A data frame with columns ENSEMBL_id, gene_name, gene_biotype.
get_gene_biotypes <- function(path, seu) {
  gene_info <- read.delim(path, skip = 1, header = FALSE)
  colnames(gene_info) <- c("ENSEMBL_id", "gene_name", "gene_biotype")
  gene_info$gene_name <- gsub("_", "-", gene_info$gene_name)
  gene_info <- gene_info[!duplicated(gene_info$gene_name), ]
  gene_info[gene_info$gene_name %in% rownames(seu), ]
}


#' Match ambient load fractions with a gene-level covariate
#'
#' @param ambient_load_fractions Named numeric vector of per-cell ambient
#'   fractions, as returned by \code{calc_ambient_load_fractions()}.
#' @param gene_biotype Named numeric vector of a per-cell covariate (e.g.
#'   exon fraction, PC gene fraction) with cell names as names.
#'
#' @return The covariate vector reordered to match ambient_load_fractions.
match_ambient_load_with_gene_biotypes <- function(ambient_load_fractions, gene_biotype) {
  gene_biotype[match(names(ambient_load_fractions), names(gene_biotype))]
}


#' Calculate per-cell protein-coding gene fraction
#'
#' @param seu A Seurat object.
#' @param gene_info Data frame from \code{get_gene_biotypes()}.
#'
#' @return Named numeric vector of per-cell protein-coding UMI fractions.
calc_fraction_pc_genes_per_cell <- function(seu, gene_info) {
  pc_genes <- gene_info[gene_info$gene_biotype == "protein_coding", "gene_name"]
  Matrix::colSums(
    seu@assays$RNA$counts[rownames(seu@assays$RNA$counts) %in% pc_genes, ]
  ) / Matrix::colSums(seu@assays$RNA$counts)
}


#' Calculate per-cell ribosomal gene fraction
#'
#' @param seu A Seurat object.
#'
#' @return Named numeric vector of per-cell ribosomal UMI fractions.
calc_fraction_ribo_genes_per_cell <- function(seu) {
  Matrix::colSums(
    seu@assays$RNA$counts[grepl("Rps|Rpl|RPS|RPL", rownames(seu)), ]
  ) / Matrix::colSums(seu@assays$RNA$counts)
}


#' Calculate variance in ambient load explained by cellular parameters
#'
#' Used for the low, medium, and high complexity datasets.
#'
#' @param seu A Seurat object.
#' @param species Character. "human" or "mouse".
#' @param ambient_fraction Named numeric vector of per-cell ambient fractions.
#' @param ambient_exon_fraction_per_cell Per-cell exon fractions.
#' @param ambient_pc_fraction_per_cell Per-cell protein-coding fractions.
#' @param ambient_ribo_fraction_per_cell Per-cell ribosomal fractions.
#' @param ambient_mt_fraction_per_cell Per-cell mitochondrial fractions.
#' @param raw_counts Raw count matrix.
#' @param endogenous_counts Endogenous count matrix.
#' @param ambient_counts Ambient count matrix.
#'
#' @return A data frame with columns parameters, percent_explained, species.
calc_var_explained_by_cellular_params <- function(seu,
                                                   species = c("human", "mouse"),
                                                   ambient_fraction,
                                                   ambient_exon_fraction_per_cell,
                                                   ambient_pc_fraction_per_cell,
                                                   ambient_ribo_fraction_per_cell,
                                                   ambient_mt_fraction_per_cell,
                                                   raw_counts,
                                                   endogenous_counts,
                                                   ambient_counts) {
  species <- match.arg(species)
  sp_cells <- rownames(seu@meta.data[seu$origin == species, ])

  af  <- ambient_fraction[names(ambient_fraction) %in% sp_cells]
  ef  <- ambient_exon_fraction_per_cell[names(ambient_exon_fraction_per_cell) %in% sp_cells]
  pcf <- ambient_pc_fraction_per_cell[names(ambient_pc_fraction_per_cell) %in% sp_cells]
  rbf <- ambient_ribo_fraction_per_cell[names(ambient_ribo_fraction_per_cell) %in% sp_cells]
  mtf <- ambient_mt_fraction_per_cell[names(ambient_mt_fraction_per_cell) %in% sp_cells]

  cor2 <- function(x, y) stats::cor(x, y, use = "pairwise.complete.obs")^2

  df <- data.frame(
    parameters = c("UMIs", "Feat", "Exons", "PC", "ribo", "mito"),
    percent_explained = 100 * c(
      cor2(log1p(Matrix::colSums(endogenous_counts)), Matrix::colSums(ambient_counts) / Matrix::colSums(raw_counts)),
      cor2(log1p(Matrix::colSums(endogenous_counts > 0)), Matrix::colSums(ambient_counts) / Matrix::colSums(raw_counts)),
      cor2(af, ef),
      cor2(af, pcf),
      cor2(af, rbf),
      cor2(af, mtf)
    )
  )
  df$species <- species
  df
}


#' Calculate variance in ambient load explained by cellular parameters (genotype)
#'
#' Variant for the genotype dataset which lacks a species metadata column.
#'
#' @inheritParams calc_var_explained_by_cellular_params
#'
#' @return A data frame with columns parameters and percent_explained.
calc_var_explained_by_cellular_params_geno <- function(seu,
                                                        ambient_fraction,
                                                        ambient_exon_fraction_per_cell,
                                                        ambient_pc_fraction_per_cell,
                                                        ambient_ribo_fraction_per_cell,
                                                        ambient_mt_fraction_per_cell,
                                                        raw_counts,
                                                        endogenous_counts,
                                                        ambient_counts) {
  cor2 <- function(x, y) stats::cor(x, y, use = "pairwise.complete.obs")^2

  df <- data.frame(
    parameters = c("UMIs", "Feat", "Exons", "PC", "ribo", "mito"),
    percent_explained = 100 * c(
      cor2(log1p(Matrix::colSums(endogenous_counts, na.rm = TRUE)),
           Matrix::colSums(ambient_counts, na.rm = TRUE) / Matrix::colSums(raw_counts, na.rm = TRUE)),
      cor2(log1p(Matrix::colSums(endogenous_counts > 0, na.rm = TRUE)),
           Matrix::colSums(ambient_counts, na.rm = TRUE) / Matrix::colSums(raw_counts, na.rm = TRUE)),
      cor2(ambient_fraction, ambient_exon_fraction_per_cell),
      cor2(ambient_fraction, ambient_pc_fraction_per_cell),
      cor2(ambient_fraction, ambient_ribo_fraction_per_cell),
      cor2(ambient_fraction, ambient_mt_fraction_per_cell)
    )
  )
  df
}


#' Estimate ambient background expression from empty droplets
#'
#' @param path Character. Path to a count matrix folder or H5 file.
#' @param object A Seurat object used to match the gene set.
#' @param threshold Integer. Max UMI count for a droplet to be considered
#'   empty (default 50).
#'
#' @return Named numeric vector of summed ambient UMI counts per gene.
calc_background_profile <- function(path, object, threshold = 50) {
  if (grepl("\\.h5$", path, ignore.case = TRUE)) {
    counts <- Seurat::Read10X_h5(path)
  } else {
    counts <- load_folder(path)
  }

  background_barcodes <- names(which(Matrix::colSums(counts) <= threshold))
  background          <- counts[, colnames(counts) %in% background_barcodes, drop = FALSE]

  rownames(background) <- gsub("---", "-", gsub("_", "-", rownames(background)))
  background <- background[rownames(background) %in% rownames(object), , drop = FALSE]
  background <- background[!duplicated(rownames(background)), ]

  missing_genes  <- rownames(object)[!(rownames(object) %in% rownames(background))]
  missing_counts <- Matrix::Matrix(0, nrow = length(missing_genes),
                                   ncol = ncol(background), sparse = TRUE,
                                   dimnames = list(missing_genes, colnames(background)))
  background <- rbind(missing_counts, background)
  background <- background[match(rownames(object), rownames(background)), ]

  Matrix::rowSums(background)
}


#' Estimate per-cell ambient RNA fraction from strain SNP data (genotype)
#'
#' For each sample, reads the intermediate SNP calling results, filters to
#' high-confidence CAST-specific SNPs (>= 2 supporting reads), and calculates
#' the ambient fraction per unique gene-in-cell combination as ALT / total.
#'
#' @param sample Character. Sample name matching the directory structure
#'   under data/processed/genotype/genotype_estimation/.
#'
#' @return A data frame with columns Barcode, Symbol, Ambient, Position.
calc_ambient_from_SNPs <- function(sample) {
  sa <- readRDS(here::here(paste0(
    "data/processed/genotype/genotype_estimation/", sample, "/csc_intermediate.RDS"
  )))

  sa <- sa[sa$Strain == "CAST", ]
  sa <- sa[sa$DP >= 2, ]
  sa <- sa[sa$Strain != sa$SNPisALT, ]

  sa$cell <- paste(sa$cell, sample, sep = "_")
  sa <- sa[sa$cell %in% gen_seu$barcode, ]

  keep <- vapply(seq_len(nrow(sa)), function(i) {
    !grepl(sa$Strain[i], sa$SNPisALT[i])
  }, logical(1))
  sa <- sa[keep, ]

  sa$ID <- paste(sa$cell, sa$SYMBOL, sep = ":")

  res <- do.call(rbind, lapply(unique(sa$ID), function(id) {
    rows <- sa[sa$ID == id, ]
    data.frame(
      Barcode  = sub(":.*", "", id),
      Symbol   = sub(".*:", "", id),
      Ambient  = sum(rows$AD) / sum(rows$DP),
      Position = rows$position[1]
    )
  }))
  res
}


#' Quartile coefficient of dispersion
#'
#' @param x Numeric vector.
#' @param na.rm Logical. Remove NA values (default TRUE).
#'
#' @return Scalar QCD value.
calc_qcd <- function(x, na.rm = TRUE) {
  q <- stats::quantile(x, c(0.25, 0.75), na.rm = na.rm, names = FALSE)
  (q[2] - q[1]) / (q[2] + q[1])
}


# Removal evaluation helpers ----------------------------------------------

#' Divide corrected counts into ambient and endogenous compartments per tool
#'
#' Extracts total, ambient, and endogenous count matrices for cells of a
#' given species from each decontamination assay.
#'
#' @param seu A Seurat object containing all decontamination assays.
#' @param extract Character. "human" or "mouse".
#'
#' @return A named list of sparse matrices (input_total, input_ambient,
#'   input_endogenous, cb_total, cb_ambient, ...).
divide_ambient_and_endogenous_per_tool <- function(seu, extract) {
  if (extract == "human") {
    ambient_pat    <- "^mm10-"
    endogenous_pat <- "^GRCh38-"
  } else {
    ambient_pat    <- "^GRCh38-"
    endogenous_pat <- "^mm10-"
  }

  cells <- rownames(seu@meta.data[seu$origin == extract, ])

  tools <- c(input = "RNA", cb = "CellBender", sx = "SoupXfull",
             sxr = "SoupXreduced", dx = "DecontXfull", dxr = "DecontXreduced",
             fc = "FastCAR", ar = "scAR", cd = "scCDC", cc = "CellClear")

  out <- list()
  for (abbr in names(tools)) {
    assay_name <- tools[[abbr]]
    total <- seu@assays[[assay_name]]$counts[, colnames(seu) %in% cells, drop = FALSE]
    out[[paste0(abbr, "_total")]]       <- total
    out[[paste0(abbr, "_ambient")]]     <- total[grepl(ambient_pat,    rownames(total)), ]
    out[[paste0(abbr, "_endogenous")]]  <- total[grepl(endogenous_pat, rownames(total)), ]
  }
  out
}


#' Collect total counts per tool (Smart-seq2 / no-species-split datasets)
#'
#' @param seu A Seurat object.
#'
#' @return Named list of endogenous count matrices per tool.
get_total_counts_per_tool <- function(seu) {
  list(
    input_endogenous = seu@assays$RNA$counts,
    cb_endogenous    = seu@assays$CellBender$counts,
    sx_endogenous    = seu@assays$SoupXfull$counts,
    sxr_endogenous   = seu@assays$SoupXreduced$counts,
    dx_endogenous    = seu@assays$DecontXfull$counts,
    dxr_endogenous   = seu@assays$DecontXreduced$counts,
    fc_endogenous    = seu@assays$FastCAR$counts,
    ar_endogenous    = seu@assays$scAR$counts,
    cd_endogenous    = seu@assays$scCDC$counts,
    cc_endogenous    = seu@assays$CellClear$counts
  )
}


#' Build a sparse matrix subset from a feature-barcode index
#'
#' @param symbols Character vector of gene names (rows to extract).
#' @param barcodes Character vector of cell barcodes (columns to extract).
#' @param rownames_full Full rowname set of the output matrix.
#' @param colnames_full Full colname set of the output matrix.
#' @param source_mat Source sparse matrix.
#'
#' @return A sparse matrix with the same dimensions as the full index.
make_sparse_subset <- function(symbols, barcodes, rownames_full, colnames_full, source_mat) {
  row_idx <- match(symbols,  rownames_full)
  col_idx <- match(barcodes, colnames_full)
  keep    <- !is.na(row_idx) & !is.na(col_idx)
  row_idx <- row_idx[keep]
  col_idx <- col_idx[keep]
  vals    <- source_mat[cbind(row_idx, col_idx)]
  Matrix::sparseMatrix(
    i = row_idx, j = col_idx, x = vals,
    dims     = c(length(rownames_full), length(colnames_full)),
    dimnames = list(rownames_full, colnames_full)
  )
}


#' Build a sparse ambient-removed matrix (raw - corrected) for a subset
#'
#' @param symbols Character vector of gene names.
#' @param barcodes Character vector of cell barcodes.
#' @param rn Full rowname set.
#' @param cn Full colname set.
#' @param raw_mat Raw count matrix.
#' @param corr_mat Corrected count matrix.
#'
#' @return A sparse matrix of (raw - corrected) counts for the given subset.
make_ambient_subset <- function(symbols, barcodes, rn, cn, raw_mat, corr_mat) {
  row_idx <- match(symbols,  rn)
  col_idx <- match(barcodes, cn)
  keep    <- !is.na(row_idx) & !is.na(col_idx)
  row_idx <- row_idx[keep]
  col_idx <- col_idx[keep]
  vals    <- raw_mat[cbind(row_idx, col_idx)] - corr_mat[cbind(row_idx, col_idx)]
  Matrix::sparseMatrix(
    i = row_idx, j = col_idx, x = vals,
    dims     = c(length(rn), length(cn)),
    dimnames = list(rn, cn)
  )
}


#' Evaluate endogenous signal retention vs total endogenous UMIs per feature
#'
#' @param metadata Data frame with orig.ident column.
#' @param dataset Character. Dataset label for output.
#' @param corrected_counts Named list of corrected count matrices from
#'   \code{divide_ambient_and_endogenous_per_tool()} or
#'   \code{get_total_counts_per_tool()}.
#' @param ngroups Integer. Number of gene UMI groups (default 1).
#'
#' @return Data frame with columns method, group, error, fraction_changed,
#'   dataset.
calc_endogenous_removal_vs_total_endogenous_UMIs_per_feature <- function(
    metadata, dataset, corrected_counts, ngroups = 1) {

  input_endogenous <- corrected_counts$input_endogenous
  tool_names <- c("sx", "sxr", "fc", "dx", "dxr", "cd", "cb", "ar", "cc")

  endogenous_results <- data.frame()

  for (samp in unique(metadata$orig.ident)) {
    samp_cells <- rownames(metadata[metadata$orig.ident == samp, ])
    tmp_rs <- Matrix::rowSums(
      input_endogenous[, colnames(input_endogenous) %in% samp_cells, drop = FALSE]
    )
    tmp_rs <- sort(tmp_rs[tmp_rs > 0])
    npergroup <- floor(length(tmp_rs) / ngroups)

    groups <- setNames(
      rep(seq_len(ngroups), c(rep(npergroup, ngroups - 1), length(tmp_rs) - npergroup * (ngroups - 1))),
      names(tmp_rs)
    )

    tmp_input <- input_endogenous[, colnames(input_endogenous) %in% samp_cells, drop = FALSE]

    error_list <- lapply(tool_names, function(tool) {
      mat <- corrected_counts[[paste0(tool, "_endogenous")]]
      mat[, colnames(mat) %in% samp_cells, drop = FALSE]
    })
    names(error_list) <- tool_names

    for (grp in seq_len(ngroups)) {
      genes_in_group <- names(groups)[groups == grp]
      for (method in tool_names) {
        tmp_err   <- error_list[[method]][rownames(error_list[[method]]) %in% genes_in_group, ]
        tmp_inp   <- tmp_input[rownames(tmp_input) %in% genes_in_group, ]
        tmp_SE    <- pmax(0, 1 - abs(rowSums(tmp_err) - rowSums(tmp_inp)) / rowSums(tmp_inp))
        endogenous_results <- rbind(endogenous_results, data.frame(
          method           = method,
          group            = grp,
          error            = mean(tmp_SE),
          fraction_changed = sum(tmp_SE) / length(genes_in_group),
          dataset          = dataset
        ))
      }
    }
  }

  endogenous_results$group <- as.factor(endogenous_results$group)
  endogenous_results
}


#' Evaluate ambient RNA removal vs total ambient UMIs per feature
#'
#' @inheritParams calc_endogenous_removal_vs_total_endogenous_UMIs_per_feature
#'
#' @return Data frame with columns method, group, error, fraction_changed,
#'   dataset.
calc_ambient_removal_vs_total_ambient_UMIs_per_feature <- function(
    metadata, dataset, corrected_counts, ngroups = 1) {

  input_ambient <- corrected_counts$input_ambient
  tool_names    <- c("sx", "sxr", "fc", "dx", "dxr", "cd", "cb", "ar", "cc")

  ambient_results <- data.frame()

  for (samp in unique(metadata$orig.ident)) {
    samp_cells <- rownames(metadata[metadata$orig.ident == samp, ])
    tmp_rs <- Matrix::rowSums(
      input_ambient[, colnames(input_ambient) %in% samp_cells, drop = FALSE]
    )
    tmp_rs <- sort(tmp_rs[tmp_rs > 0])
    npergroup <- floor(length(tmp_rs) / ngroups)

    groups <- setNames(
      rep(seq_len(ngroups), c(rep(npergroup, ngroups - 1), length(tmp_rs) - npergroup * (ngroups - 1))),
      names(tmp_rs)
    )

    tmp_input <- input_ambient[, colnames(input_ambient) %in% samp_cells, drop = FALSE]

    error_list <- lapply(tool_names, function(tool) {
      mat <- corrected_counts[[paste0(tool, "_ambient")]]
      mat[, colnames(mat) %in% samp_cells, drop = FALSE]
    })
    names(error_list) <- tool_names

    for (grp in seq_len(ngroups)) {
      genes_in_group <- names(groups)[groups == grp]
      for (method in tool_names) {
        tmp_err <- error_list[[method]][rownames(error_list[[method]]) %in% genes_in_group, ]
        tmp_inp <- tmp_input[rownames(tmp_input) %in% genes_in_group, ]
        tmp_SE  <- pmin(1, abs(rowSums(tmp_inp) - rowSums(tmp_err)) / rowSums(tmp_inp))
        ambient_results <- rbind(ambient_results, data.frame(
          method           = method,
          group            = grp,
          error            = mean(tmp_SE),
          fraction_changed = sum(tmp_SE) / length(genes_in_group),
          dataset          = dataset
        ))
      }
    }
  }

  ambient_results$group <- as.factor(ambient_results$group)
  ambient_results
}


#' Run ambient or endogenous removal analysis across all tools
#'
#' Generates per-tool scatter plots of fraction removed vs log UMIs, saves
#' PDFs, and writes a combined xlsx table.
#'
#' @param seu A Seurat object.
#' @param dataset_name Character. Label used in output filenames.
#' @param compartment Character. "ambient" or "endogenous".
#' @param plot_dir Character. Directory for output PDF plots.
#' @param table_path Character. File path for the output xlsx table.
#'
#' @return A named list with elements plots, data, per_tool_species_dfs.
run_compartment_removal_analysis <- function(seu,
                                              dataset_name,
                                              compartment = c("ambient", "endogenous"),
                                              plot_dir,
                                              table_path) {
  compartment <- match.arg(compartment)

  tools      <- c("CellBender", "CellClear", "FastCAR", "scAR", "scCDC",
                  "SoupX", "SoupX_r", "DecontX", "DecontX_r")
  tools_abbr <- c("cb", "cc", "fc", "ar", "cd", "sx", "sxr", "dx", "dxr")
  tools_map  <- setNames(tools_abbr, tools)
  input_name <- paste0("input_", compartment)

  corrected <- list(
    human = divide_ambient_and_endogenous_per_tool(seu, "human"),
    mouse = divide_ambient_and_endogenous_per_tool(seu, "mouse")
  )

  compartment_dfs <- list()
  for (species in names(corrected)) {
    for (tool in tools) {
      abbr      <- tools_map[[tool]]
      tool_name <- paste0(abbr, "_", compartment)
      compartment_dfs[[paste(tool, species, sep = "_")]] <- data.frame(
        log_UMIs = log1p(rowSums(corrected[[species]][[input_name]])),
        fraction_removed = (rowSums(corrected[[species]][[input_name]]) -
                              rowSums(corrected[[species]][[tool_name]])) /
          rowSums(corrected[[species]][[input_name]]),
        species = species
      )
    }
  }

  get_plot_df <- function(tool) {
    rbind(
      compartment_dfs[[paste(tool, "human", sep = "_")]],
      compartment_dfs[[paste(tool, "mouse", sep = "_")]]
    )
  }

  plots <- setNames(
    lapply(seq_along(tools), function(i) {
      ggplot2::ggplot(get_plot_df(tools[[i]]),
                      ggplot2::aes(x = log_UMIs, y = fraction_removed)) +
        ggrastr::geom_point_rast(colour = custom_colors[[i]], size = 0.5,
                                  stroke = 0, alpha = 0.5, raster.dpi = 1500) +
        plot_theme() +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::ylim(-1, 1)
    }),
    paste0(dataset_name, "_", compartment, "_removal_", tools)
  )

  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  for (plot_name in names(plots)) {
    ggplot2::ggsave(
      filename = file.path(plot_dir, paste0(plot_name, ".pdf")),
      plot     = plots[[plot_name]],
      dpi = 300, width = 35, height = 35, units = "mm"
    )
  }

  final_df <- do.call(rbind, lapply(tools, function(tool) {
    df <- get_plot_df(tool)
    df$decontamination_method <- tool
    df
  }))
  final_df$dataset     <- dataset_name
  final_df$compartment <- compartment
  writexl::write_xlsx(final_df, path = table_path)

  list(plots = plots, data = final_df, per_tool_species_dfs = compartment_dfs)
}


# Figure helpers ----------------------------------------------------------

#' Base ggplot2 theme for publication figures
#'
#' @param fixed_aspect Logical. Enforce aspect.ratio = 1 (default TRUE).
#'
#' @return A ggplot2 theme object.
plot_theme <- function(fixed_aspect = TRUE) {
  base <- ggplot2::theme(
    panel.background  = ggplot2::element_rect(fill = "white", colour = NA),
    plot.background   = ggplot2::element_rect(fill = "white", colour = NA),
    plot.title        = ggplot2::element_blank(),
    axis.text.y       = ggplot2::element_text(margin = ggplot2::margin(r = 0.5, unit = "mm")),
    axis.text.x       = ggplot2::element_text(margin = ggplot2::margin(t = 0.6, unit = "mm")),
    axis.title.x      = ggplot2::element_text(margin = ggplot2::margin(t = 0.3, unit = "mm")),
    axis.title.y      = ggplot2::element_text(margin = ggplot2::margin(r = 0.6, unit = "mm")),
    axis.title        = ggplot2::element_text(size = 7),
    axis.text         = ggplot2::element_text(size = 5, colour = "black"),
    axis.line         = ggplot2::element_line(linewidth = 0.4),
    axis.ticks        = ggplot2::element_line(lineend = "square", linewidth = 0.4,
                                               color = "black"),
    axis.ticks.length = ggplot2::unit(0.5, "mm"),
    legend.key.size   = ggplot2::unit(0.3, "cm"),
    legend.title      = ggplot2::element_text(size = 5,
                                               margin = ggplot2::margin(l = 1)),
    legend.text       = ggplot2::element_text(size = 5,
                                               margin = ggplot2::margin(l = 1)),
    legend.position   = "right",
    legend.box.margin = ggplot2::margin(-15, -15, -15, -15),
    legend.margin     = ggplot2::margin(0, 0, 0, 0),
    legend.justification = "left",
    legend.spacing.x  = ggplot2::unit(0, "mm"),
    legend.key.width  = ggplot2::unit(1, "mm"),
    text              = ggplot2::element_text(colour = "black")
  )
  if (fixed_aspect) base <- base + ggplot2::theme(aspect.ratio = 1)
  base
}

#' Add a text annotation at publication font size
#'
#' @param ... Arguments passed to \code{ggplot2::annotate()}.
#' @param size Numeric. Text size in ggplot units (default 5/.pt).
general_annotation <- function(..., size = 5 / .pt) {
  ggplot2::annotate(..., size = size)
}

#' Remove the stroke from the first layer of a ggplot
#'
#' @param figure A ggplot object.
#'
#' @return The modified ggplot object.
remove_stroke <- function(figure) {
  figure$layers[[1]]$aes_params$stroke <- 0
  figure
}

# Standard UMAP axis labels
UMAP_labs <- ggplot2::labs(x = "UMAP 1", y = "UMAP 2")

#' Guide settings for a continuous colour legend
guide_continuous <- function(legend_title = "Title",
                              barwidth     = ggplot2::unit(0.2, "lines"),
                              barheight    = 4) {
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      title          = legend_title,
      title.position = "left",
      barwidth       = barwidth,
      barheight      = barheight,
      direction      = "vertical",
      title.hjust    = 0.5,
      label.theme    = ggplot2::element_text(margin = ggplot2::margin(l = 1)),
      title.theme    = ggplot2::element_text(angle = 90, vjust = 0.5)
    )
  )
}

#' Guide settings for a discrete colour / fill legend
discrete_guide_settings <- function(ncol = 1) {
  ggplot2::guide_legend(
    direction      = "vertical",
    label.theme    = ggplot2::element_text(margin = ggplot2::margin(r = -3, l = 1),
                                            vjust = 0.5),
    override.aes   = list(size = 1),
    keyheight      = ggplot2::unit(0.5, "lines"),
    keywidth       = ggplot2::unit(0.2, "lines"),
    label.position = "right",
    label.hjust    = 0,
    label.vjust    = 0.5,
    ncol           = ncol
  )
}

guide_discrete <- function(ncol = 1) {
  ggplot2::guides(
    color = discrete_guide_settings(ncol),
    fill  = discrete_guide_settings(ncol)
  )
}

#' Combined continuous colour + discrete size guide
guide_combined <- function(legend_title = "Title",
                            barwidth     = ggplot2::unit(0.2, "lines"),
                            barheight    = 4,
                            size_ncol    = 1) {
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      title          = legend_title,
      title.position = "left",
      barwidth       = barwidth,
      barheight      = barheight,
      direction      = "vertical",
      title.hjust    = 0.5,
      label          = FALSE,
      title.theme    = ggplot2::element_text(angle = 90, vjust = 0.5)
    ),
    size = ggplot2::guide_legend(
      direction      = "vertical",
      label.theme    = ggplot2::element_text(margin = ggplot2::margin(r = -3, l = 1),
                                              vjust = 0.5),
      override.aes   = list(size = 1),
      label.position = "right",
      ncol           = size_ncol
    )
  )
}

# Colour palette (colourblind-friendly, tested across all figures)
custom_colors <- c(
  "#0D0887", "#ED7953", "#FDB32F", "#81D8D0", "#C8A2C8",
  "#66C2A5", "#8DA0CB", "#9C179E", "#E78AC3", "#A6D854",
  "#FFD92F", "#5E4FA2", "#3288BD", "#F46D43", "#FDAE61",
  "#1A9850", "#7570B3", "#D95F02", "#E6AB02", "#66A61E",
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#BC80BD",
  "#FB8072", "#80B1D3", "#FDB462", "#B3DE69"
)


# Plotting functions -------------------------------------------------------

#' Ridge plot of ambient load per cell across samples
ambient_load_ridge_plot <- function(ambient_load, avg_ambient_load) {
  ggplot2::ggplot(ambient_load, ggplot2::aes(x = ambient_percent, y = group,
                                              fill = group)) +
    ggridges::geom_density_ridges(alpha = 0.6, linewidth = 0.1) +
    ggplot2::geom_text(
      data = avg_ambient_load,
      ggplot2::aes(x = 0, y = group,
                   label = paste0(value, "%")),
      hjust = 0, vjust = -0.4, size = 2
    ) +
    plot_theme() +
    ggplot2::scale_fill_manual(values = custom_colors) +
    ggplot2::labs(x = "Ambient load (%)", y = "Percent of cells") +
    ggplot2::theme(legend.position = "none")
}

#' Histogram of ambient load for single-sample datasets
make_histogram_of_ambient_load <- function(ambient_load, avg_ambient_load) {
  hist_obj <- graphics::hist(ambient_load$ambient_percent, breaks = 100, plot = FALSE)
  hist_df  <- data.frame(
    mids   = hist_obj$mids,
    counts = hist_obj$counts / sum(hist_obj$counts)
  )
  ggplot2::ggplot(hist_df, ggplot2::aes(x = mids, y = counts)) +
    ggplot2::geom_col(fill = custom_colors[1]) +
    ggplot2::labs(x = "Ambient load (%)", y = "Percent of cells") +
    ggplot2::geom_vline(xintercept = 5) +
    general_annotation("text", x = 4, y = 0.2,
                        label = paste0("Avg. ambient load:\n", avg_ambient_load$value, "%"),
                        hjust = 1, vjust = 1) +
    plot_theme() +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 0, unit = "mm"))
    )
}

make_df_of_ambient_load_vs_total_background_expression <- function(
    species = c("human", "mouse"), background, ambient_counts) {
  species <- match.arg(species)
  data.frame(
    log_background_expression = log1p(background[match(names(rowSums(ambient_counts)), names(background))]),
    log_ambient_expression    = log1p(rowSums(ambient_counts)),
    species                   = species
  )
}

make_df_of_ambient_load_vs_total_endogenous_expression <- function(
    species = c("human", "mouse"), endogenous_counts, ambient_counts) {
  species <- match.arg(species)
  data.frame(
    log_endogenous_expression = log1p(rowSums(endogenous_counts)),
    log_ambient_expression    = log1p(rowSums(ambient_counts)),
    species                   = species
  )
}

make_df_of_amb_load_vs_endogenous_expression_frequency <- function(
    raw_human_endogenous, raw_mouse_endogenous,
    raw_human_ambient,   raw_mouse_ambient) {
  df_h <- data.frame(
    endogenous_expression_frequency = rowSums(raw_mouse_endogenous > 0) / ncol(raw_mouse_endogenous),
    log_ambient_expression          = log1p(rowSums(raw_human_ambient)),
    species                         = "human"
  )
  df_m <- data.frame(
    endogenous_expression_frequency = rowSums(raw_human_endogenous > 0) / ncol(raw_human_endogenous),
    log_ambient_expression          = log1p(rowSums(raw_mouse_ambient)),
    species                         = "mouse"
  )
  rbind(df_h, df_m)
}

make_df_of_ambient_fraction_vs_endogenous_UMIs <- function(
    species = c("human", "mouse"), endogenous_counts, ambient_counts, raw_counts) {
  species <- match.arg(species)
  data.frame(
    log_endogenous_UMIs = log1p(colSums(endogenous_counts)),
    ambient_fraction    = colSums(ambient_counts) / colSums(raw_counts),
    species             = species
  )
}

make_df_of_ambient_fraction_vs_total_endogenous_features <- function(
    endogenous_counts, ambient_counts, raw_counts,
    species = c("human", "mouse")) {
  species <- match.arg(species)
  data.frame(
    log_endogenous_features = log1p(colSums(endogenous_counts > 0)),
    ambient_fraction        = colSums(ambient_counts) / colSums(raw_counts),
    species                 = species
  )
}

make_barplot_of_explained_variance <- function(var_explained) {
  ggplot2::ggplot(var_explained, ggplot2::aes(x = parameters, y = percent_explained)) +
    ggplot2::geom_bar(stat = "identity", fill = custom_colors[[1]]) +
    plot_theme() +
    ggplot2::labs(x = "", y = "Explained variance (%)") +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

make_scatterplot_of_ambient_load_vs_total_background_expression <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = log_background_expression,
                                    y = log_ambient_expression)) +
    ggrastr::geom_point_rast(size = 0.4, stroke = 0,
                              color = custom_colors[[1]], raster.dpi = 1500) +
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    plot_theme() +
    ggplot2::labs(x = "log(background expr.)", y = "log(ambient expr.)")
}

make_scatterplot_of_amb_load_vs_total_endo_expression <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = log_endogenous_expression,
                                    y = log_ambient_expression)) +
    ggrastr::geom_point_rast(size = 0.4, stroke = 0,
                              color = custom_colors[[1]], raster.dpi = 1500) +
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    plot_theme() +
    ggplot2::labs(x = "log(endogenous expr.)", y = "log(ambient expr.)")
}

make_scatterplot_of_amb_load_vs_endo_expression_frequency <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = endogenous_expression_frequency,
                                    y = log_ambient_expression)) +
    ggrastr::geom_point_rast(size = 0.4, stroke = 0,
                              color = custom_colors[[1]], raster.dpi = 1500) +
    plot_theme() +
    ggplot2::labs(x = "Endogenous expr. freq.", y = "log(ambient expr.)")
}

make_scatterplot_of_amb_load_vs_exon_fraction <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = exon_fraction, y = log_ambient_expression)) +
    ggrastr::geom_point_rast(size = 0.4, stroke = 0,
                              color = custom_colors[[1]], raster.dpi = 1500) +
    plot_theme() +
    ggplot2::labs(x = "Exon fraction", y = "log(ambient expr.)")
}

make_scatterplot_of_ambient_fraction_vs_endogenous_UMIs <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = log_endogenous_UMIs, y = ambient_fraction)) +
    ggrastr::geom_point_rast(size = 0.4, stroke = 0,
                              color = custom_colors[[1]], raster.dpi = 1500) +
    plot_theme() +
    ggplot2::labs(x = "log UMIs endogenous", y = "Ambient fraction")
}

make_scatterplot_of_ambient_fraction_vs_total_endogenous_features <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = log_endogenous_features, y = ambient_fraction)) +
    ggrastr::geom_point_rast(size = 0.4, stroke = 0,
                              color = custom_colors[[1]], raster.dpi = 1500) +
    plot_theme() +
    ggplot2::labs(x = "log detected features (endogenous)", y = "Ambient fraction")
}


# Biology evaluation ------------------------------------------------------

#' Collect and summarise biology evaluation results across datasets
#'
#' Normalises each metric relative to the RNA (uncorrected) assay, averages
#' across datasets, and produces a dot plot.
#'
#' @param ds Named list of data frames, one per dataset, each with Assay and
#'   metric columns.
#' @param columns Character vector of metric column names to include.
#' @param remove_rna Logical. Drop the RNA assay from the output (default TRUE).
#' @param subsetTo Character. Subset rows to this value of subsetColumn.
#' @param subsetColumn Character. Column name for subsetting.
#' @param reverse_score Character vector of metrics to multiply by -1.
#'
#' @return A list with elements: a ggplot dot plot and a supplemental data frame.
collectResults <- function(ds, columns = c("LISI", "ASW", "PCR", "ARI"),
                            remove_rna    = TRUE,
                            subsetTo      = NULL,
                            subsetColumn  = NULL,
                            reverse_score = NULL) {
  supplemental <- ds

  for (column in columns) {
    for (dsidx in seq_along(ds)) {
      data <- ds[[dsidx]]
      if (!is.null(subsetTo)) data <- data[data[, subsetColumn] == subsetTo, ]

      if ("Sample" %in% colnames(data)) {
        columns_keep <- "Sample"
        for (sample in unique(data$Sample)) {
          if (sum(is.na(data[data$Sample == sample, column])) < nrow(data[data$Sample == sample, ])) {
            div <- mean(data[data$Assay == "RNA" & data$Sample == sample, column], na.rm = TRUE)
            data[data$Sample == sample, column] <- log2(data[data$Sample == sample, column] / div)
            supplemental[[dsidx]][supplemental[[dsidx]]$Sample == sample, column] <-
              log2(supplemental[[dsidx]][supplemental[[dsidx]]$Sample == sample, column] / div)
          }
        }
      } else {
        columns_keep <- c()
        if (sum(is.na(data[, column])) < nrow(data)) {
          div <- mean(data[data$Assay == "RNA", column], na.rm = TRUE)
          data[, column] <- log2(data[, column] / div)
          supplemental[[dsidx]][, column] <- log2(supplemental[[dsidx]][, column] / div)
        }
      }

      rel <- tapply(data[, column], data$Assay, mean, na.rm = TRUE)
      rel <- rel[match(unique(data$Assay), names(rel))]

      if (dsidx == 1) dat <- t(as.matrix(rel)) else dat <- rbind(dat, rel)
    }

    df_tmp <- data.frame(apply(dat, 2, mean, na.rm = TRUE))
    colnames(df_tmp) <- column
    if (column == columns[1]) df_mean <- df_tmp else df_mean <- cbind(df_mean, df_tmp)
  }

  if (remove_rna) df_mean <- df_mean[rownames(df_mean) != "RNA", ]

  if (!is.null(reverse_score)) {
    for (sc in reverse_score) df_mean[, sc] <- -1 * df_mean[, sc]
  }

  method_order <- names(sort(rank(-rowMeans(df_mean), ties.method = "first")))

  df_long <- df_mean |>
    tibble::rownames_to_column("Method") |>
    tidyr::pivot_longer(-Method, names_to = "Metric", values_to = "Mean")
  df_long$Metric <- factor(df_long$Metric, levels = columns)
  df_long$Method <- factor(df_long$Method, levels = rev(method_order))

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Metric, y = Method)) +
    ggplot2::geom_point(ggplot2::aes(color = Mean, size = abs(Mean))) +
    ggplot2::scale_size(range = c(2, 8)) +
    ggplot2::scale_color_gradient2(low = "#0571b0", mid = "grey", high = "#ca0020",
                                    midpoint = 0) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = "Metric", y = "Method", color = "Mean", size = "Mean")

  for (dsidx in seq_along(supplemental)) {
    supplemental[[dsidx]]$Dataset <- names(supplemental)[dsidx]
    supplemental[[dsidx]] <- supplemental[[dsidx]][,
      c(columns_keep, "Assay", "Dataset", columns)
    ]
    if (remove_rna)
      supplemental[[dsidx]] <- supplemental[[dsidx]][supplemental[[dsidx]]$Assay != "RNA", ]
  }
  supplemental <- do.call("rbind", supplemental)

  list(p, supplemental)
}


#' Evaluate label transfer quality using Azimuth
#'
#' @param object A Seurat object.
#' @param reference_path Character. Subdirectory under
#'   data/raw/reference_sources/Azimuth/ containing the reference.
#' @param convertTo Character. "human" or "mouse" (ortholog conversion).
#'   NULL to skip conversion.
#' @param subsetTo Character. Value of subset_label to retain. NULL = no subset.
#' @param sample_label Character. Metadata column for sample IDs (default
#'   "orig.ident").
#' @param subset_label Character. Metadata column for subsetting (default
#'   "origin").
#' @param verbose Logical (default TRUE).
#' @param split_genes_by_origin Logical. Separate human/mouse gene prefixes
#'   before ortholog conversion (default TRUE).
#'
#' @return A data frame of per-sample, per-assay evaluation metrics.
evaluateLabelTransfer <- function(object,
                                   reference_path,
                                   convertTo            = "human",
                                   subsetTo             = NULL,
                                   sample_label         = "orig.ident",
                                   subset_label         = "origin",
                                   verbose              = TRUE,
                                   split_genes_by_origin = TRUE) {
  azimuth_path <- here::here("data/raw/reference_sources/Azimuth", reference_path)

  first   <- TRUE
  assays  <- names(object@assays)
  samples <- unique(object@meta.data[, sample_label])

  for (sample in samples) {
    if (verbose) message("Processing sample: ", sample)
    if (sum(object@meta.data[, sample_label] == sample) < 200) next

    for (assay in assays) {
      if (verbose) message("\tProcessing assay: ", assay)
      query <- object[, object@meta.data[, sample_label] == sample]
      if (!is.null(subsetTo))
        query <- query[, query@meta.data[, subset_label] == subsetTo]

      query[["RNA"]] <- Seurat::CreateAssayObject(
        Seurat::GetAssayData(query, assay = assay, layer = "counts")
      )
      for (del in assays[assays != "RNA"]) query[[del]] <- NULL

      # Ortholog conversion
      if (!is.null(convertTo)) {
        if (split_genes_by_origin) {
          counts      <- Seurat::GetAssayData(query, layer = "counts")
          all_genes   <- rownames(query)
          human_genes <- all_genes[grep("^GRCh38-", all_genes)]
          mouse_genes <- all_genes[grep("^mm10-",   all_genes)]
          human_cnts  <- counts[human_genes, ]
          mouse_cnts  <- counts[mouse_genes, ]
          rownames(human_cnts) <- substr(human_genes, 8, nchar(human_genes))
          rownames(mouse_cnts) <- substr(mouse_genes, 6, nchar(mouse_genes))

          if (convertTo == "mouse") {
            conv <- suppressMessages(orthogene::convert_orthologs(
              data.frame(Gene = rownames(mouse_cnts)),
              gene_input = "Gene", gene_output = "columns",
              input_species = "mouse", output_species = "human", verbose = FALSE
            ))
            human_cnts <- human_cnts[rownames(human_cnts) %in% conv$ortholog_gene, ]
            miss <- conv[!(conv$ortholog_gene %in% rownames(human_cnts)), "ortholog_gene"]
            if (length(miss) > 0) {
              human_cnts <- Matrix::rbind2(human_cnts,
                Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                     dims = c(length(miss), ncol(human_cnts)),
                                     dimnames = list(miss, colnames(human_cnts))))
            }
            human_cnts <- human_cnts[match(conv$ortholog_gene, rownames(human_cnts)), ]
            rownames(human_cnts) <- conv$input_gene
            mouse_cnts <- mouse_cnts[match(conv$input_gene, rownames(mouse_cnts)), ]
            converted_counts <- human_cnts + mouse_cnts
          } else {
            conv <- suppressMessages(orthogene::convert_orthologs(
              data.frame(Gene = rownames(human_cnts)),
              gene_input = "Gene", gene_output = "columns",
              input_species = "human", output_species = "mouse", verbose = FALSE
            ))
            mouse_cnts <- mouse_cnts[rownames(mouse_cnts) %in% conv$ortholog_gene, ]
            miss <- conv[!(conv$ortholog_gene %in% rownames(mouse_cnts)), "ortholog_gene"]
            if (length(miss) > 0) {
              mouse_cnts <- Matrix::rbind2(mouse_cnts,
                Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                     dims = c(length(miss), ncol(mouse_cnts)),
                                     dimnames = list(miss, colnames(mouse_cnts))))
            }
            mouse_cnts <- mouse_cnts[match(conv$ortholog_gene, rownames(mouse_cnts)), ]
            rownames(mouse_cnts) <- conv$input_gene
            human_cnts <- human_cnts[match(conv$input_gene, rownames(human_cnts)), ]
            converted_counts <- human_cnts + mouse_cnts
          }
        } else {
          converted_counts <- Seurat::GetAssayData(query, layer = "counts")
          if (convertTo == "mouse") {
            conv <- suppressMessages(orthogene::convert_orthologs(
              data.frame(Gene = rownames(query)), gene_input = "Gene",
              gene_output = "columns", input_species = "human",
              output_species = "mouse", verbose = FALSE
            ))
          } else {
            conv <- suppressMessages(orthogene::convert_orthologs(
              data.frame(Gene = rownames(query)), gene_input = "Gene",
              gene_output = "columns", input_species = "mouse",
              output_species = "human", verbose = FALSE
            ))
          }
          converted_counts <- converted_counts[rownames(converted_counts) %in% conv$input_gene, ]
          converted_counts <- converted_counts[match(conv$input_gene, rownames(converted_counts)), ]
          rownames(converted_counts) <- conv$ortholog_gene
        }
        suppressWarnings(query[["RNA"]] <- Seurat::CreateAssayObject(counts = converted_counts))
      }

      query     <- query[, Matrix::colSums(query@assays$RNA@counts) > 0]
      reference <- suppressWarnings(Azimuth::LoadReference(azimuth_path)$map)
      suppressWarnings(reference[["SCT"]] <- reference@assays$refAssay)

      query <- suppressMessages(Azimuth::ConvertGeneNames(
        object          = query,
        reference.names = rownames(reference),
        homolog.table   = "https://seurat.nygenome.org/azimuth/references/homologs.rds"
      ))
      suppressWarnings(query <- SeuratObject::CreateSeuratObject(
        Seurat::CreateAssayObject(Seurat::GetAssayData(query, assay = "RNA", layer = "counts")),
        meta.data = query@meta.data
      ))
      query <- query[, Matrix::colSums(query@assays$RNA@counts) > 0]
      query <- Seurat::SCTransform(query, verbose = FALSE)

      calcn <- as.data.frame(Seurat:::CalcN(query[["RNA"]]))
      colnames(calcn) <- paste(colnames(calcn), "RNA", sep = "_")
      query <- Seurat::AddMetaData(query, calcn)
      query <- Seurat::PercentageFeatureSet(query, pattern = "^MT-",
                                             col.name = "percent.mt", assay = "RNA")

      dims    <- as.double(slot(reference, "neighbors")$refdr.annoy.neighbors@alg.info$ndim)
      anchors <- suppressMessages(suppressWarnings(Seurat::FindTransferAnchors(
        reference            = reference,
        query                = query,
        k.filter             = NA,
        reference.neighbors  = "refdr.annoy.neighbors",
        reference.assay      = "SCT",
        query.assay          = "SCT",
        reference.reduction  = "refDR",
        normalization.method = "SCT",
        features             = rownames(Seurat::Loadings(reference[["refDR"]])),
        dims                 = 1:dims,
        n.trees              = 20,
        mapping.score.k      = min(100, floor(ncol(query) / 4)),
        verbose              = FALSE
      )))

      query <- suppressWarnings(Seurat::TransferData(
        reference     = reference,
        query         = query,
        query.assay   = "SCT",
        dims          = 1:dims,
        anchorset     = anchors,
        refdata       = list(annotation.l1 = reference[["annotation.l1", drop = TRUE]]),
        n.trees       = 20,
        store.weights = TRUE,
        k.weight      = min(50, ceiling(nrow(anchors@anchors) / 4)),
        verbose       = FALSE
      ))
      query <- Seurat::AddMetaData(
        query,
        metadata = Seurat::MappingScore(anchors = anchors, ndim = dims, verbose = FALSE,
                                         kanchors = min(50, ceiling(nrow(anchors@anchors) / 4)),
                                         ksmooth  = min(100, floor(ncol(query) / 4))),
        col.name = "mapping.score"
      )

      map_scores    <- sum(query$mapping.score >= 0.9) / ncol(query)
      annot_fractions <- sum(query$predicted.annotation.l1.score >= 0.9) / ncol(query)

      query <- query[, query$predicted.annotation.l1 %in%
                       names(which(table(query$predicted.annotation.l1) >=
                                     max(3, ceiling(ncol(query) / 100)))), ]
      Seurat::DefaultAssay(query) <- "RNA"
      query <- Seurat::NormalizeData(query, verbose = FALSE)
      query <- suppressWarnings(Seurat::FindVariableFeatures(query, verbose = FALSE))
      query <- Seurat::ScaleData(query, verbose = FALSE)
      query <- suppressWarnings(Seurat::RunPCA(query, verbose = FALSE))

      if (length(unique(query$predicted.annotation.l1)) > 1) {
        query <- Seurat::FindNeighbors(query, reduction = "pca", verbose = FALSE)
        aris  <- sapply(seq(0.05, 1.0, by = 0.05), function(res) {
          query <- Seurat::FindClusters(query, res = res, verbose = FALSE)
          aricode::ARI(query$seurat_clusters, query$predicted.annotation.l1)
        })
        ARI <- max(aris)

        query_sce <- Seurat::as.SingleCellExperiment(query)
        query_sce <- suppressWarnings(suppressMessages(miloR::Milo(query_sce)))
        out <- capture.output(query_sce <- suppressWarnings(suppressMessages(
          miloDE::assign_neighbourhoods(query_sce, k = 20, order = 2,
                                        filtering = TRUE, reducedDim_name = "PCA",
                                        verbose = FALSE)
        )))
        query_nhoods <- suppressWarnings(miloR::nhoods(query_sce))

        purity <- purity_oe <- numeric(ncol(query_nhoods))
        for (nhood in seq_len(ncol(query_nhoods))) {
          obs_p  <- max(table(query[, query_nhoods[, nhood] == 1]$predicted.annotation.l1)) /
            sum(table(query[, query_nhoods[, nhood] == 1]$predicted.annotation.l1))
          ns     <- sum(query_nhoods[, nhood] == 1)
          labels <- query$predicted.annotation.l1
          maxes  <- vapply(seq_len(100), function(r) {
            set.seed(r); max(table(labels[sample.int(length(labels), ns)]))
          }, numeric(1))
          totals <- rep(ns, 100)
          purity[nhood]    <- obs_p
          purity_oe[nhood] <- log2(obs_p / mean(maxes / totals))
        }

        dd   <- distances::distances(query@reductions$pca@cell.embeddings[, 1:20])
        dd   <- distances::distance_matrix(dd)
        sil  <- cluster::silhouette(
          as.numeric(factor(query$predicted.annotation.l1)), dd
        )
        cASWs <- (mean(sil[, "sil_width"]) + 1) / 2
      } else {
        ARI <- cASWs <- purity <- purity_oe <- NA
      }

      query <- Seurat::NormalizeData(query, verbose = FALSE)
      query <- Seurat::ScaleData(query, verbose = FALSE, features = rownames(query))
      query <- subset(query, features = rownames(reference))
      cnts  <- Seurat::GetAssayData(query, layer = "RNA", slot = "data")
      RMSE  <- mean(vapply(unique(query$predicted.annotation.l1), function(ct) {
        ct_mat <- cnts[, query$predicted.annotation.l1 == ct, drop = FALSE]
        ct_mat <- ct_mat[rowSums(ct_mat > 0) >= ncol(ct_mat) * 0.1, , drop = FALSE]
        if (nrow(ct_mat) == 0) return(NA_real_)
        mean(sqrt(colMeans((t(ct_mat) - rowMeans(ct_mat))^2)))
      }, numeric(1)), na.rm = TRUE)

      df_tmp <- data.frame(
        Sample     = sample,
        Assay      = assay,
        nCells     = ncol(query),
        Mapping    = map_scores,
        Annotation = annot_fractions,
        Purity     = mean(purity,    na.rm = TRUE),
        Purity_OE  = mean(purity_oe, na.rm = TRUE),
        RMSE       = RMSE,
        cASW       = cASWs,
        ARI        = ARI,
        nLabels    = length(unique(query$predicted.annotation.l1))
      )
      if (first) { df_res <- df_tmp; first <- FALSE } else df_res <- rbind(df_res, df_tmp)
    }
  }
  df_res
}


#' Evaluate sample integration quality across decontamination methods
#'
#' @param object A Seurat object with all decontamination assays.
#' @param sample_label Character. Metadata column for sample IDs (default
#'   "orig.ident").
#' @param max_cells Integer. Max cells to use for ASW (default 20000).
#' @param verbose Logical (default TRUE).
#' @param known_markers Character vector of known marker gene names.
#'
#' @return A data frame of per-assay, per-integration-state metrics.
evaluateIntegration <- function(object,
                                 sample_label  = "orig.ident",
                                 max_cells     = 20000,
                                 verbose       = TRUE,
                                 known_markers = c()) {
  first  <- TRUE
  assays <- names(object@assays)

  for (assay in assays) {
    if (verbose) message("Processing assay: ", assay)
    obj_assay <- SeuratObject::CreateSeuratObject(
      Seurat::GetAssayData(object, assay = assay, layer = "counts"),
      meta.data = object@meta.data
    )
    obj_assay <- Seurat::NormalizeData(obj_assay, verbose = FALSE)
    obj_list  <- Seurat::SplitObject(obj_assay, split.by = sample_label)
    Seurat::VariableFeatures(obj_assay) <- Seurat::SelectIntegrationFeatures(
      obj_list, verbose = FALSE
    )
    obj_assay <- Seurat::ScaleData(obj_assay, verbose = FALSE)
    obj_assay <- Seurat::RunPCA(obj_assay, verbose = FALSE)
    obj_assay <- harmony::RunHarmony(obj_assay, sample_label, verbose = FALSE)

    for (method in c("Unintegrated", "Integrated")) {
      if (verbose) message("\tEvaluating ", method)
      reduction <- if (method == "Unintegrated") "pca" else "harmony"

      obj_assay <- Seurat::FindNeighbors(obj_assay, reduction = reduction,
                                          verbose = FALSE)
      aris <- sapply(seq(0.05, 1.0, by = 0.05), function(res) {
        obj_assay <- Seurat::FindClusters(obj_assay, res = res, verbose = FALSE)
        aricode::ARI(obj_assay$seurat_clusters, obj_assay@meta.data[, sample_label])
      })
      ARI <- max(aris)

      embeds <- obj_assay@reductions[[reduction]]@cell.embeddings[, 1:20]
      LISI   <- mean(
        (lisi::compute_lisi(embeds, obj_assay@meta.data, sample_label)[, 1] - 1) /
          (length(unique(obj_assay@meta.data[, sample_label])) - 1)
      )

      pca_data <- list(
        x    = embeds,
        sdev = obj_assay@reductions[[reduction]]@stdev
      )
      PCR <- kBET::pcRegression(
        pca_data, batch = as.numeric(factor(obj_assay@meta.data[, sample_label]))
      )$R2Var

      if (ncol(obj_assay) >= max_cells) {
        bASWs <- vapply(seq_len(5), function(rep) {
          set.seed(rep)
          sel     <- sample(colnames(obj_assay), max_cells)
          obj_sub <- subset(obj_assay, cells = sel)
          dd  <- distances::distance_matrix(
            distances::distances(obj_sub@reductions[[reduction]]@cell.embeddings[, 1:20])
          )
          sil <- cluster::silhouette(
            as.numeric(factor(obj_sub@meta.data[, sample_label])), dd
          )
          1 - mean(abs(sil[, "sil_width"]))
        }, numeric(1))
        bASW <- mean(bASWs)
      } else {
        dd   <- distances::distance_matrix(distances::distances(embeds))
        sil  <- cluster::silhouette(
          as.numeric(factor(obj_assay@meta.data[, sample_label])), dd
        )
        bASW <- 1 - mean(abs(sil[, "sil_width"]))
      }

      # Marker gene metrics
      W_obj <- suppressMessages(
        Seurat::FindNeighbors(obj_assay, verbose = FALSE, return.neighbor = TRUE,
                               reduction = reduction)@neighbors$RNA.nn
      )
      W <- Matrix::sparseMatrix(
        i = rep(seq_len(nrow(W_obj)), each = ncol(W_obj@nn.idx)),
        j = as.vector(t(W_obj@nn.idx)),
        x = 1, dims = c(nrow(W_obj), nrow(W_obj))
      )
      W_sym <- (W + Matrix::t(W)) / 2
      W_sum <- sum(W_sym)

      obj_markers <- subset(obj_assay, features = known_markers)
      exprs       <- Seurat::GetAssayData(obj_markers, layer = "data")
      genes       <- rownames(obj_markers)

      morans <- entropy <- ginis <- numeric(0)
      for (marker in known_markers) {
        if (!(marker %in% genes)) next
        expr <- as.matrix(exprs[genes == marker, , drop = FALSE])
        z    <- expr - mean(expr)
        den  <- sum(z^2)
        num  <- as.numeric(Matrix::t(z) %*% (W_sym %*% z))
        m_obs <- (length(z) / W_sum) * (num / den)
        m_null <- vapply(seq_len(100), function(s) {
          z_s <- sample(z, length(z))
          (length(z_s) / W_sum) * (as.numeric(Matrix::t(z_s) %*% (W_sym %*% z_s)) / den)
        }, numeric(1))
        morans <- c(morans, (m_obs - mean(m_null)) / sd(m_null))

        WX      <- Matrix::Diagonal(x = expr) %*% W
        col_s   <- Matrix::colSums(WX)
        entr_v  <- vapply(which(col_s > 0), function(i) {
          p <- expr[as.logical(W[, i])]
          p <- p / sum(p)
          -sum(p * log2(p + 1e-12)) / log2(length(p))
        }, numeric(1))
        entropy <- c(entropy, mean(entr_v[!is.infinite(entr_v)], na.rm = TRUE))
        ginis   <- c(ginis, DescTools::Gini(sort(expr)))
      }

      df_tmp <- data.frame(
        Assay   = assay,
        State   = method,
        LISI    = LISI,
        ASW     = bASW,
        PCR     = 1 - PCR,
        ARI     = 1 - ARI,
        Gini    = mean(ginis,   na.rm = TRUE),
        Entropy = mean(entropy, na.rm = TRUE),
        Morans  = mean(morans,  na.rm = TRUE)
      )
      if (first) { df_res <- df_tmp; first <- FALSE } else df_res <- rbind(df_res, df_tmp)
    }
  }
  df_res
}


#' Evaluate per-sample biology across decontamination methods
#'
#' @param object A Seurat object with all decontamination assays.
#' @param sample_label Character. Metadata column for sample IDs (default
#'   "orig.ident").
#' @param species_label Character. Metadata column for species / cell type
#'   labels (default "origin").
#' @param known_markers Character vector of known marker genes.
#' @param evalKnown Logical. Evaluate known marker metrics (default TRUE).
#' @param evalNhoods Logical. Evaluate Milo neighbourhoods (default TRUE).
#' @param evalDeNovo Logical. Evaluate de novo marker discovery (default TRUE).
#' @param verbose Logical (default TRUE).
#'
#' @return A data frame of per-sample, per-assay biology metrics.
evaluatePerSampleBiology <- function(object,
                                      sample_label  = "orig.ident",
                                      species_label = "origin",
                                      known_markers = c(),
                                      evalKnown     = TRUE,
                                      evalNhoods    = TRUE,
                                      evalDeNovo    = TRUE,
                                      verbose       = TRUE) {
  first   <- TRUE
  assays  <- names(object@assays)
  samples <- unique(object@meta.data[, sample_label])

  for (sample in samples) {
    if (verbose) message("Processing sample: ", sample)
    obj_sub <- object[, object@meta.data[, sample_label] == sample]

    for (assay in assays) {
      if (verbose) message("\tProcessing assay: ", assay)
      obj_a <- SeuratObject::CreateSeuratObject(
        Seurat::CreateAssayObject(
          counts = Seurat::GetAssayData(obj_sub, assay = assay, layer = "counts")
        ),
        meta.data = obj_sub@meta.data
      )
      obj_a <- Seurat::NormalizeData(obj_a, verbose = FALSE)
      obj_a <- suppressWarnings(Seurat::FindVariableFeatures(obj_a, verbose = FALSE))
      obj_a <- Seurat::ScaleData(obj_a, verbose = FALSE)
      obj_a <- Seurat::RunPCA(obj_a, verbose = FALSE)

      if (evalNhoods) {
        obj_sce <- Seurat::as.SingleCellExperiment(obj_a)
        obj_sce <- suppressWarnings(suppressMessages(miloR::Milo(obj_sce)))
        out <- capture.output(obj_sce <- suppressWarnings(suppressMessages(
          miloDE::assign_neighbourhoods(obj_sce, k = 20, order = 2,
                                        filtering = TRUE, reducedDim_name = "PCA",
                                        verbose = FALSE)
        )))
        obj_nhoods <- suppressWarnings(miloR::nhoods(obj_sce))
        rm(obj_sce); gc()

        LISI       <- lisi::compute_lisi(obj_a@reductions$pca@cell.embeddings,
                                          obj_a@meta.data, species_label)
        n_species  <- length(unique(obj_a@meta.data[, species_label]))
        species_LISI <- mean((n_species - LISI[, 1]) / (n_species - 1))

        sp_max <- sp_tot <- sp_emax <- sp_etot <- numeric(ncol(obj_nhoods))
        mk_max <- mk_tot <- numeric(ncol(obj_nhoods))

        for (nhood in seq_len(ncol(obj_nhoods))) {
          obj_a$nhood <- obj_nhoods[, nhood]
          sp_max[nhood] <- max(table(obj_a[, obj_a$nhood == 1]@meta.data[, species_label]))
          sp_tot[nhood] <- sum(table(obj_a[, obj_a$nhood == 1]@meta.data[, species_label]))

          ns     <- sum(obj_nhoods[, nhood] == 1)
          labels <- obj_a@meta.data[, species_label]
          maxes  <- vapply(seq_len(100), function(r) {
            set.seed(r); max(table(labels[sample.int(length(labels), ns)]))
          }, numeric(1))
          sp_emax[nhood] <- mean(maxes)
          sp_etot[nhood] <- ns

          if (evalDeNovo) {
            nhood_sp <- names(which.max(
              table(obj_a[, obj_a$nhood == 1]@meta.data[, species_label])
            ))
            markers <- presto::wilcoxauc(obj_a, group_by = "nhood")
            markers <- markers[markers$group == 1 & markers$padj <= 0.05 &
                                 markers$logFC >= log2(1.1) & markers$pct_in >= 1, ]
            pat <- if (nhood_sp == "mouse") "^mm10-" else "^GRCh38-"
            mk_max[nhood] <- sum(grepl(pat, markers$feature))
            mk_tot[nhood] <- nrow(markers)
          }
        }

        nhood_purity    <- sum(sp_max) / sum(sp_tot)
        nhood_purity_oe <- log2(nhood_purity / (sum(sp_emax) / sum(sp_etot)))
        marker_purity   <- if (evalDeNovo) sum(mk_max) / sum(mk_tot) else NA
      } else {
        species_LISI <- nhood_purity <- nhood_purity_oe <- marker_purity <- NA
      }

      if (evalKnown) {
        W_obj <- suppressMessages(
          Seurat::FindNeighbors(obj_a, verbose = FALSE, return.neighbor = TRUE)@neighbors$RNA.nn
        )
        W     <- Matrix::sparseMatrix(
          i = rep(seq_len(nrow(W_obj)), each = ncol(W_obj@nn.idx)),
          j = as.vector(t(W_obj@nn.idx)),
          x = 1, dims = c(nrow(W_obj), nrow(W_obj))
        )
        W_sym <- (W + Matrix::t(W)) / 2
        W_sum <- sum(W_sym)

        obj_mk <- subset(obj_a, features = known_markers)
        exprs  <- Seurat::GetAssayData(obj_mk, layer = "data")
        genes  <- rownames(obj_mk)
        shuffles <- lapply(seq_len(100), function(s) {
          set.seed(s); sample(seq_len(ncol(obj_a)))
        })

        morans <- entropy <- ginis <- numeric(0)
        for (marker in known_markers) {
          if (!(marker %in% genes)) next
          expr  <- as.matrix(exprs[genes == marker, , drop = FALSE])
          z     <- expr - mean(expr)
          den   <- sum(z^2)
          num   <- as.numeric(Matrix::t(z) %*% (W_sym %*% z))
          m_obs <- (length(z) / W_sum) * (num / den)
          m_null <- vapply(seq_len(100), function(i) {
            z_s <- z[shuffles[[i]], , drop = FALSE]
            (length(z_s) / W_sum) *
              (as.numeric(Matrix::t(z_s) %*% (W_sym %*% z_s)) / den)
          }, numeric(1))
          morans <- c(morans, (m_obs - mean(m_null)) / sd(m_null))

          WX     <- Matrix::Diagonal(x = expr) %*% W
          col_s  <- Matrix::colSums(WX)
          entr_v <- vapply(which(col_s > 0), function(i) {
            p <- expr[as.logical(W[, i])]
            p <- p / sum(p)
            -sum(p * log2(p + 1e-12)) / log2(length(p))
          }, numeric(1))
          entropy <- c(entropy, mean(entr_v[!is.infinite(entr_v)], na.rm = TRUE))
          ginis   <- c(ginis, DescTools::Gini(sort(expr)))
        }
        morans_mean  <- mean(morans,  na.rm = TRUE)
        entropy_mean <- mean(entropy, na.rm = TRUE)
        gini_mean    <- mean(ginis,   na.rm = TRUE)
      } else {
        morans_mean <- entropy_mean <- gini_mean <- NA
      }

      df_tmp <- data.frame(
        Sample    = sample,
        Assay     = assay,
        Purity    = nhood_purity,
        Purity_OE = nhood_purity_oe,
        Species   = species_LISI,
        Markers   = marker_purity,
        Entropy   = entropy_mean,
        Morans    = morans_mean,
        Gini      = gini_mean
      )
      if (first) { df_res <- df_tmp; first <- FALSE } else df_res <- rbind(df_res, df_tmp)
    }
  }
  df_res
}


#' Define marker genes that are consistently clean across simulation samples
#'
#' @param object A Seurat object (without ambient RNA, e.g. synt_seu_wo_ambient).
#' @param min_lfc Numeric. Minimum log2 fold-change (default log2(1.5)).
#' @param max_pct Numeric. Max % expression in other cells (default 50).
#' @param min_target_auc Numeric. Minimum AUC in target cell type (default 0.7).
#' @param max_background_auc Numeric. Max AUC in background (default 0.6).
#' @param label_key Character. Metadata column for cell type (default "celltype").
#'
#' @return Named list of marker gene vectors per cell type.
defineCleanMarkers <- function(object,
                                min_lfc          = log2(1.5),
                                max_pct          = 50,
                                min_target_auc   = 0.7,
                                max_background_auc = 0.6,
                                label_key        = "celltype") {
  first   <- TRUE
  samples <- unique(object$orig.ident)

  for (sample in samples) {
    obj <- subset(object, orig.ident == sample)
    seu <- Seurat::CreateSeuratObject(
      Seurat::CreateAssayObject(
        counts = Seurat::GetAssayData(obj, layer = "counts", assay = "RNA")
      ),
      meta.data = obj@meta.data
    )
    seu <- Seurat::NormalizeData(seu, verbose = FALSE)
    markers_tmp <- presto::wilcoxauc(seu, group_by = label_key)
    markers_tmp <- markers_tmp[
      markers_tmp$logFC  >= min_lfc &
        markers_tmp$pct_out <= max_pct &
        markers_tmp[, 1]  %in% markers_tmp[markers_tmp$auc >= min_target_auc, 1] &
        markers_tmp[, 1]  %in% names(which(
          table(markers_tmp[markers_tmp$auc >= max_background_auc, 1]) == 1
        )), 
    ]
    if (first) { marker_res <- markers_tmp; first <- FALSE } else
      marker_res <- rbind(marker_res, markers_tmp)
  }

  marker_res$ID <- paste(marker_res[, 1], marker_res[, 2], sep = "_")
  marker_res     <- marker_res[
    marker_res$ID %in% names(which(table(marker_res$ID) == length(samples))), 
  ]
  celltypes <- unique(object@meta.data[, label_key])
  setNames(
    lapply(celltypes, function(ct)
      unique(marker_res[marker_res$group == ct, 1])),
    celltypes
  )
}


#' Evaluate decontamination on the synthetic dataset (ground truth comparison)
#'
#' @param object Named list of two Seurat objects: "clean" (without ambient)
#'   and "unclean" (with ambient), both from sample 3.
#' @param label_key Character. Metadata column for cell type (default "celltype").
#' @param marker_list Named list of marker genes from \code{defineCleanMarkers()}.
#'
#' @return A list with elements: a summary data frame, a list of UMAP ggplots,
#'   and a vector of cluster counts.
evaluateSynthetic <- function(object, label_key = "celltype", marker_list = NULL) {
  types <- names(object)
  first <- TRUE
  umaps <- list()
  nc    <- c()

  for (type in types) {
    obj    <- object[[type]]
    assays <- names(obj@assays)
    celltypes <- unique(obj@meta.data[, label_key])

    for (assay in assays) {
      seu <- Seurat::CreateSeuratObject(
        Seurat::CreateAssayObject(
          counts = Seurat::GetAssayData(obj, layer = "counts", assay = assay)
        ),
        meta.data = obj@meta.data
      )
      seu <- Seurat::NormalizeData(seu, verbose = FALSE)
      seu <- Seurat::FindVariableFeatures(seu, verbose = FALSE)
      seu <- Seurat::ScaleData(seu, verbose = FALSE)
      seu <- Seurat::RunPCA(seu, verbose = FALSE)
      seu <- Seurat::FindNeighbors(seu, verbose = FALSE)
      seu <- UCell::AddModuleScore_UCell(seu, marker_list)
      seu <- Seurat::RunUMAP(seu, 1:20, verbose = FALSE)

      aris <- sapply(seq(0.05, 1, by = 0.05), function(res) {
        seu <- Seurat::FindClusters(seu, res = res, verbose = FALSE)
        aricode::ARI(seu$seurat_clusters, seu$celltype)
      })
      ari_max <- max(aris)
      best_res <- seq(0.05, 1, by = 0.05)[min(which(aris == ari_max))]
      seu <- Seurat::FindClusters(seu, res = best_res, verbose = FALSE)

      # Map clusters to cell types
      seu$seurat_clusters <- as.character(seu$seurat_clusters)
      tbl     <- table(seu$celltype, seu$seurat_clusters)
      mapping <- apply(tbl, 2, which.max)
      for (idx in seq_along(mapping)) {
        seu@meta.data[seu$seurat_clusters == names(mapping)[idx], "seurat_clusters"] <-
          rownames(tbl)[mapping[idx]]
      }
      nc <- c(nc, ncol(tbl))
      seu$label <- ifelse(seu$seurat_clusters == seu$celltype, seu$celltype, "False")

      umaps[[length(umaps) + 1]] <-
        Seurat::DimPlot(seu, group.by = "label") +
        ggplot2::ggtitle(paste0(type, " - ", assay)) +
        ggplot2::scale_color_manual(values = custom_colors) +
        plot_theme()

      # Metrics
      pca_emb         <- seu@reductions$pca@cell.embeddings[, 1:20]
      distances_list  <- c()
      module_scores   <- c()
      second_module   <- c()

      for (ct in celltypes) {
        X      <- pca_emb[seu$celltype == ct, ]
        X_norm <- X / sqrt(rowSums(X^2))
        med    <- pcaPP::l1median(X_norm)
        vn     <- med / sqrt(sum(med^2))
        sim    <- X_norm %*% vn
        distances_list <- c(distances_list, (1 - sim)[, 1])
        module_scores  <- c(module_scores,
          seu@meta.data[seu$celltype == ct, paste0(ct, "_UCell")])
        second_module  <- c(second_module,
          apply(seu@meta.data[seu$celltype == ct,
                               paste0(celltypes[celltypes != ct], "_UCell"),
                               drop = FALSE], 1, max))
      }

      df_tmp <- data.frame(
        Method      = assay,
        Type        = type,
        Gini        = DescTools::Gini(sort(distances_list)),
        ARI         = ari_max,
        Distance    = stats::median(distances_list),
        Module      = mean(module_scores),
        deltaModule = mean(module_scores - second_module)
      )
      if (first) { df_res <- df_tmp; first <- FALSE } else df_res <- rbind(df_res, df_tmp)
    }
  }

  list(df_res, umaps, nc)
}
