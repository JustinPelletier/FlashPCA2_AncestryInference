#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(FNN)
})

# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
# Args:
# 1) ref_pca_file      : FlashPCA2 reference eigenvectors
# 2) meta_file         : TSV, col1 = SampleID, col2 = Population/Group
# 3) study_id_file     : text file of Study sample IDs, or "NONE"
# 4) mode              : reference_only OR projection
# 5) k                 : kNN neighbors
# 6) n_pcs             : number of PCs to use for kNN
# 7) threshold_N       : min neighbor count to assign label
# 8) study_pca_file    : FlashPCA2 projected study PCs, required if mode=projection
#
# Example reference-only:
# Rscript PCA_flashpca.R reference_flashpca.eigenvectors meta.tsv NONE reference_only 10 4 6
#
# Example projection:
# Rscript PCA_flashpca.R reference_flashpca.eigenvectors meta.tsv study_ids.txt projection 10 4 6 study_flashpca.projected.txt

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript PCA_flashpca.R <ref_pca_file> <meta_file.tsv> <study_id_file|NONE> <mode> [k] [n_pcs] [threshold_N] [study_pca_file]")
}

ref_pca_file   <- args[1]
meta_file      <- args[2]
study_id_file  <- args[3]
mode           <- args[4]

k             <- if (length(args) >= 5) as.integer(args[5]) else 10L
n_pcs         <- if (length(args) >= 6) as.integer(args[6]) else 10L
threshold_N   <- if (length(args) >= 7) as.integer(args[7]) else 7L
study_pca_file <- if (length(args) >= 8) args[8] else NA_character_

if (!mode %in% c("reference_only", "projection")) {
  stop("mode must be either 'reference_only' or 'projection'")
}

if (mode == "projection" && (is.na(study_pca_file) || study_pca_file == "")) {
  stop("Projection mode requires a study_pca_file as argument 8.")
}

if (any(is.na(c(k, n_pcs, threshold_N)))) {
  stop("Bad numeric arguments for k, n_pcs, or threshold_N.")
}

dir.create("plot_PCA", showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
get_big_palette <- function(n) {
  okabe_ito <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
  kelly20 <- c("#F2F3F4","#222222","#F3C300","#875692","#F38400","#A1CAF1","#BE0032","#C2B280",
               "#848482","#008856","#E68FAC","#0067A5","#F99379","#604E97","#F6A600","#B3446C",
               "#DCD300","#882D17","#8DB600","#654522")
  tableau20 <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948",
                 "#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
                 "#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
                 "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")
  base <- unique(c(okabe_ito, kelly20, tableau20))
  if (n <= length(base)) return(base[seq_len(n)])
  c(base, grDevices::hcl(h = seq(15, 375, length.out = n - length(base) + 1)[-1], c = 65, l = 55))
}

clean_id <- function(x) {
  x <- as.character(x)
  x <- sub("^0_", "", x)
  x <- sub("^.*:", "", x)
  x
}

read_flashpca_coords <- function(file, dataset_label) {
  dt <- fread(file, data.table = TRUE)

  # Common FlashPCA2/PLINK-like formats:
  # FID IID PC1 PC2...
  # ID PC1 PC2...
  # sample PC1 PC2...
  nms <- names(dt)

  pc_cols <- grep("^(PC|V)[0-9]+$", nms, value = TRUE)

  # If PC columns are not named, assume first 2 cols are IDs and rest are PCs
  if (length(pc_cols) < 2) {
    if (ncol(dt) >= 4) {
      setnames(dt, 1:2, c("FID", "IID"))
      pc_n <- ncol(dt) - 2
      setnames(dt, 3:ncol(dt), paste0("PC", seq_len(pc_n)))
    } else {
      stop(sprintf("Could not infer PC columns in %s", file))
    }
  }

  nms <- names(dt)

  if (all(c("FID", "IID") %in% nms)) {
    dt[, FID := clean_id(FID)]
    dt[, IID := clean_id(IID)]
  } else if (all(c("#FID", "IID") %in% nms)) {
    setnames(dt, "#FID", "FID")
    dt[, FID := clean_id(FID)]
    dt[, IID := clean_id(IID)]
  } else if ("IID" %in% nms) {
    dt[, FID := IID]
    dt[, IID := clean_id(IID)]
    dt[, FID := clean_id(FID)]
  } else {
    setnames(dt, 1, "IID")
    dt[, FID := IID]
    dt[, IID := clean_id(IID)]
    dt[, FID := clean_id(FID)]
  }

  pc_cols <- grep("^PC\\d+$", names(dt), value = TRUE)

  if (length(pc_cols) < 2) {
    non_id <- setdiff(names(dt), c("FID", "IID"))
    numeric_cols <- non_id[sapply(dt[, ..non_id], is.numeric)]
    setnames(dt, numeric_cols, paste0("PC", seq_along(numeric_cols)))
    pc_cols <- grep("^PC\\d+$", names(dt), value = TRUE)
  }

  dt[, Dataset := dataset_label]
  dt[, c("FID", "IID", pc_cols, "Dataset"), with = FALSE]
}

# ----------------------------------------------------------------------
# Metadata
# ----------------------------------------------------------------------
meta <- fread(meta_file, sep = "\t", header = TRUE)
if (ncol(meta) < 2) stop("Metadata TSV must have at least two columns: ID and Group.")

setnames(meta, names(meta)[1:2], c("ID", "Group"))
meta[, ID := as.character(ID)]
meta[, Group := as.character(Group)]

# ----------------------------------------------------------------------
# Read FlashPCA2 coordinates
# ----------------------------------------------------------------------
ref_dt <- read_flashpca_coords(ref_pca_file, "Reference")

if (mode == "projection") {
  study_dt <- read_flashpca_coords(study_pca_file, "Study")
  dt <- rbindlist(list(ref_dt, study_dt), use.names = TRUE, fill = TRUE)
} else {
  dt <- ref_dt
}

pc_cols <- grep("^PC\\d+$", names(dt), value = TRUE)
n_pcs <- max(1L, min(n_pcs, length(pc_cols)))
use_pcs <- paste0("PC", seq_len(n_pcs))

# Optional ID filtering
if (mode == "projection" && study_id_file != "NONE") {
  study_ids <- fread(study_id_file, header = FALSE, col.names = "IID")[, as.character(IID)]
  study_ids <- clean_id(study_ids)
  dt <- dt[Dataset == "Reference" | IID %chin% study_ids]
}

dt <- merge(dt, meta[, .(ID, Group)], by.x = "IID", by.y = "ID", all.x = TRUE)
dt[Dataset == "Study", Group := "Study"]

ref_groups <- sort(unique(na.omit(dt[Dataset == "Reference", Group])))
pal <- get_big_palette(length(ref_groups))
names(pal) <- ref_groups

xlim <- range(dt$PC1, na.rm = TRUE)
ylim <- range(dt$PC2, na.rm = TRUE)

# ----------------------------------------------------------------------
# Reference-only plot
# ----------------------------------------------------------------------
if (mode == "reference_only") {
  p_ref <- ggplot(dt[Dataset == "Reference"], aes(PC1, PC2, color = Group)) +
    geom_point(size = 2.5, alpha = 0.95) +
    scale_color_manual(values = pal, drop = FALSE) +
    theme_classic(base_size = 20) +
    theme(legend.position = "bottom") +
    labs(
      title = "Reference PCA",
      x = "PC1",
      y = "PC2",
      color = "Group"
    )

  ggsave("plot_PCA/PCA_PC1_vs_PC2_reference_only.png", p_ref, width = 12, height = 10, dpi = 300)
}

# ----------------------------------------------------------------------
# Projection plots + kNN
# ----------------------------------------------------------------------
if (mode == "projection") {

  p_overlay <- ggplot() +
    geom_point(data = dt[Dataset == "Reference"], aes(PC1, PC2, color = Group), size = 2.5, alpha = 0.95) +
    geom_point(data = dt[Dataset == "Study"], aes(PC1, PC2), color = "black", size = 2.2, alpha = 0.20) +
    scale_color_manual(values = pal, drop = FALSE) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_classic(base_size = 20) +
    theme(legend.position = "bottom") +
    labs(
      title = "PCA: Reference FlashPCA2 + Study projection",
      x = "PC1",
      y = "PC2",
      color = "Group"
    )

  ggsave("plot_PCA/PCA_PC1_vs_PC2_reference_plus_study.png", p_overlay, width = 12, height = 10, dpi = 300)

  ref_only <- dt[Dataset == "Reference"][, Facet := "Reference only"]
  ref_plus <- copy(dt)[, Facet := "Reference + Study"]
  df_plot <- rbindlist(list(ref_only, ref_plus), use.names = TRUE, fill = TRUE)
  df_plot[, Facet := factor(Facet, levels = c("Reference only", "Reference + Study"))]

  p_facets <- ggplot() +
    geom_point(data = df_plot[Facet == "Reference only"], aes(PC1, PC2, color = Group), size = 2.3, alpha = 0.95) +
    geom_point(data = df_plot[Facet == "Reference + Study" & Dataset == "Reference"], aes(PC1, PC2, color = Group), size = 2.3, alpha = 0.95) +
    geom_point(data = df_plot[Facet == "Reference + Study" & Dataset == "Study"], aes(PC1, PC2), color = "black", size = 2.0, alpha = 0.20) +
    scale_color_manual(values = pal, drop = FALSE) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_classic(base_size = 18) +
    theme(legend.position = "bottom") +
    labs(
      title = "PCA facets: Reference only vs Reference + Study",
      x = "PC1",
      y = "PC2",
      color = "Group"
    ) +
    facet_wrap(~ Facet, nrow = 1)

  ggsave("plot_PCA/PCA_PC1_vs_PC2_facets.png", p_facets, width = 16, height = 8, dpi = 300)

  ref <- dt[Dataset == "Reference"]
  qry <- dt[Dataset == "Study"]

  if (nrow(ref) > 0 && nrow(qry) > 0) {

    ref_mat <- as.matrix(ref[, ..use_pcs])
    qry_mat <- as.matrix(qry[, ..use_pcs])

    k <- min(k, nrow(ref))
    threshold_N <- min(threshold_N, k)

    nn <- FNN::get.knnx(data = ref_mat, query = qry_mat, k = k)

    ref_levels <- sort(unique(na.omit(ref$Group)))
    lbl_mat <- matrix(
      as.character(ref$Group)[as.vector(nn$nn.index)],
      nrow = nrow(nn$nn.index),
      byrow = FALSE
    )

    eps <- 1e-9
    w <- 1 / (nn$nn.dist + eps)

    counts_mat <- sapply(ref_levels, function(g) rowSums(lbl_mat == g, na.rm = TRUE))
    counts_mat <- as.matrix(counts_mat)
    colnames(counts_mat) <- paste0(ref_levels, "_count")

    wsum_mat <- sapply(ref_levels, function(g) rowSums((lbl_mat == g) * w, na.rm = TRUE))
    wsum_mat <- as.matrix(wsum_mat)
    colnames(wsum_mat) <- paste0(ref_levels, "_wsum")

    pred <- character(nrow(qry))

    for (i in seq_len(nrow(qry))) {
      row_counts <- counts_mat[i, , drop = TRUE]
      hits <- sub("_count$", "", names(row_counts))[row_counts >= threshold_N]

      if (length(hits) == 0L) {
        pred[i] <- "Undefined"
      } else if (length(hits) == 1L) {
        pred[i] <- hits
      } else {
        ws <- wsum_mat[i, paste0(hits, "_wsum"), drop = TRUE]
        hits2 <- hits[ws == max(ws)]
        pred[i] <- if (length(hits2) == 1L) hits2 else "Undefined"
      }
    }

    counts_dt <- as.data.table(counts_mat)
    props_dt <- counts_dt[, lapply(.SD, function(x) x / k)]
    setnames(props_dt, names(props_dt), sub("_count$", "_prop", names(props_dt)))

    result <- cbind(
      qry[, .(FID, IID, PC1, PC2)],
      Predicted_Group = pred,
      counts_dt,
      props_dt
    )

    out_tsv <- sprintf("plot_PCA/Study_KNN_breakdown_k%d_thresh%d_PC%d.tsv", k, threshold_N, n_pcs)
    fwrite(result, out_tsv, sep = "\t")

    all_groups <- sort(unique(c(ref_groups, "Undefined")))
    pal_cag <- c(pal, "Undefined" = "grey70")[all_groups]

    qry_plot <- merge(qry, result[, .(IID, Predicted_Group)], by = "IID", all.x = TRUE)
    qry_plot[is.na(Predicted_Group), Predicted_Group := "Undefined"]

    p_study <- ggplot(qry_plot, aes(PC1, PC2, color = Predicted_Group)) +
      geom_point(size = 2.3, alpha = 0.95) +
      scale_color_manual(values = pal_cag, drop = TRUE) +
      coord_cartesian(xlim = xlim, ylim = ylim) +
      theme_classic(base_size = 20) +
      theme(legend.position = "bottom") +
      labs(
        title = sprintf("Study only — inferred group (k=%d, PCs=1-%d, thresh=%d)", k, n_pcs, threshold_N),
        x = "PC1",
        y = "PC2",
        color = "Inferred group"
      )

    ggsave("plot_PCA/PCA_Study_by_inferred_group.png", p_study, width = 12, height = 10, dpi = 300)

    max_scatter <- min(5L, n_pcs)
    sm_cols <- paste0("PC", seq_len(max_scatter))
    cols_map <- setNames(c(pal, "Undefined" = "grey70"), c(ref_groups, "Undefined"))
    pt_col <- cols_map[as.character(qry_plot$Predicted_Group)]

    png("plot_PCA/PCA_Study_scattermatrix_PC1toPC5.png", width = 1200, height = 1200, res = 120)
    pairs(qry_plot[, ..sm_cols], col = pt_col, pch = 16, upper.panel = NULL, diag.panel = NULL, cex = 0.6)
    legend("topright", inset = 0.02, legend = names(cols_map), col = cols_map, pch = 16, cex = 0.9, bty = "n")
    dev.off()

  } else {
    message("Skipped kNN: need both Reference and Study samples present.")
  }
}
