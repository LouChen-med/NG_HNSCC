library(Seurat)
library(tidyverse)
library(patchwork)
library(scalop)
library(here)
source(here("analysis_scripts/functions/scRNA_Functions.R"))
source(here("analysis_scripts/functions/ST_Functions.R"))
source(here("analysis_scripts/functions/VisiumHD_Functions.R"))


# Load Metadata and Signatures ---------------------------------------------------------------

# Load per-sample parameter table
sample_parameters <- readRDS(file = here("metadata/per_sample_parameters_df.rds"))
# Load sample level metadata
samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))
# Load single-cell cell-type and cellular-states programs
sc_progs <- xlsx::read.xlsx(here("aux_data/HNSCC_Cancer_and_TME_Progs.xlsx"), sheetIndex = 1) %>%
  dplyr::select(-c("G1S", "G2M"))
sc_progs <- lapply(sc_progs, function(prog) {
  prog[!is.na(prog)]
})
names(sc_progs)[names(sc_progs) == "Immunosuppr_Macrophage"] <- "Imm_supp_Mac"
sc_progs$B_cell <- c(sc_progs$B_cell, "MZB1", "DERL3", "XBP1", "JCHAIN", "SDC1", "PRDM1", "IGHG1")
sc_progs$Plasma_cell <- NULL
# Load Epithelial vs Stromal signatures
epi_stroma_progs <- readRDS(file = here("results/Generated_Data/epi_nonEpi_gene_programs.rds"))
# Load state colors
state_cols <- read.delim(file = here("aux_data/state_col_tab_extend.tsv"), sep = "\t", header = TRUE)
state_cols <- setNames(state_cols$V2, state_cols$V1)
state_cols[names(state_cols) == "LowQ"] <- "grey95"
# Create a complete color palette
cell_cols <- c(state_cols[names(state_cols) %in% names(sc_progs)],
               setNames(c("#999933", "#117733", "#008EA0FF", "#44AA99", "#703092", "#FF2400", "#DD1C77", "grey95", "transparent"),
                        c("Epi.diff", "Secretory", "Stress", "Myocyte", "Plasma_cell", "Imm_supp_Mac", "Neutrophil", "LowQ", "Unresolved")))

state_cols <- c(state_cols, cell_cols[names(cell_cols) %in% c("Imm_supp_Mac", "Neutrophil", "Unresolved")])


# Score and Plot Per-Sample Visium HD Spots --------------------------------------------------

# Restrict to the Visium HD samples used in the analysis
sample_parameters <- sample_parameters %>% dplyr::filter(Sample %in% c("P5624", "P5684"))
visiumHD_dir <- here("data", "VisiumHD")

for(i in seq_along(sample_parameters$Sample)) {
  message(paste0("Processing sample: ", sample_parameters$Sample[[i]]))

  # Load expression matrix
  samp_dir <- grep(sample_parameters$Sample[[i]], list.dirs(visiumHD_dir, recursive = FALSE), value = TRUE)
  data_dir <- paste0(samp_dir, "/outs/binned_outputs/square_016um")
  m <- as.matrix(Matrix::readMM(paste(data_dir, grep("matrix(_[A-Za-z0-9]+)?\\.mtx$", list.files(data_dir), value = TRUE), sep = "/")))
  colnames(m) <- read.table(paste(data_dir, grep("barcodes(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 1]
  rownames(m) <- read.table(paste(data_dir, grep("features(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 2]
  if (sum(grep("GRCh38", rownames(m))) == 0) {
    m <- m
  } else {
    # Remove HPV genes and GRCh38 prefix
    m <- m[-grep("^hpv", rownames(m)), ]
    rownames(m) <- sub("GRCh38_", "", rownames(m))
  }
  colnames(m) <- paste(sample_parameters$Sample[[i]], colnames(m), sep = "_")

  # Load spatial image object
  spatial_image <- load_spatial_image(data_dir, image_res = "lowres")

  # Process raw data (no centering)
  proc_mat <- m[-grep("^MT-|^RPL|^RPS", rownames(m)), ]
  if(min(colSums(proc_mat)) == 0) {proc_mat <- proc_mat[, colSums(proc_mat) > 0]}
  scaling_factor <- 1000000/colSums(proc_mat)
  proc_mat <- sweep(proc_mat, MARGIN = 2, STATS = scaling_factor, FUN = "*")
  proc_mat <- log2(1 + (proc_mat/10))
  # proc_mat <- proc_mat[rowMeans(proc_mat) > 0.1, ]
  # top_expr_genes <- proc_mat %>% rowMeans() %>% sort(decreasing = T) %>% head(7500) %>% names
  # proc_mat <- proc_mat[top_expr_genes, ]

  # Score for Epithelial vs Stromal signatures and Metaprogram identity and plot per-sample spatial heatmap
  # epi_stroma_score <- metaprog_score(proc_mat, epi_stroma_progs, score_diff = 0.4, score_cut = 0.5, bin_size = 100, bin_num = 30, conserved.genes = 0, center_rows = TRUE, scale = FALSE)
  # spot_type_scores <- metaprog_score(proc_mat, sc_progs, score_diff = 0, score_cut = 0, bin_size = 100, bin_num = 30, conserved.genes = 0, center_rows = TRUE, scale = FALSE, n_cores = 8)
  # spot_type_scores <- metaprog_score(proc_mat, sc_progs, score_diff = 0, score_cut = 0, bin_size = 100, bin_num = 30, conserved.genes = 0.5, center_rows = TRUE, scale = TRUE, n_cores = 8)
  # setdiff(names(sc_progs), unique(unname(spot_type_scores$meta_program)))
  spot_type_scores <- scalop::sigScores(m = proc_mat, sigs = sc_progs, conserved.genes = 0.5)

  # Assign top scoring cell-type/state
  mat <- as.matrix(spot_type_scores)
  mat[is.na(mat)] <- -Inf
  idx <- max.col(mat, ties.method = "first")
  row_max <- apply(mat, 1, max)
  assign_vec <- colnames(mat)[idx]
  assign_vec[row_max <= 0 | is.infinite(row_max)] <- "Unresolved"
  names(assign_vec) <- rownames(mat)

  cell_stats <- cell_qc(m) %>% dplyr::filter(CellID %in% colnames(proc_mat))
  samp_meta <- tibble(SpotID = stringr::str_remove(colnames(proc_mat), "^[^_]*_"),
                      Sample = stringr::str_extract(colnames(proc_mat), "^[^_]+"),
                      Key = colnames(proc_mat),
                      Complexity = unname(cell_stats$Complexity),
                      Percent_Mito = unname(cell_stats$Percent_mt),
                      Percent_HK = unname(cell_stats$HK),
                      # Type_State = unname(spot_type_scores$meta_program[match(colnames(proc_mat), names(spot_type_scores$meta_program))]),
                      Type_State = unname(assign_vec[match(Key, names(assign_vec))]),
                      Inclusion = ifelse(Complexity > 150 & Percent_Mito < 20, "Include", "Exclude"))

  # Plot QC parameters
  p <- ggpubr::ggviolin(samp_meta$Complexity, fill = "lightskyblue", xlab = "Complexity", ylab = FALSE) + ggplot2::theme(axis.text.x = element_blank()) +
    ggpubr::ggviolin(samp_meta$Percent_Mito, fill = "lightskyblue", xlab = "Percent_Mito", ylab = FALSE) + ggplot2::theme(axis.text.x = element_blank()) +
    ggpubr::ggviolin(samp_meta$Percent_HK, fill = "lightskyblue", xlab = "Percent_HouseKeeping", ylab = FALSE) + ggplot2::theme(axis.text.x = element_blank())
  ggplot2::ggsave(filename = paste0(here("results/Per_Sample_Results/"), sample_parameters$Sample[[i]], "/", sample_parameters$Sample[[i]], "_VisiumHD_Per-Spot_QC_Parameters.png"),
                  plot = p, height = 8, width = 10)

  p <- ggpubr::gghistogram(samp_meta$Complexity, xlab = "Number of detected genes per spot", ylab = "Number of spots", fill = "lightsteelblue3",
                           alpha = 0.4, bins = 40) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) + ggplot2::scale_x_continuous(breaks = round(seq(0, max(samp_meta$Complexity), by = 200), 1), limits = c(0, max(samp_meta$Complexity)), expand = c(0, 0)) +
    ggplot2::geom_vline(xintercept = 150, color = "red") +
    ggplot2::theme(panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"), aspect.ratio = 0.7)
  ggplot2::ggsave(filename = paste0(here("results/Per_Sample_Results/"), sample_parameters$Sample[[i]], "/", sample_parameters$Sample[[i]], "_VisiumHD_Complexity_Distribution.png"),
                  plot = p, height = 8, width = 10)

  # Plot spatial map
  samp_meta$Plot <- ifelse(samp_meta$Inclusion == "Include", samp_meta$Type_State, "LowQ")
  p <- plot_spatial_features(metadata = samp_meta, image_obj = spatial_image, color_by = "Plot", cols = cell_cols[names(cell_cols) %in% unique(samp_meta$Plot)], image.alpha = 0, pt.size = 4) +
    guides(fill = guide_legend(ncol = 1, override.aes = list(size = 5)))
  ggplot2::ggsave(filename = paste0(here("results/Paper_Figures/Revision_Figs/"), sample_parameters$Sample[[i]], "_VisiumHD_Dominant_scRNAseq_State_Map.png"),
                  plot = p, width = 35, height = 35, units = "cm", dpi = 600)

  samp_meta <- samp_meta %>% dplyr::mutate(Plot = dplyr::case_when(Complexity > 150 & Percent_Mito < 20 & Type_State %in% c("Stress", "Senescence", "Hypoxia", "Epi.diff", "Cell_Cycle") ~ "non-pEMT_Malig",
                                                                  Complexity > 150 & Percent_Mito < 20 & Type_State == "pEMT" ~ "pEMT",
                                                                  Complexity > 150 & Percent_Mito < 20 & Type_State %in% c("Imm_supp_Mac", "Neutrophil") ~ "Pro-tumor myeloid",
                                                                  TRUE ~ "none")) %>%
    dplyr::mutate(Plot = factor(Plot, levels = c("pEMT", "Pro-tumor myeloid", "non-pEMT_Malig", "none"))) %>% droplevels()
  tmp_cols <- setNames(c("#FF2400", "#88CCEE", "#999933", "grey95"), c("pEMT", "Pro-tumor myeloid", "non-pEMT_Malig", "none"))
  p <- plot_spatial_features(metadata = samp_meta, image_obj = spatial_image, color_by = "Plot", cols = tmp_cols[names(tmp_cols) %in% unique(samp_meta$Plot)], pt.size = 4, image.alpha = 0, rm_spots = NULL, na.value = "gray", stroke = 0) +
    guides(fill = guide_legend(override.aes = list(size = 5), title = "State-pair"))
  ggplot2::ggsave(filename = paste0(here("results/Paper_Figures/Revision_Figs/"), sample_parameters$Sample[[i]], "_VisiumHD_pEMT_ISM_Neutrophil_scRNAseq_progs.png"),
                  plot = p, width = 35, height = 35, units = "cm", dpi = 600)

  samp_meta <- samp_meta %>% dplyr::mutate(Plot = dplyr::case_when(Complexity > 150 & Percent_Mito < 20 & Type_State %in% c("Stress", "Senescence", "Hypoxia", "Epi.diff", "Cell_Cycle") ~ "non-pEMT_Malig",
                                                                  Complexity > 150 & Percent_Mito < 20 & Type_State %in% c("pEMT", "Imm_supp_Mac", "Neutrophil") ~ as.character(Type_State), .default = "none")) %>%
    dplyr::mutate(Plot = factor(Plot, levels = c("pEMT", "Imm_supp_Mac", "Neutrophil", "non-pEMT_Malig", "none"))) %>% droplevels()
  tmp_cols <- setNames(c("#053061", "#FF2400", "#DD1C77", "#d1e5f0", "grey95"), c("pEMT", "Imm_supp_Mac", "Neutrophil", "non-pEMT_Malig", "none"))
  p <- plot_spatial_features(metadata = samp_meta, image_obj = spatial_image, color_by = "Plot", cols = tmp_cols[names(tmp_cols) %in% unique(samp_meta$Plot)], pt.size = 4, image.alpha = 0, rm_spots = NULL, na.value = "gray", stroke = 0) +
    guides(fill = guide_legend(override.aes = list(size = 5), title = "State-pair"))

  microns_per_pixel <- dplyr::case_when(sample_parameters$Sample[[i]] == "P5624" ~ 0.5550283554379033,
                                        sample_parameters$Sample[[i]] == "P5684" ~ 0.4914280829307512)
  p <- add_spatial_scale_bar(p = p, image_obj = spatial_image, microns_per_pixel = microns_per_pixel, bar_um = 1000)
  ggplot2::ggsave(filename = paste0(here("results/Paper_Figures/Revision_Figs/"), sample_parameters$Sample[[i]], "_VisiumHD_pEMT_ISM_Neutrophil_scRNAseq_progs_V2.png"),
                  plot = p, width = 35, height = 35, units = "cm", dpi = 600)
}
