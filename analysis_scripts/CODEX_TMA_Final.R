library(here)
library(tidyverse)
library(Seurat)
library(pheatmap)
library(ggforce)
library(mclust)
library(ComplexHeatmap)
library(sf)
library(randomForest)
library(FNN)
source(here("analysis_scripts/functions/CODEX_Functions.R"))

# Load sample level metadata
samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))
# Load state colors
state_cols <- read.delim(file = here("aux_data/state_col_tab_extend.tsv"), sep = "\t", header = TRUE)
state_cols <- setNames(state_cols$V2, state_cols$V1)


flip_map <- c(A2 = TRUE, A3 = TRUE, A4 = FALSE, A5 = TRUE,
              B2 = FALSE, B5 = TRUE, B6 = TRUE, B8 = TRUE, B10 = TRUE,
              C2 = TRUE, C3 = TRUE, C4 = TRUE, C6 = TRUE, C7 = TRUE, C9 = TRUE, C10 = TRUE,
              D2 = TRUE, D3 = TRUE, D4 = TRUE, D5 = FALSE, D6 = FALSE, D7 = FALSE, D8 = FALSE, D9 = FALSE, D10 = FALSE,
              E3 = TRUE, E5 = TRUE)


# Load HNSCC CODEX TMA ----------------------------------------------------------

# Load InstanSeg results
Instanseg_final <- read_csv(here("data", "CODEX", "Instanseg_final.csv"))
# Instanseg_final <- read_csv(here("Analysis/Revision_Figs/CODEX_TMA/InstanSeg_Full.csv"))
colnames(Instanseg_final) <- gsub("^Cell\\.\\.([^.]+)\\.\\.Mean$", "\\1", colnames(Instanseg_final))
# Load Nimbus scores
Nimbus <- read_csv(here("data", "CODEX", "nimbus_cell_table_all.csv"))
# Nimbus <- read_csv(here("Analysis/Revision_Figs/CODEX_TMA/Nimbus_Full.csv"))
# Load marker info table
marker_info <- read.csv(file = here("aux_data/marker_info.csv"))

# Prepare data for downstream analysis
rmat <- rbind(Nimbus %>% mutate(label = paste0(fov, "_", label)))
cells <- as.data.frame(Instanseg_final) %>% 
  dplyr::select(cell_name = Image,
                image_id = Image_1,
                centroid_x = `Centroid.X.µm`,
                centroid_y = `Centroid.Y.µm`,
                cell_id = `Object.ID`,
                classification = `Classification`)
cells <- dplyr::left_join(cells, rmat, by = c("cell_name" = "label"))
# Remove undesired tissue regions (folds / low quality) - will present as NA in the FOV column.
cells_all <- cells[!is.na(cells$fov), ]

# Create matrix data 
rmat <- cells_all %>% 
  dplyr::select(-cell_name, -image_id, -centroid_x, -centroid_y, -cell_id, -fov, -classification) %>% 
  as.matrix() %>% t()
colnames(rmat) <- cells_all$cell_name

# Create annotation data
cells <- cells_all%>% 
  dplyr::select(cell_name, image_id, centroid_x, centroid_y, cell_id)



# QC and Batch Assessment -------------------------------------------------

# Number of cells per TMA core
table(sub("^((?:fov0_[A-Za-z][0-9])).*$", "\\1", colnames(rmat)))
core_id <- sub("^fov0_([A-Za-z0-9]+)_.*$", "\\1", colnames(rmat)) %>% factor()

### ============ QC histograms ============
plot_df <- as.data.frame(rmat) %>% 
  rownames_to_column(var = "Markers") %>% 
  pivot_longer(cols = -Markers, names_to = "Cell", values_to = "Scores")
valleys <- plot_df %>% 
  group_by(Markers) %>% 
  summarize(valley = {
      d <- density(Scores, na.rm = TRUE)
      x <- d$x; y <- d$y
      peak_inds <- which(diff(sign(diff(y))) < 0) + 1
      if(length(peak_inds) < 2) return(NA_real_)
      
      top2 <- peak_inds[order(y[peak_inds], decreasing = TRUE)[1:2]]
      lo <- min(x[top2]); hi <- max(x[top2])
      between <- which(x > lo & x < hi)
      if(length(between) == 0) return(NA_real_)
      
      x[ between[ which.min(y[between]) ] ]
    }, .groups="drop")
ncol <- 4
nrow <- 4
perPage <- ncol * nrow
nPages <- ceiling(n_distinct(plot_df$Markers) / perPage)

base <- ggplot(plot_df, aes(x = Scores)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Markers Nimbus Scores", x = "Score", y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5, size = 14), strip.text = element_text(size = 8))

for(page in seq_len(nPages)) {
  p <- base +
    facet_wrap_paginate(~ Markers, ncol = ncol, nrow = nrow, scales = "free_y", page = page) +
    geom_vline(data = valleys, aes(xintercept = valley), color = "red", linetype = "dashed") +
    geom_text(data = valleys, aes(x = valley, y = Inf, label = round(valley, 2)), color = "red", hjust = 1.1, vjust = 2, size = 3)
  
  ggsave(filename = paste0(here("results/Diagnostic_Plots/CODEX_TMA/"), sprintf("marker_histograms_page%02d.png", page)), plot = p, width = 12, height = 12)
  # ggsave(filename = paste0(here("Analysis/Revision_Figs/CODEX_TMA/extend2/"), sprintf("marker_histograms_page%02d.png", page)), plot = p, width = 12, height = 12)
}


### ============ Marker correlation heatmap and clustering ============
marker_cor <- cor(t(rmat), use = "pairwise.complete.obs", method = "pearson")
m_dist <- as.dist(1 - (marker_cor))
m_hc  <- hclust(m_dist, method = "ward.D2")
plot(m_hc)
rect.hclust(m_hc, k = 12, border = 2:6)
p <- ComplexHeatmap::Heatmap(marker_cor, 
                             name = "pearson",
                             col = circlize::colorRamp2(breaks = seq(-1, 1, length.out = 11), hotmap),
                             row_names_side = "left",
                             row_dend_side = "right",
                             row_title_side = "right",
                             show_column_names = T,
                             show_row_names = T,
                             column_split = 12,
                             row_split = 12,
                             cluster_rows = make_hclust(),
                             cluster_columns = make_hclust(),
)
pdf(file = here("results/Diagnostic_Plots/CODEX_TMA/Marker_Correlation_Heatmap_Unfiltered.pdf"), width = 16, height = 14)
# pdf(file = here("Analysis/Revision_Figs/CODEX_TMA/extend2/Marker_Correlation_Heatmap_Unfiltered.pdf"), width = 16, height = 14)
draw(p)
dev.off()


### ============ Assess batch effect between the different cores ============
core_id <- sub("^fov0_([A-Za-z0-9]+)_.*$", "\\1", colnames(rmat)) %>% factor()
metadata <- data.frame(cell_id = colnames(rmat),
                       core_id = core_id,
                       stringsAsFactors = FALSE)

dens_df <- as.data.frame(t(rmat))
dens_df <- dens_df %>% 
  mutate(core_id = metadata$core_id) %>% 
  pivot_longer(cols = -core_id, names_to  = "marker", values_to = "intensity")

ncol <- 4
nrow <- 4
perPage <- ncol * nrow
nPages <- ceiling(n_distinct(dens_df$marker) / perPage)
base <- ggplot(dens_df, aes(intensity, colour = core_id, group = core_id)) +
  geom_density(size = .3, alpha = .25, show.legend = FALSE) +
  labs(title = "Raw intensity densities by core",
       x = "intensity", y = "density") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(hjust = .5))

for (page in seq_len(nPages)) {
  p <- base +
    facet_wrap_paginate(~ marker, ncol = ncol, nrow = nrow, scales = "free", page = page)
  
  ggsave(filename = here("results/Diagnostic_Plots/CODEX_TMA/", sprintf("per-core_marker_density_page%02d.png", page)), plot = p, width = 12, height = 12, dpi = 300)
  # ggsave(filename = here("Analysis/Revision_Figs/CODEX_TMA/extend2/", sprintf("per-core_marker_density_page%02d.png", page)), plot = p, width = 12, height = 12, dpi = 300)
}




# Phenotype Cell Types in a Supervised Manner (gating) ----------------------------------

# Remove DAPI stain and undesired markers
markers_to_rm <- c("CD138", "HLAA", "cMyc", "SPRR1B", "SPRR1Bv2", "DAPI", "TIMP1", "PDL1","IL1B", "ZEB1", "Col4a1", 
              "TIMP1v2", "ISG20", "FGB", "CLDN4", "CD3", "CTLA4", "Tbet", "PD1", "CD19v2", "CD4", "CD11c", "LAG3", 
              "CD45", "TFF1", "LRRC15","Caspase3")
rmat <- rmat[!rownames(rmat) %in% markers_to_rm, ]

### ======= Programmatically derive thersholds for gating using GMM with assignment confidence (posterior probabilities) =======
valid_markers <- rownames(rmat)

# FIRST, assess the bimodality and separation of peaks of each of the cell-type markers
for (marker in valid_markers) {
  # Pull out the numeric vector of scores for this marker
  scores <- as.numeric(rmat[marker, ])
  cell_ids <- colnames(rmat)
  
  # Check bimodality using BIC - which has a better fit 1 or 2 component Gaussians and Ashman's (which requires separation - guarding against a tiny shoulder)
  mod_sel <- mclust::Mclust(scores, G = 1:2, modelNames = "V")
  is_bimodal <- (mod_sel$G == 2)
  
  if (is_bimodal) {
    # check separation: Ashman's D > 2 ⇒ well-separated peaks
    mu <- mod_sel$parameters$mean
    sig2 <- mod_sel$parameters$variance$sigmasq
    ashD <- abs(mu[1] - mu[2]) / sqrt((sig2[1] + sig2[2]) / 2)
  }
  print(paste0(marker, ": ", ashD))
}


### ============ Compute thresholds and posterior probabilities for each marker ============
thresholds_list  <- list()
posteriors_list <- list()

for (marker in valid_markers) {
  # 1a. extract scores & cell IDs
  scores   <- as.numeric(rmat[marker, ])
  cell_ids <- colnames(rmat)
  
  # 1b. fit GMMs: 1- vs 2-components
  fit1 <- Mclust(scores, G = 1, modelNames = "V", verbose = FALSE)
  fit2 <- Mclust(scores, G = 2, modelNames = "V", verbose = FALSE)
  
  # 1c. compute AICs
  ll1   <- fit1$loglik; ll2 <- fit2$loglik
  df1   <- 3*1 - 1;      df2 <- 3*2 - 1  # params = 3G-1
  aic1  <- 2*df1 - 2*ll1
  aic2  <- 2*df2 - 2*ll2
  
  # 1d. Ashman’s D
  mu    <- fit2$parameters$mean                # length-2
  sig2   <- fit2$parameters$variance$sigmasq   # length-2
  ashD  <- abs(diff(mu)) / sqrt(sum(sig2)/2)
  
  # 1e. bimodality test
  bim   <- (aic2 < aic1) && (ashD > 2)
  
  # 1f. order components by mean: low = neg, high = pos
  ord        <- order(mu)
  neg_comp   <- mu[ord[1]]
  pos_comp   <- mu[ord[2]]
  props      <- fit2$parameters$pro            # mixing weights
  weight_neg <- props[ord[1]]
  weight_pos <- props[ord[2]]
  
  # 1g. per-cell posteriors
  #    (use the same ord[] to pick which component is “positive”)
  p_pos <- fit2$z[, ord[2]]
  p_neg <- fit2$z[, ord[1]]
  
  # 1h. collect into our lists
  thresholds_list[[marker]] <- tibble(
    marker      = marker,
    AIC1        = aic1,
    AIC2        = aic2,
    AshD        = ashD,
    bimodal     = bim,
    neg_comp    = neg_comp,
    pos_comp    = pos_comp,
    weight_neg  = weight_neg,
    weight_pos  = weight_pos
  )
  
  posteriors_list[[marker]] <- tibble(
    cell_id     = cell_ids,
    marker      = marker,
    p_pos       = p_pos,
    p_neg       = p_neg
  )
}

thresholds_df  <- bind_rows(thresholds_list)
posteriors_df  <- bind_rows(posteriors_list)
saveRDS(thresholds_df, file = here("results/Generated_Data/CODEX/markers_thresholds_df.rds"))
saveRDS(posteriors_df, file = here("results/Generated_Data/CODEX/markers_posterior_prob_df.rds"))


### ============ Gate cells based on their probability being cell-type marker positive ============
thresholds_df  <- readRDS(file = here("results/Generated_Data/CODEX/markers_thresholds_df.rds"))
posteriors_df  <- readRDS(file = here("results/Generated_Data/CODEX/markers_posterior_prob_df.rds"))
broad_cut <- 0.5
high_cut <- 0.9

p_wide <- posteriors_df %>%
  dplyr::select(cell_id, marker, p_pos) %>%
  tidyr::pivot_wider(names_from = marker,
                     values_from = p_pos) %>% 
  dplyr::mutate(CellType_broad = NA_character_,
                CellType_high_conf = NA_character_)

# Define helper to apply one gate
apply_gate <- function(df, cut, outcol) {
  df %>%
    mutate(!!outcol := case_when(
      # Epithelium: ANY of PanCK/EGFR/TP63 above cut
      (PanCK >= cut | EGFR >= cut | TP63 >= cut) ~ "Epithelial",
      
      # Fibroblasts: aSMA + PDGFRAB
      (aSMA >= cut & PDGFRAB >= cut)             ~ "Fibroblast",
      
      # Myeloid: both CD68 and CD11b, or CD15 above cut
      # (CD68 >= cut & CD11b >= cut | CD15 >= cut) ~ "Myeloid",
      (CD68 >= cut | CD11b >= cut | CD15 >= cut) ~ "Myeloid",
      
      # Endothelial: CD31
      (CD31 >= cut)                              ~ "Endothelial",
      
      # T cytotoxic cells: CD8 + GranzymeB
      (CD8   >= cut & GranzymeB >= cut)          ~ "Tctx-cell",
      
      # B-cells: CD19
      (CD19 >= cut)                              ~ "B-cell",
      
      # 3.7 fallback
      TRUE                                       ~ "Undefined"
    ))
}

# Apply both gates
p_wide <- p_wide %>%
  apply_gate(broad_cut, "CellType_broad") %>%
  apply_gate(high_cut, "CellType_high_conf")

cell_assign <- p_wide %>%
  transmute(CellID = cell_id, cell_type = CellType_high_conf)
order_markers <- unique(c(marker_info$Marker[marker_info$Cell_Type == "Epithelial"],
                          marker_info$Marker[marker_info$Cell_Type == "Fibroblast"],
                          marker_info$Marker[marker_info$Cell_Type %in% c("Endothelial", "Lymphatic_vessel")],
                          marker_info$Marker[marker_info$Cell_Type == "Myeloid"],
                          marker_info$Marker[marker_info$Cell_Type == "T_cell"],
                          marker_info$Marker[marker_info$Cell_Type == "B_cell"],
                          "ATP1A1", "LAMC2", "Podoplanin", "Ki67", "CA9", "GLUT1", "pNRP2", "pH2Ax", "Cytochromec", "SLPI", "S100A9",
                          "MMP9", "FN1", "Perlecan", "SPARC", "C3", "HAPLN1",
                          "cMyc", "PDL1", "IL1B", "SPRR1Bv2", "CD163", "MMP14", "VCAN", "TIMP1", "Col1a1", "PDGFRAB", "Col4a1", "Vimentin", 
                          "ZEB1", "LRRC15", "TFF1", "CTLA4", "CD4", "CD45", "CD3", "CLDN4", "HLAA", "LAG3", "CD11c", "PD1", "Tbet", 
                          "TNC", "MMP2", "FGBv2", "HA", "CD56", "Tryptase", "Osteoponin",
                          "DAPI", "ISG20", "p16", "Caspase3", "ATF3", "p21", "IFNG", "CD19v2", "CD45v2", "FGB", "SPRR1B", "TIMP1v2"))
order_celltype <- c("Epithelial", "Fibroblast", "Endothelial", "Myeloid", "Tctx-cell", "B-cell", "Undefined")
p <- plot_marker_by_celltype(
  rmat,
  cell_assignments = cell_assign,
  expression_threshold = 0.3,
  marker_order = order_markers[order_markers %in% rownames(rmat)],
  celltype_order = order_celltype,
  cluster_rows = FALSE,
  cluster_cols = FALSE)
ggsave(filename = here("results/Diagnostic_Plots/CODEX_TMA/GMM_Cell_Type_Marker_Scores_Dotplot.png"), plot = p, width = 22, height = 10)
# ggsave(filename = here("Analysis/Revision_Figs/CODEX_TMA/extend2/GMM_Cell_Type_Marker_Scores_Dotplot.png"), plot = p, width = 22, height = 10)

p_wide %>% 
  dplyr::count(CellType_high_conf) %>% 
  dplyr::arrange(desc(n)) %>% 
  dplyr::mutate(frac = n / sum(n))

ggsave(filename = here("results/Diagnostic_Plots/CODEX_TMA/GMM_Cell_Type_Marker_Scores_Dotplot_Not_Ordered.png"), plot = p, width = 22, height = 10)


### ============ Deal with cells that were initially undefined by gating ============
# try semi-supervised classification using Random Forest in the assigned cells (by the gating)
library(randomForest)
library(FNN)

# extract marker names
marker_cols <- setdiff(colnames(p_wide),
                       c("cell_id","CellType_broad","CellType_high_conf"))

# define training & undefined sets
train_df <- p_wide %>%
  filter(CellType_high_conf != "Undefined") %>%
  select(cell_id, all_of(marker_cols), cell_type = CellType_high_conf)

undef_df <- p_wide %>%
  filter(CellType_high_conf == "Undefined") %>%
  select(cell_id, all_of(marker_cols))

# ── 1. Semi-supervised RF prediction ──────────────────────────────────────────

# retrain RF on high-conf cells (if not already)
rf <- randomForest(
  x    = train_df %>% select(-cell_id, -cell_type),
  y    = factor(train_df$cell_type),
  ntree= 200
)

# predict class label for undefined cells
rf_probs  <- predict(rf,
                     newdata = undef_df %>% select(-cell_id),
                     type = "prob")
rf_label  <- colnames(rf_probs)[apply(rf_probs,1,which.max)]

# ── 2. Unsupervised clustering + centroid‐matching ──────────────────────────

# compute high_conf centroids of marker-posteriors
centroids <- train_df %>%
  group_by(cell_type) %>%
  summarise(across(all_of(marker_cols), mean), .groups="drop") %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

# k-means cluster undefined cells into K clusters (K = # types)
K  <- nrow(centroids)
set.seed(5410)
km <- kmeans(undef_df %>% select(-cell_id), centers = K, nstart = 10)

# for each cluster, find nearest centroid
cluster_centers <- km$centers          # K × markers
dmat <- as.matrix(dist(rbind(cluster_centers, centroids)))
# rows 1:K = clusters, cols (K+1):(K+T) = true types
cluster_to_type <- apply(
  dmat[1:K, (K+1):(K+nrow(centroids))],
  1, function(r) rownames(centroids)[which.min(r)]
)

cl_label <- cluster_to_type[km$cluster]

# ── 3. Spatial neighborhood voting ───────────────────────────────────────────

# prepare high_conf coords & labels
hc_coords <- cells %>%
  filter(cell_name %in% train_df$cell_id) %>%
  arrange(match(cell_name, train_df$cell_id))

hc_labels <- train_df$cell_type

# prepare undefined coords
undef_coords <- cells %>%
  filter(cell_name %in% undef_df$cell_id)

# find 10 nearest high-conf neighbors for each undefined
knn_out <- get.knnx(
  data  = hc_coords %>% dplyr::select(centroid_x, centroid_y) %>% as.matrix(),
  query = undef_coords %>% dplyr::select(centroid_x, centroid_y) %>% as.matrix(),
  k     = 10
)

# majority vote with threshold 6/10
sp_label <- apply(knn_out$nn.index, 1, function(idxs) {
  tb  <- table(hc_labels[idxs])
  top <- names(tb)[which.max(tb)]
  if (max(tb) >= 6) top else "Undefined"
})

# ── 4. Consensus 2-of-3 ──────────────────────────────────────────────────────

consensus_df <- data.frame(
  cell_id = undef_df$cell_id,
  rf      = rf_label,
  cl      = cl_label,
  sp      = sp_label,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    votes      = list(c(rf,cl,sp)),
    vote_tab   = list(table(votes)),
    maj_label  = names(vote_tab)[which.max(vote_tab)],
    maj_count  = max(vote_tab),
    final_call = ifelse(maj_count >= 2, maj_label, "Undefined")
  ) %>%
  ungroup() %>%
  select(cell_id, final_call)

# ── 5. Merge back into p_wide ────────────────────────────────────────────────

p_wide2 <- p_wide %>%
  left_join(consensus_df, by = "cell_id") %>%
  rename(CellType_tier2 = final_call) %>% 
  mutate(CellType_Final = case_when(CellType_high_conf != "Undefined" ~ CellType_high_conf,
                                    !is.na(CellType_tier2) ~ CellType_tier2,
                                    TRUE ~ "Undefined"))
saveRDS(p_wide2, file = here("results/Generated_Data/CODEX/cell_type_annotations_by_gating_V2.rds"))
# saveRDS(p_wide2, file = here("Analysis/Revision_Figs/CODEX_TMA/extend2/cell_type_annotations_by_gating_V2.rds"))

p_wide2 %>% 
  dplyr::count(CellType_Final) %>% 
  dplyr::arrange(desc(n)) %>% 
  dplyr::mutate(frac = n / sum(n))


### ============ Final cell-type dotplot ============
cell_assign2 <- p_wide2 %>%
  transmute(CellID = cell_id, cell_type = CellType_Final)

order_markers <- unique(c(marker_info$Marker[marker_info$Cell_Type == "Epithelial"],
                          marker_info$Marker[marker_info$Cell_Type == "Fibroblast"],
                          marker_info$Marker[marker_info$Cell_Type %in% c("Endothelial", "Lymphatic_vessel")],
                          marker_info$Marker[marker_info$Cell_Type == "Myeloid"],
                          marker_info$Marker[marker_info$Cell_Type == "T_cell"],
                          marker_info$Marker[marker_info$Cell_Type == "B_cell"]))
order_markers <- order_markers[order_markers != "Foxp3"]
order_celltype <- c("Epithelial", "Fibroblast", "Endothelial", "Myeloid", "Tctx-cell", "B-cell", "Undefined")

p <- plot_marker_by_celltype(
  rmat,
  cell_assignments = cell_assign2,
  expression_threshold = 0.3,
  marker_order = order_markers[order_markers %in% rownames(rmat)],
  celltype_order = order_celltype,
  cluster_rows = FALSE,
  cluster_cols = FALSE)
ggsave(filename = here("results/Diagnostic_Plots/CODEX_TMA/GMM_Cell_Type_Marker_Scores_Dotplot_Final.png"), plot = p, width = 10, height = 8)
# ggsave(filename = here("Analysis/Revision_Figs/CODEX_TMA/extend2/GMM_Cell_Type_Marker_Scores_Dotplot_Final.png"), plot = p, width = 10, height = 8)

cuts_after <- c("TP63", "PDGFRAB", "Lyve1", "CD15", "CD8")
x_lvls <- levels(p$data$Marker)
vpos <- match(cuts_after, x_lvls) + 0.5
p <- p + geom_vline(xintercept = vpos, linetype = "dashed", linewidth = 0.4, colour = "black")
ggsave(filename = here("results/Diagnostic_Plots/CODEX_TMA/Cell_Type_Marker_Scores_Dotplot.pdf"), plot = p, device = "pdf", dpi = 300, width = 10, height = 8)


### ============ Export .csv for visualization in QuPath ============
p_wide2 <- readRDS(file = here("results/Generated_Data/CODEX/cell_type_annotations_by_gating_V2.rds"))
cell_assign2 <- p_wide2 %>%
  transmute(CellID = cell_id, cell_type = CellType_Final)
cells$assignment <- cell_assign2$cell_type[match(cells$cell_name, cell_assign2$CellID)]
write_csv(cells, file = here("results/Generated_Data/CODEX/cell_type_assign.csv"), col_names = TRUE)
# write_csv(cells, file = here("Analysis/Revision_Figs/CODEX_TMA/extend2/cell_type_assign.csv"), col_names = TRUE)

core_id <- sub("^fov0_([A-Za-z0-9]+)_.*$", "\\1", colnames(rmat)) %>% factor()
metadata <- data.frame(cell_id = colnames(rmat),
                       core_id = core_id,
                       cell_type = cell_assign2$cell_type[match(colnames(rmat), cell_assign2$CellID)],
                       stringsAsFactors = FALSE)




# Phenotype Cell States in a Supervised Manner (Gating) ----------------------------------
cells <- read_csv(file = here("results/Generated_Data/CODEX/cell_type_assign.csv"))
p_wide2 <- readRDS(file = here("results/Generated_Data/CODEX/cell_type_annotations_by_gating_V2.rds"))

cell_assign2 <- p_wide2 %>%
  transmute(CellID = cell_id, cell_type = CellType_Final)

core_id <- sub("^fov0_([A-Za-z0-9]+)_.*$", "\\1", colnames(rmat)) %>% factor()
metadata <- data.frame(cell_id = colnames(rmat),
                       core_id = core_id,
                       cell_type = cell_assign2$cell_type[match(colnames(rmat), cell_assign2$CellID)],
                       stringsAsFactors = FALSE)

### ============ Gate Epithelial cells ============
cut <- 1
p_wide2 <- p_wide2 %>%
  mutate(CellState = case_when(CellType_Final == "Epithelial" & LAMC2 >= 1 & p16 < cut ~ "pEMT",
                               CellType_Final == "Epithelial" & ((pH2Ax >= cut | p16 >= cut | p21 >= cut) & (SLPI  >= cut | S100A9 >= cut) & (Ki67 <= 0.9)) ~ "Senescence",
                               CellType_Final == "Epithelial" & CA9 >= cut & GLUT1 >= cut & p16 < cut ~ "Hypoxia",
                               CellType_Final == "Epithelial" & Ki67 >= cut ~ "Cell-Cycle",
                               # CellType_Final == "Epithelial" & Ki67 > quantile(p_wide2$Ki67, probs = 0.25) & SLPI < 1 ~ "Epithelial",
                               CellType_Final == "Epithelial" & SLPI < 1 & p16 < 1 ~ "Epithelial",
                               CellType_Final == "Epithelial" ~ "Senescence",
                               TRUE ~ CellType_Final))

### ============ Gate Myeloid cells into their subtypes ============
# define myeloid gating markers
my_markers <- c("CD15", "S100A9", "CD68", "CD163", "HLADR", "CD11b")
cutoff <- 1    # posterior probability cutoff

p_wide2 <- p_wide2 %>%
  mutate(CellState = case_when(CellType_Final == "Epithelial" ~ CellState,
      CellType_Final == "Myeloid" &
        # (1) Neutrophil = CD15 & S100A9
        CD11b >= cutoff & if_any(all_of(c("CD15", "S100A9")), ~ . >= cutoff)
      ~ "Neutrophil",
      
      CellType_Final == "Myeloid" &
        # (2) M2‐like = CD68 & CD163
        CD163  >= cutoff & if_any(all_of(c("CD68", "CD11b")), ~ . >= cutoff)
      ~ "M2-like_Macrophage",
      
      CellType_Final == "Myeloid" &
        # (3) M1-like = CD68 & HLA-DR & IL-1β
        CD68   >= cutoff & HLADR   >= cutoff
      ~ "M1-like_Macrophage",
      
      CellType_Final == "Myeloid" &
        # (4) Monocyte = CD11b & HLA-DR (CD68 low by exclusion)
        CD11b  >= cutoff & HLADR   >= cutoff
      ~ "Monocyte",
      
      # fallback: leave whatever broad type you had
      TRUE ~ CellType_Final))


cells <- cells %>% dplyr::mutate(assignment = p_wide2$CellState)
write_csv(cells, file = here("results/Generated_Data/CODEX/cell_type_and_state_assign_v2.csv"), col_names = TRUE)


round(prop.table(table(codex_parse_core_id(cells$cell_name), cells$assignment), margin = 1) * 100, 3)

cell_assign3 <- cells %>%
  transmute(CellID = cell_name, cell_type = assignment)
order_markers <- c("PanCK", "EGFR", "TP63", "GLUT1", "Cytochromec", "pH2Ax", "p16", "SLPI", "Ki67", "CA9", "LAMC2", "TNC", "Podoplanin", 
                   "Postn", "aSMA", "SPARC", "FN1", "CD31",
                   "CD15", "CD11b", "S100A9", "MMP9", "CD68", "HLADR", "CD163", 
                   "GranzymeB", "CD8")
order_celltype <- c("Epithelial", "Senescence", "Cell-Cycle", "Hypoxia", "pEMT", "Fibroblast", "Endothelial", "Neutrophil", "Monocyte", "M1-like_Macrophage", "M2-like_Macrophage", "Tctx-cell", "Undefined")
p <- plot_marker_by_celltype(
  rmat,
  cell_assignments = cell_assign3,
  expression_threshold = 0.3,
  marker_order = order_markers[order_markers %in% rownames(rmat)],
  celltype_order = order_celltype,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  mean_score_max = 0.5)
ggsave(filename = here("results/Diagnostic_Plots/CODEX_TMA/Cell_State_Marker_Scores_Dotplot_v3.pdf"), plot = p, device = "pdf", dpi = 300, width = 15, height = 8)

marker_sets <- list(
  `Tctx-cell` = list(pos = c("CD8","GranzymeB")),
  `M2-like_Macrophage` = list(pos = c("CD163","CD68","CD11b"),
                              neg = "CD15"),
  `M1-like_Macrophage` = list(pos = c("CD68","HLADR")),
  Monocyte = list(pos = c("CD11b","HLADR"),
                  neg = "CD68"),
  Neutrophil = list(pos = c("CD11b","CD15","S100A9")),
  Endothelial = list(pos = c("CD31")),
  Fibroblast = list(pos = c("aSMA","Postn","SPARC","FN1")),
  
  # states (optionally “negate other states” too, depends how strict you want)
  pEMT = list(
    pos = c("LAMC2","PanCK","EGFR","TP63"),
    neg = c("p16","SLPI")
  ),
  Hypoxia = list(
    pos = c("CA9","PanCK","EGFR","TP63"),
    neg = c("LAMC2","Ki67","p16","SLPI","S100A9")
  ),
  `Cell-Cycle` = list(
    pos = c("Ki67","PanCK","EGFR","TP63"),
    neg = c("LAMC2","CA9","p16","SLPI","S100A9")
  ),
  Senescence = list(
    pos = c("p16","SLPI","S100A9","PanCK","EGFR","TP63"),
    neg = c("LAMC2","CA9","Ki67")
  ),
  
  # epithelial baseline
  Epithelial = list(
    pos = c("PanCK","EGFR","TP63"),
    neg = c("LAMC2","CA9","p16","SLPI")
  )
)

order_celltype <- c("Epithelial", "Senescence", "Cell-Cycle", "Hypoxia", "pEMT", 
                    "Fibroblast", "Endothelial", "Neutrophil", "Monocyte", "M1-like_Macrophage", "M2-like_Macrophage", "Tctx-cell")

p <- plot_marker_sets_by_celltype(
  rmat = rmat,
  cell_assignments = cell_assign3,
  marker_sets = marker_sets,
  expression_threshold = 0.3,
  expression_threshold_neg = 0.5,
  fraction_rule = "pos_and_neg",
  fill_mode = "softmax",
  softmax_temp = 0.05,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  celltype_order = order_celltype,
  set_order = names(marker_sets),
  mean_score_max = 0.7
)

ggsave(filename = here("results/Diagnostic_Plots/CODEX_TMA/Cell_State_Marker_Sets_Dotplot_v2.pdf"), plot = p, device = "pdf", dpi = 300, width = 8.5, height = 7)





# Plot Cell-type/state Abundance Per Core ---------------------------------
cells <- read_csv(file = here("results/Generated_Data/CODEX/cell_type_and_state_assign_v2.csv"))
my_pal <- c("Fibroblast" = "#88CCEE", "Endothelial" = "#DDCC77", "B-cell" = "#AA4499", "Tctx-cell" = "sandybrown", 
            "Neutrophil" = "#DD1C77", "M2-like_Macrophage" = "#FF2400", "M1-like_Macrophage" = "pink2", "Monocyte" = "#1b9e77",
            "Myeloid" = "#ADDD8E", "Epithelial" = "#999933", "pEMT" = "#E4CAB5", "Senescence" = "#CC6677", "Cell-Cycle" = "#117733", "Hypoxia" = "#888888", "Undefined" = "grey85")

codex_parse_core_id <- function(x) sub("^.*?_([A-Z](?:10|[1-9])).*$", "\\1", x)

plot_df <- cells %>%
  mutate(core = codex_parse_core_id(cell_name)) %>%      # add core column
  count(core, assignment, name = "n_cells") %>%          # counts
  group_by(core) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%              # to fraction
  ungroup()

# optional: order cores
# core_levels <- c(pEMT_core, pEMT_edge)
# plot_df$core <- factor(plot_df$core, levels = core_levels)

plot_df$assignment <- factor(plot_df$assignment, levels = rev(c("pEMT", "Hypoxia", "Cell-Cycle", "Senescence", "Epithelial", "Myeloid", "Neutrophil", "M2-like_Macrophage", "M1-like_Macrophage", 
                                                                "Monocyte", "Tctx-cell", "B-cell", "Endothelial", "Fibroblast", "Undefined")))

ggplot(plot_df, aes(x = core, y = prop, fill = assignment)) +
  geom_col(width = 0.85, colour = "black", linewidth = 0.15) +
  scale_fill_manual(values = my_pal, name = "Cell type / state") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Composition per core (%)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position    = "right",
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )

per_sample_states <- as.data.frame(table(metadata$MPid[grep("pEMT|Hypoxia|Inflammatory", metadata$MPid), drop = TRUE], metadata$Sample[grep("pEMT|Hypoxia|Inflammatory", metadata$MPid)])) %>% 
  dplyr::mutate(Site = factor(samples_metadata$Ext_Site[match(Var2, samples_metadata$Sample)], levels = c("Laryngeal", "HPV+ OP", "HPV- OP", "Oral Cavity")), .before = "Freq") %>% 
  group_by(Var2) %>% mutate(Prop = Freq / sum(Freq), OrderBy = Prop[Var1 == "pEMT"])
p <- ggplot(data = per_sample_states, aes(x = reorder(Var2, OrderBy), y = Prop, fill = Var1)) + geom_col() + 
  scale_fill_manual(name = "States", values = state_cols[levels(per_sample_states$Var1)]) + scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0), limits = c(0, NA)) + 
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9, size = 16), plot.title = element_text(hjust = 0.5, size = 16),
                          legend.text = element_text(size = 14), legend.title = element_text(size = 16),
                          strip.background = element_blank(), strip.text = element_text(size = 14)) +
  facet_grid(.~Site, scales = 'free', space = 'free', drop = F) + guides(fill = guide_legend(override.aes = list(size = 7), reverse = TRUE))


per_core_abund <- cells %>%
  mutate(core = codex_parse_core_id(cell_name),
         assignment = case_when(
           assignment %in% c("Myeloid", "Neutrophil", "M2-like_Macrophage") ~ "Myeloid",
           TRUE ~ assignment
         )) %>%
  count(core, assignment, name = "n_cells") %>%
  group_by(core) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup()
by_assign <- per_core_abund %>%
  group_by(assignment) %>%
  summarise(
    Q1   = quantile(prop, 0.25, na.rm = TRUE),
    Q3   = quantile(prop, 0.75, na.rm = TRUE),
    IQR  = IQR(prop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(lower = pmax(0, Q1 - 1.5*IQR),   # clamp at 0
         upper = pmin(1, Q3 + 1.5*IQR))   # clamp at 1

f_plot <- per_core_abund %>%
  left_join(by_assign, by = "assignment") %>%
  mutate(is_outlier = prop < lower | prop > upper)

# Plot distribution with outlier cores labeled
ggplot(f_plot, aes(x = assignment, y = prop)) +
  geom_violin(aes(fill = assignment), width = 0.9, alpha = 0.25, colour = NA) +
  geom_boxplot(width = 0.4, outlier.shape = NA, colour = "grey30") +
  geom_jitter(aes(colour = assignment), width = 0.08, size = 2, alpha = 0.9) +
  ggrepel::geom_text_repel(
    data = dplyr::filter(f_plot, is_outlier),
    aes(label = core),
    size = 3, min.segment.length = 0.1, box.padding = 0.15, seed = 7
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Cell-type abundance per core") +
  theme_classic() +
  theme(legend.position = "none")




# Spatial Analysis --------------------------------------------------------
source(here("analysis_scripts/functions/CODEX_tools.R"))
my_pal <- c("Fibroblast" = "#88CCEE", "Endothelial" = "#DDCC77", "B-cell" = "#AA4499", "Tctx-cell" = "sandybrown", 
            "Neutrophil" = "#DD1C77", "M2-like_Macrophage" = "#FF2400", "M1-like_Macrophage" = "pink2", "Monocyte" = "#1b9e77",
            "Myeloid" = "#ADDD8E", "Epithelial" = "#999933", "pEMT" = "#E4CAB5", "Senescence" = "#CC6677", "Cell-Cycle" = "#117733", "Hypoxia" = "#888888", "Undefined" = "grey85")
codex_defaults$myeloid_levels <- c("Myeloid", "Neutrophil", "M2-like_Macrophage")
px2um <- 0.5067
flip_map <- c(A2 = TRUE, A3 = TRUE, A4 = FALSE, A5 = TRUE, 
              B2 = FALSE, B5 = TRUE, B6 = TRUE, B8 = TRUE, B10 = TRUE, 
              C2 = TRUE, C3 = TRUE, C4 = TRUE, C6 = TRUE, C7 = TRUE, C9 = TRUE, C10 = TRUE, 
              D2 = TRUE, D3 = TRUE, D4 = TRUE, D5 = FALSE, D6 = FALSE, D7 = FALSE, D8 = FALSE, D9 = FALSE, D10 = FALSE, 
              E3 = TRUE, E5 = TRUE)

### ============ Process each core and generate per-core results ============
## Run external scripts that run per-core processing and computation of myeloid infiltration and distance of pEMT cells from Fibroblasts and Myeloid 
## and saves the results to here("results/Generated_Data/CODEX/per_core_results")
# system2("bash", args = c(here("analysis_scripts/submit_full_pipeline_per_core.sh")))


### ============ Merge per-core results and analyze together ============
# cells <- read_csv(file = here("results/Generated_Data/CODEX/cell_type_and_state_assign_v2.csv"))
# rds_files <- list.files(here("results/Generated_Data/CODEX/per_core_results"), pattern = "^Core_.*\\.rds$", full.names = TRUE)
cells <- read_csv(file = here("results/Generated_Data/CODEX/cell_type_and_state_assign_v2.csv"))
rds_files <- list.files(here("results/Generated_Data/CODEX/per_core_results"), pattern = "^Core_.*\\.rds$", full.names = TRUE)

per_samp_pemt_abund <- cells %>% 
  dplyr::mutate(core_id = codex_parse_core_id(cell_name)) %>% 
  dplyr::group_by(core_id) %>% 
  dplyr::summarise(total = n(),
                   pEMT = sum(assignment == "pEMT"),
                   pemt_frac = pEMT / total,
                   .groups = "drop") %>% 
  dplyr::arrange(desc(pemt_frac))
low_qulity_cores <- c("A2", "D4", "E5", "A5", "C5")
low_pEMT_cores <- c("A3", "B10", "B2", "D8")

merge_res <- map_dfr(rds_files, function(f) {
  x <- readRDS(f)
  
  core <- x$core_id %||% sub("^Core_(.+)\\.rds$", "\\1", basename(f))
  
  dens <- tryCatch(
    x$Myeloid_infilt$summary$density[[1]],
    error = function(e) NA_real_
  )
  
  # Distance medians (rename to desired column names)
  dist_wide <- tryCatch({
    x$Distance_table %>%
      dplyr::filter(target %in% c("Myeloid","Fibroblast","Boundary")) %>%
      dplyr::select(target, median_um) %>%
      tidyr::pivot_wider(names_from = target, values_from = median_um) %>%
      dplyr::rename(
        med_pEMT_to_myeloid_um    = Myeloid,
        med_pEMT_to_fibroblast_um = Fibroblast,
        med_pEMT_to_boundary_um   = Boundary
      )
  }, error = function(e) tibble(
    med_pEMT_to_myeloid_um    = NA_real_,
    med_pEMT_to_fibroblast_um = NA_real_,
    med_pEMT_to_boundary_um   = NA_real_
  ))
  
  tibble(core_id = core, myeloid_density = dens) %>% bind_cols(dist_wide)
})

merge_res <- merge_res %>% 
  dplyr::mutate(group = dplyr::case_when(core_id %in% c("C3","B5","A4", "B6", "C2", "D5") ~ "pEMT_edge",
                                         core_id %in% c("C6","D7","E3") ~ "pEMT_core",
                                         core_id %in% c(low_qulity_cores, low_pEMT_cores) ~ "excluded",
                                         TRUE ~ "other")) %>% 
  dplyr::mutate(group = factor(group, levels = c("pEMT_core", "pEMT_edge", "other", "excluded")),
                myeloid_density = myeloid_density * 1000000) %>%                                # Convert density units from µm² to mm²
  dplyr::filter(group != "excluded") %>% 
  dplyr::arrange(desc(med_pEMT_to_boundary_um))



### ============ Plot the results - create figure panels ============
plot_df <- merge_res %>% dplyr::filter(group %in% c("pEMT_core", "pEMT_edge"))

# Panel 1 - Tumor Infiltrating Myeloid Cells
p1 <- codex_make_boxdot(plot_df, yvar = "myeloid_density", x_var = "group", ylab = "Infiltrating Myeloid Cells Per mm² Tumour Area", test_method = "wilcox.test", test_method_args = list(exact = TRUE),
                        binwidth = 10, expend_mult_y = 0.005, box_width = 0.6, dot_mm = 2, box_palette = c(Core = "#B24745FF", Edge = "#DF8F44FF"), dot_fill_fixed = "#0c2a50ff", dot_outline = "#E7E9EB")
ggsave(plot = p1, file = here("results/Paper_Figures/Revision_Figs/R1_C9_CODEX_pEMT_Edge_vs_Core_Myeloid_Infiltration_Boxplot.pdf"), device = "pdf", dpi = 300, width = 2, height = 3)

# Panel 2 - Average Distance From pEMT-cell to The Tumor Boundary
p2 <- codex_make_boxdot(plot_df, yvar = "med_pEMT_to_boundary_um", x_var = "group", ylab = "Mean pEMT to Tumor Boundary Distance (µm)", test_method = "t.test", test_method_args = list(var.equal = FALSE),
                        binwidth = 1.5, expend_mult_y = 0.01, box_width = 0.6, dot_mm = 2, box_palette = c(Core = "#B24745FF", Edge = "#DF8F44FF"), dot_fill_fixed = "#0c2a50ff", dot_outline = "#E7E9EB")
ggsave(plot = p2, file = here("results/Paper_Figures/Revision_Figs/R1_C9_CODEX_pEMT_Edge_vs_Core_pEMT_Distance_to_Tumor_Boundary.pdf"), device = "pdf", dpi = 300, width = 2, height = 3)

# Panel 3 - Distance From pEMT-cell to Either Fibroblast or Myeloid cells in pEMT-core vs pEMT-edge Tumors
# Filter out cores with extreme low abundance of Fibroblast or pro-tumor Myeloid (keep only cores where both fractions are >= 1st percentile of their respective distributions)
# Build wide table + joint abundance
core_abunds <- cells %>%
  mutate(core = codex_parse_core_id(cell_name),
         assignment = case_when(
           assignment %in% c("Myeloid", "Neutrophil", "M2-like_Macrophage") ~ "Myeloid",
           TRUE ~ assignment
         )) %>%
  count(core, assignment, name = "n_cells") %>%
  group_by(core) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup()

filt_cores <- core_abunds %>%
  dplyr::filter(assignment %in% c("Myeloid","Fibroblast")) %>%
  dplyr::select(core, assignment, prop) %>%
  tidyr::pivot_wider(names_from = assignment, values_from = prop) %>%
  dplyr::mutate(Myeloid = coalesce(Myeloid, 0),
                Fibroblast = coalesce(Fibroblast, 0)) %>% 
  dplyr::filter(Myeloid >= quantile(Myeloid, probs = 0.01) & Fibroblast >= quantile(Fibroblast, probs = 0.01)) %>% 
  dplyr::filter(core %in% plot_df$core_id)


# Plot panel 3
long_df <- plot_df %>%
  dplyr::select(core_id, group,
                Myeloid   = med_pEMT_to_myeloid_um,
                Fibroblast = med_pEMT_to_fibroblast_um) %>%
  tidyr::pivot_longer(c(Myeloid, Fibroblast),
                      names_to = "neighbor", values_to = "distance_um") %>%
  dplyr::mutate(neighbor = factor(neighbor, levels = c("Myeloid","Fibroblast"))) %>% 
  droplevels() %>% 
  filter(core_id %in% filt_cores$core)

# Plot: paired points + connecting lines - facet by group
strip_cols  <- c(pEMT_core = "#B24745FF", pEMT_edge = "#DF8F44FF")
df_in  <- filter(long_df, distance_um <= 80)
df_out <- filter(long_df, distance_um >  80)

p3 <- ggplot(long_df, aes(x = neighbor, y = distance_um)) +
  # boxes
  geom_boxplot(aes(fill = neighbor), width = 0.5,
               colour = "grey35", alpha = 0.85, outlier.shape = NA) +
  # paired connectors (will be clipped to the top)
  geom_line(aes(group = core_id), colour = "#0c2a50ff", linewidth = 0.6, alpha = 0.85) +
  # inlier points (circles)
  geom_point(data = df_in, shape = 21, size = 2, fill = "#0c2a50ff", colour = "#E7E9EB", stroke = 0.9) +
  # outlier points re-drawn as triangles pinned to the top
  geom_point(data = df_out, aes(y = 80),
             shape = 24, size = 2, fill = "#0c2a50ff", colour = "#E7E9EB", stroke = 0.9) +
  # outlier labels (print true value)
  geom_text(data = df_out, aes(y = 80, label = paste0(round(distance_um), "µm")),
            # position = position_nudge(x = -0.2), 
            vjust = -0.9, size = 3, colour = "#0c2a50ff") +
  # # stats at fixed height
  # ggpubr::stat_pvalue_manual(
  #   stats_df, xmin = "xmin", xmax = "xmax",
  #   y.position = "y.position", label = "label",
  #   tip.length = 0.01, size = 4, remove.bracket = TRUE 
  # ) +
  scale_fill_manual(values = c(Myeloid = "tan3", Fibroblast = "#88CCEE"), guide = "none") +
  scale_y_continuous(limits = c(20, 80), oob = scales::squish, expand = expansion(mult = c(0.02, 0.12))) +  
  labs(x = NULL, y = "Median distance from pEMT (µm)") +
  ggh4x::facet_wrap2(
    ~ group, nrow = 1,
    strip = ggh4x::strip_themed(
      background_x = lapply(strip_cols[levels(long_df$group)],
                            function(col) element_rect(fill = col, colour = "black", linewidth = 0.6)),
      text_x = lapply(strip_cols[levels(long_df$group)],
                      function(col) element_text(colour = "black"))
    )
  ) +
  (if (exists("common_theme")) common_theme else ggplot2::theme_classic()) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6), axis.line = element_blank()) +
  coord_cartesian(clip = "off") +
  ggh4x::force_panelsizes(
    rows = unit(60, "mm"),
    cols = rep(unit(30, "mm"), length(unique(long_df$group))),
    respect = TRUE
  )
ggsave(plot = p3, file = here("results/Paper_Figures/Revision_Figs/R1_C9_CODEX_pEMT_Edge_vs_Core_pEMT_Distance_to_Fibroblast_and_Myeloid_v2.pdf"), device = "pdf", dpi = 300, width = 4, height = 3)

# Compute statistics (difference to Myeloid and Fibfroblast in pEMT-core vs pEMT-edge samples)
# Myeloid: pEMT_core vs pEMT_edge
myeloid_test <- wilcox.test(
  distance_um ~ group,
  data = filter(long_df, neighbor == "Myeloid"),
  exact = FALSE, conf.int = TRUE
)

# Fibroblast: pEMT_core vs pEMT_edge
fibro_test <- wilcox.test(
  distance_um ~ group,
  data = filter(long_df, neighbor == "Fibroblast"),
  exact = FALSE, conf.int = TRUE
)

### ============ Plot all TMA cores with full annotation & boundary in a loop ============
processed_dir <- here("results/Generated_Data/CODEX/processed_cores")
out_dir <- here("results/Paper_Figures/Revision_Figs/CODEX_TMA/Cores_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
core_files <- list.files(processed_dir, pattern = "^Core_.*\\.rds$", full.names = TRUE)

my_pal <- c("Fibroblast" = "#88CCEE", "Endothelial" = "#DDCC77", "Tctx-cell" = "sandybrown", 
            "Neutrophil" = "#DD1C77", "M2-like_Macrophage" = "#FF2400", "M1-like_Macrophage" = "pink2", "Monocyte" = "#1b9e77",
            "Epithelial" = "#999933", "pEMT" = "#E4CAB5", "Senescence" = "#CC6677", "Cell-Cycle" = "#117733", "Hypoxia" = "#888888", "Undefined" = "grey85")

for (f in core_files) {
  core <- readRDS(f)
  core_id <- core$core_id
  core$prep$cells_sf <- core$prep$cells_sf %>% 
    dplyr::mutate(assignment = dplyr::case_when(assignment == "B-cell" ~ "Undefined",
                                                assignment == "Myeloid" ~ "M2-like_Macrophage",
                                                TRUE ~ assignment))
  p <- codex_plot_core(core$prep$cells_sf,
                       cell_palette = my_pal,
                       boundaries = core$nest_bound,
                       boundary_linewidths = 0.8,
                       boundary_colours = "white")
  ggplot2::ggsave(file.path(out_dir, sprintf("Core_%s_with_Boundary.png", core_id)), p, width = 14, height = 12, dpi = 300)
}


### ============ Group similar cell-types and export to QuPath for complete plotting ============
cells <- read_csv(file = here("results/Generated_Data/CODEX/cell_type_and_state_assign_v2.csv"))
cells <- cells %>%
  mutate(assignment = case_when(assignment %in% c("Neutrophil", "Myeloid", "M2-like_Macrophage") ~ "Pro-tumor Myeloid",
                                assignment %in% c("M1-like_Macrophage", "Monocyte") ~ "Anti-tumor Myeloid",
                                assignment == "Epithelial" ~ "Non-pEMT Malignant",
                                assignment == "pEMT" ~ "pEMT Malignant",
                                TRUE ~ assignment))
write_csv(cells, file = here("results/Generated_Data/CODEX/final_cell_type_assign.csv"), col_names = TRUE)




# Visium-like Spot Simulation -----------------------------------------------
processed_dir <- here("results/Generated_Data/CODEX/processed_cores")
out_dir <- here("results/Generated_Data/CODEX/visium_sim")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pix_per_um <- 1/0.5067

core_files <- list.files(processed_dir, pattern = "^Core_.*\\.rds$", full.names = TRUE)

for (f in core_files) {
  core <- readRDS(f)
  
  sim <- codex_simulate_visium_grid_core(
    prep_obj = core,
    pix_per_um = pix_per_um,
    spot_diameter_um = 55,
    spot_pitch_um = 100,
    roi_method = "convex_centroids",
    roi_buffer_um = 0,
    rotations_deg = seq(0, 30, 5),
    offset_steps = 6,
    require_full_spot_inside = TRUE
  )
  
  # cell → spot assignment + summaries
  c2s <- codex_assign_cells_to_spots(core$prep$cells_tbl, sim$spots_tbl)
  spot_sum <- codex_summarize_spots(core$prep$cells_tbl, sim$spots_tbl, c2s)
  
  # save RDS for downstream R analysis
  saveRDS(
    list(core_id = core$core_id, params = sim$params, roi = sim$roi,
         spots = spot_sum, cell_to_spot = c2s),
    file = file.path(out_dir, paste0("VisiumSim_", core$core_id, ".rds"))
  )
  
  # export QuPath CSV (centers in QuPath coordinate frame)
  codex_export_visium_grid_for_qupath(
    prep_obj = core,
    spots_tbl = spot_sum,
    out_csv = file.path(out_dir, paste0("VisiumGrid_", core$core_id, "_forQuPath.csv")),
    include_measurements = TRUE
  )
}

vis_sim_files <- list.files(out_dir, pattern = "^VisiumSim_.*\\.rds$", full.names = TRUE)
low_qulity_cores <- c("A2", "D4", "E5", "A5", "C5")
vis_sim_files <- vis_sim_files[!grepl(paste0("VisiumSim_(", paste(low_qulity_cores, collapse="|"), ")\\.rds$"), basename(vis_sim_files))]

# Read all cores and bind all spots
top_mps_tab_all <- lapply(vis_sim_files, extract_top_mps_from_visium_sim) %>% 
  purrr::compact() %>% 
  dplyr::bind_rows()

# Fraction of spots above a given proportion
grid_props <- seq(0, 1, 0.02)
percent_over_prop <- lapply(c("First_Score", "Second_Score", "Top_Two"), function(group) {
  data.frame(Cell_Prop = grid_props,
             Spot_Percent = sapply(grid_props, function(p) mean(top_mps_tab_all[[group]] >= p)),
             Group = group,
             stringsAsFactors = FALSE)}) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(Group = stringr::str_replace_all(Group, pattern = "Score", replacement = "MP"),
                Group = factor(Group, levels = c("Top_Two", "First_MP", "Second_MP")))

intersect_coords <- lapply(levels(percent_over_prop$Group), function(g) {
  group_data <- percent_over_prop[percent_over_prop$Group == g, ]
  ap <- approx(group_data$Cell_Prop, group_data$Spot_Percent, xout = 0.5)
  if (is.na(ap$y) || ap$y == 0) return(NULL)
  data.frame(x_int = c(-Inf, ap$x, ap$x),
             y_int = c(ap$y,  ap$y, -Inf),
             group_int = g,
             stringsAsFactors = FALSE)})
intersect_coords <- Filter(Negate(is.null), intersect_coords)

p <- ggplot(percent_over_prop, aes(x = Cell_Prop, y = Spot_Percent)) +
  stat_smooth(aes(x = Cell_Prop, y = Spot_Percent, colour = Group), method = "lm", formula = y ~ poly(x, 21), se = FALSE, linewidth = 2) + 
  geom_point(size = 3, color = "#0c2a50ff", stroke = 1, alpha = 0.6) + 
  sapply(seq_along(intersect_coords), function(i) geom_line(aes(x = x_int, y = y_int), linetype = 2, data = intersect_coords[[i]])) +
  xlab("Cell Type/State Proportion") + ylab("Spots Above a Given Proportion") + 
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(c(0.01, 0.01)), breaks = sort(c(0, 0.25, 0.50, 0.75, 1, sapply(seq_along(intersect_coords), function(i) intersect_coords[[i]]$y_int[[1]])))) +
  scale_x_continuous(expand = c(0, 0)) +
  scalop::theme_scalop(legend.position = "right") +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 2), axis.ticks.length.y = grid::unit(0.15, "cm"),
        axis.text.y = element_text(margin = margin(r = 0.1, unit = "cm")), panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(linewidth = 5)))
ggsave(filename = here("results/Paper_Figures/Revision_Figs/CODEX_TMA/Pseudo_Spot_Cell_Type_Dominancy_v2.pdf"), plot = p, device = "pdf", dpi = 300, width = 10, height = 8)
