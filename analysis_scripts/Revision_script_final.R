library(tidyverse)
library(patchwork)
library(scalop)
library(here)
source(here("analysis_scripts/functions/scRNA_Functions.R"))
source(here("analysis_scripts/functions/ST_Functions.R"))


# HPV Directed Analysis (R1C2) -------------------------------------------------

# Load per-sample parameter table
sample_parameters <- readRDS(file = here("metadata/per_sample_parameters_df.rds"))
# Load sample level metadata
samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))
# Load merged metadata tibble
metadata <- readRDS(file = here("metadata/merged_metadata_Extend.rds"))
# Load state colors 
state_cols <- read.delim(file = here("aux_data/state_col_tab_extend.tsv"), sep = "\t", header = TRUE)
state_cols <- setNames(state_cols$V2, state_cols$V1)


### ============ Intra-tumor - differences in total HPV reads between genetic subclones ============

include_classes <- c("Epithelial", "Mixed")
N_BINS  <- 30
N_REP   <- 100
SD_MULT <- 1
# Use Wilson-derived SE for observed p (TRUE), or simple binomial SE (FALSE)
USE_WILSON_SE <- FALSE

cpm_no_log <- function(mat) {
  lib <- colSums(mat); lib[lib == 0] <- 1
  sweep(mat, 2, lib, "/") * 1e6
}
safe_quantile_breaks <- function(x, nbins) {
  br <- stats::quantile(x, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE)
  unique(br)
}
make_bin_sampler <- function(bin2cands) {
  bins_avail <- as.integer(names(bin2cands))
  function(bin_id) {
    bid <- as.integer(bin_id)
    if (!length(bins_avail)) return(NA_character_)
    delta <- 0L
    repeat {
      cand_bins <- bins_avail[bins_avail %in% c(bid - delta, bid + delta)]
      cands <- unlist(bin2cands[as.character(cand_bins)], use.names = FALSE)
      if (length(cands)) return(sample(cands, size = 1, replace = TRUE))
      delta <- delta + 1L
      if (delta > 1000L) return(NA_character_)
    }
  }
}

# HPV+ samples
sample_parameters_hpvp <- sample_parameters %>%
  left_join(samples_metadata %>% dplyr::select(Sample, HPV), by = "Sample") %>%
  filter(HPV == "HPVpos")

# Per-sample loop (NA subclones removed)
obs_list <- list()
exp_list <- list()
hpv_summary_sample <- list()
exp_sample_list <- list()   # per-replicate expected overall fractions

for (i in seq_along(sample_parameters_hpvp$Sample)) {
  samp_id <- sample_parameters_hpvp$Sample[[i]]
  site_id <- sample_parameters_hpvp$Site[[i]]
  message("Processing sample: ", samp_id)
  
  # Load counts
  data_dir <- paste(here("data"), site_id, samp_id, sep = "/")
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

  # Keep malignant(+mixed) spots with NON-NA subclone
  samp_meta <- metadata %>%
    filter(Sample == samp_id, EpiStroma %in% include_classes, !is.na(Subclone)) %>%
    dplyr::select(SpotID, Subclone, Sample)
  
  spots_use <- intersect(colnames(m), samp_meta$SpotID)
  if (!length(spots_use)) next
  m <- m[, spots_use, drop = FALSE]
  
  # Per-sample expression bins
  cpm_s <- cpm_no_log(m)
  gene_mean_cpm_s <- rowMeans(cpm_s)
  names(gene_mean_cpm_s) <- rownames(m)
  
  pos_idx <- gene_mean_cpm_s > 0
  if (!any(pos_idx)) { warning("No positive-mean genes in ", samp_id); next }
  breaks  <- safe_quantile_breaks(gene_mean_cpm_s[pos_idx], N_BINS)
  if (length(breaks) < 2) breaks <- range(gene_mean_cpm_s[pos_idx])
  bins_pos <- cut(gene_mean_cpm_s[pos_idx], breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  gene_bins_s <- rep(1L, length(gene_mean_cpm_s)); names(gene_bins_s) <- names(gene_mean_cpm_s)
  gene_bins_s[pos_idx] <- bins_pos
  
  all_hpv_genes_s <- grep("^hpv", names(gene_bins_s), value = TRUE)
  non_hpv_s <- setdiff(names(gene_bins_s), all_hpv_genes_s)
  bin2cands_s <- split(non_hpv_s, gene_bins_s[non_hpv_s])
  sample_from_bin <- make_bin_sampler(bin2cands_s)
  
  # HPV genes present
  hpv_idx_all <- grep("^hpv", rownames(m))
  if (!length(hpv_idx_all)) next
  hpv_present <- rowSums(m[hpv_idx_all, , drop = FALSE]) > 0
  hpv_idx_present <- hpv_idx_all[hpv_present]
  hpv_genes_present <- rownames(m)[hpv_idx_present]
  if (!length(hpv_genes_present)) next
  
  # Observed per spot/subclone
  hpv_sum_per_spot <- colSums(m[hpv_idx_present, , drop = FALSE])
  hpv_pos <- hpv_sum_per_spot > 0
  
  hpv_summary_sample[[samp_id]] <- tibble(
    Sample = samp_id,
    n_spots = length(spots_use),
    n_HPVpos = sum(hpv_pos),
    obs_frac = mean(hpv_pos)
  )
  
  obs_sc <- samp_meta %>%
    left_join(tibble(SpotID = spots_use, HPVpos = as.logical(hpv_pos)), by = "SpotID") %>%
    group_by(Subclone) %>%
    summarise(n_spots = n(),
              n_pos   = sum(HPVpos, na.rm = TRUE),
              prop    = ifelse(n_spots > 0, n_pos / n_spots, NA_real_),
              .groups = "drop") %>%
    mutate(Sample = samp_id, Type = "Observed")
  
  # Expected (per-sample bin-matched) — N_REP replicates
  hpv_bins <- gene_bins_s[hpv_genes_present]
  gene_rows <- rownames(m)
  
  control_res <- vector("list", N_REP)
  rep_fracs <- numeric(N_REP)
  
  for (r in seq_len(N_REP)) {
    control_genes <- vapply(hpv_bins, sample_from_bin, FUN.VALUE = character(1))
    control_genes <- control_genes[!is.na(control_genes)]
    control_genes_in_mat <- intersect(control_genes, gene_rows)
    
    if (length(control_genes_in_mat) == 0) {
      ctrl_pos <- rep(FALSE, length(spots_use))
    } else {
      ctrl_sum <- colSums(m[control_genes_in_mat, , drop = FALSE])
      ctrl_pos <- ctrl_sum > 0
    }
    
    rep_fracs[r] <- mean(ctrl_pos)
    
    control_res[[r]] <- samp_meta %>%
      left_join(tibble(SpotID = spots_use, CTRLpos = as.logical(ctrl_pos)), by = "SpotID") %>%
      group_by(Subclone) %>%
      summarise(n_spots = n(),
                frac     = mean(CTRLpos, na.rm = TRUE),
                .groups  = "drop") %>%
      mutate(rep = r)
  }
  
  exp_sc <- bind_rows(control_res) %>%
    group_by(Subclone) %>%
    summarise(n_spots = max(n_spots),
              prop    = mean(frac, na.rm = TRUE),
              prop_sd = sd(frac,   na.rm = TRUE),
              .groups = "drop") %>%
    mutate(Sample = samp_id, Type = "Expected")
  
  obs_list[[samp_id]] <- obs_sc
  exp_list[[samp_id]] <- exp_sc
  exp_sample_list[[samp_id]] <- tibble(Sample = samp_id, rep = seq_len(N_REP), exp_frac = rep_fracs)
}

# collapse
obs_df <- bind_rows(obs_list)
exp_df <- bind_rows(exp_list)
hpv_summary_sample <- bind_rows(hpv_summary_sample) %>% arrange(desc(obs_frac))
exp_sample_df <- bind_rows(exp_sample_list)

# Build data for plotting
# 1) Subclone ordering within each sample (no 'Total' anywhere)
obs_order <- obs_df %>%
  dplyr::mutate(Subclone = as.character(Subclone)) %>%
  dplyr::arrange(Sample, prop) %>%
  dplyr::mutate(Subclone_key = paste(Sample, Subclone, sep = "|"))

levs <- obs_order$Subclone_key %>% unique()  # y-axis order for all panels

# 2) Observed: add Wilson CI per subclone
obs_plot_intra <- obs_df %>%
  dplyr::mutate(
    Subclone   = as.character(Subclone),
    Subclone_key = paste(Sample, Subclone, sep = "|"),
    Subclone_f = factor(Subclone_key, levels = levs)
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ci = list(binom::binom.confint(x = n_pos, n = n_spots, methods = "wilson"))) %>%
  tidyr::unnest_wider(ci) %>%
  dplyr::rename(prop_low = lower, prop_high = upper) %>%
  dplyr::ungroup()

# 3) Expected: points only (no ±SD bars)
exp_plot_intra <- exp_df %>%
  dplyr::mutate(
    Subclone   = as.character(Subclone),
    Subclone_key = paste(Sample, Subclone, sep = "|"),
    Subclone_f = factor(Subclone_key, levels = levs)
  )

# 4) Sample order (descending by overall observed fraction)
sample_order <- hpv_summary_sample %>%
  dplyr::arrange(dplyr::desc(obs_frac)) %>%
  dplyr::pull(Sample) %>% unique()

obs_plot_intra <- obs_plot_intra %>% dplyr::mutate(Sample_f = factor(Sample, levels = sample_order))
exp_plot_intra <- exp_plot_intra %>% dplyr::mutate(Sample_f = factor(Sample, levels = sample_order))

# 5) Δ and 95% CI for subclones only
# (uses N_REP from your earlier workflow; if not defined, set N_REP <- 10)
if (!exists("N_REP")) N_REP <- 10

delta_intra <- obs_df %>%
  dplyr::select(Sample, Subclone, n_spots, n_pos) %>%
  dplyr::mutate(p_obs = ifelse(n_spots > 0, n_pos / n_spots, NA_real_)) %>%
  dplyr::left_join(
    exp_df %>% dplyr::select(Sample, Subclone, mu_exp = prop, sd_exp = prop_sd),
    by = c("Sample","Subclone")
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    se_obs   = sqrt(pmax(p_obs * (1 - p_obs) / pmax(n_spots, 1L), 0)),
    se_exp   = sd_exp / sqrt(N_REP),
    se_delta = sqrt(se_obs^2 + se_exp^2),
    delta    = p_obs - mu_exp,
    lo_delta = delta - qnorm(0.975) * se_delta,
    hi_delta = delta + qnorm(0.975) * se_delta
  ) %>% dplyr::ungroup() %>%
  dplyr::mutate(
    Subclone   = as.character(Subclone),
    Subclone_key = paste(Sample, Subclone, sep = "|"),
    Subclone_f = factor(Subclone_key, levels = levs),
    Sample_f   = factor(Sample, levels = sample_order)
  )

# 6) Plots (left: observed vs expected; right: Δ with CI)
pos_obs <- ggplot2::position_nudge(y = -0.2)
pos_exp <- ggplot2::position_nudge(y =  +0.2)

p_left_intra <- ggplot2::ggplot() +
  # Expected (points only)
  ggplot2::geom_point(
    data = dplyr::filter(exp_plot_intra, !is.na(prop)),
    ggplot2::aes(x = prop, y = Subclone_f, shape = "Expected", color = "Expected"),
    size = 2.2, position = pos_exp
  ) +
  # Observed (Wilson CI + point)
  ggplot2::geom_errorbarh(
    data = dplyr::filter(obs_plot_intra, !is.na(prop)),
    ggplot2::aes(y = Subclone_f, xmin = prop_low, xmax = prop_high),
    height = 0.18, color = "black", position = pos_obs
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(obs_plot_intra, !is.na(prop)),
    ggplot2::aes(x = prop, y = Subclone_f, shape = "Observed", color = "Observed"),
    size = 2.2, position = pos_obs
  ) +
  ggplot2::facet_grid(rows = ggplot2::vars(Sample_f), scales = "free_y", space  = "free_y") +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_y_discrete(labels = function(k) sub(".*\\|", "", k)) +
  ggplot2::scale_shape_manual(values = c("Observed" = 17, "Expected" = 16), name = "Type") +
  ggplot2::scale_color_manual(values = c("Observed" = "#1b9e77", "Expected" = "#d95f02"), name = "Type") +
  ggplot2::labs(
    title = "HPV-positivity by subclone",
    subtitle = "observed vs per-sample bin-matched expected",
    x = "HPV+ spot proportion", y = ""
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "bottom", panel.spacing.y = grid::unit(0, "pt"))

# Δ axis limits/ticks
# break_i <- if (lim_i <= 0.05) 0.01 else 0.04
lim <- max(0.03, ceiling(max(abs(c(delta_intra$lo_delta, delta_intra$hi_delta)), na.rm = TRUE) * 100) / 100)
break_by   <- if (lim <= 0.05) 0.01 else 0.04
breaks_pct <- seq(-lim, lim, by = break_by)

p_right_intra <- ggplot2::ggplot(delta_intra) +
  ggplot2::geom_vline(xintercept = c(-0.02, 0.02), linetype = 2, color = "grey55") +
  ggplot2::geom_vline(xintercept = 0, linetype = 1, color = "#0c2a50ff") +
  ggplot2::geom_errorbarh(
    ggplot2::aes(y = Subclone_f, xmin = lo_delta, xmax = hi_delta, color = "Subclone"),
    height = 0.18
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = delta, y = Subclone_f, color = "Subclone", shape = "Subclone"),
    size = 2.2
  ) +
  ggplot2::facet_grid(rows = ggplot2::vars(Sample_f), scales = "free_y", space = "free_y") +
  ggplot2::scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = breaks_pct, limits = c(-lim, lim)
  ) +
  ggplot2::scale_y_discrete(labels = function(k) sub(".*\\|", "", k)) +
  ggplot2::scale_color_manual(values = c(Subclone = "#0072B2"), name = expression(Corrected~(Delta))) +
  ggplot2::scale_shape_manual(values = c(Subclone = 15), name = expression(Corrected~(Delta))) +
  ggplot2::labs(
    title = "Corrected difference",
    subtitle = expression(Delta == "Observed" - "Expected"),
    x = expression(Delta~"(percent)"), y = ""
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "bottom", panel.spacing.y = grid::unit(0, "pt"))

# Combined intra-tumor figure
p_intra <- (p_left_intra | p_right_intra) + patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(filename = here("results/Paper_Figures/Revision_Figs/PANEL_Per-sample_HPV_Positivity_Rate_Across_Subclones.png"), plot = p_intra, width = 12, height = 10, dpi = 300)
ggplot2::ggsave(filename = here("results/Paper_Figures/Revision_Figs/PANEL_Per-sample_HPV_Positivity_Rate_Across_Subclones.pdf"), plot = p_intra, width = 12, height = 10, dpi = 300, device = "pdf")



### ============ Inter-tumor - differences in total HPV reads between the samples ============

hpvp <- sample_parameters %>%
  dplyr::left_join(samples_metadata %>% dplyr::select(Sample, HPV), by = "Sample") %>%
  dplyr::filter(HPV == "HPVpos")

gene_sum_cpm <- NULL
total_spots  <- 0L

for (i in seq_along(hpvp$Sample)) {
  samp_id <- hpvp$Sample[[i]]
  site_id <- hpvp$Site[[i]]
  data_dir <- paste(here("data"), site_id, samp_id, sep = "/")
  m_i <- as.matrix(Matrix::readMM(paste(data_dir, grep("matrix(_[A-Za-z0-9]+)?\\.mtx$", list.files(data_dir), value = TRUE), sep = "/")))
  colnames(m_i) <- read.table(paste(data_dir, grep("barcodes(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 1]
  rownames(m_i) <- read.table(paste(data_dir, grep("features(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 2]
  if (sum(grep("GRCh38", rownames(m_i))) == 0) {
    m_i <- m_i
  } else {
    # Remove HPV genes and GRCh38 prefix
    m_i <- m_i[-grep("^hpv", rownames(m_i)), ]
    rownames(m_i) <- sub("GRCh38_", "", rownames(m_i))
  }

  spots_use_i <- metadata %>%
    filter(Sample == samp_id, EpiStroma %in% include_classes) %>%
    pull(SpotID) %>%
    intersect(colnames(m_i))
  if (!length(spots_use_i)) next
  
  cpm_i <- cpm_no_log(m_i[, spots_use_i, drop = FALSE])
  
  if (is.null(gene_sum_cpm)) {
    gene_sum_cpm <- setNames(rowSums(cpm_i), rownames(cpm_i))
  } else {
    common <- intersect(names(gene_sum_cpm), rownames(cpm_i))
    gene_sum_cpm[common] <- gene_sum_cpm[common] + rowSums(cpm_i[common, , drop = FALSE])
    new_genes <- setdiff(rownames(cpm_i), names(gene_sum_cpm))
    if (length(new_genes))
      gene_sum_cpm[new_genes] <- rowSums(cpm_i[new_genes, , drop = FALSE])
  }
  total_spots <- total_spots + ncol(cpm_i)
}

gene_mean_cpm <- gene_sum_cpm / max(total_spots, 1L)
names(gene_mean_cpm) <- names(gene_sum_cpm)

pos_idx <- gene_mean_cpm > 0
breaks  <- safe_quantile_breaks(gene_mean_cpm[pos_idx], N_BINS)
if (length(breaks) < 2) breaks <- range(gene_mean_cpm[pos_idx])
bins_pos <- cut(gene_mean_cpm[pos_idx], breaks = breaks, include.lowest = TRUE, labels = FALSE)

gene_bins <- rep(1L, length(gene_mean_cpm))           # zeros → lowest bin
names(gene_bins) <- names(gene_mean_cpm)
gene_bins[pos_idx] <- bins_pos

all_hpv_genes <- grep("^hpv", names(gene_bins), value = TRUE)
non_hpv_genes <- setdiff(names(gene_bins), all_hpv_genes)
bin2cands     <- split(non_hpv_genes, gene_bins[non_hpv_genes])
sample_from_bin <- make_bin_sampler(bin2cands)

# 2) Per-sample observed & expected (TOTAL only, global control)
get_expected_total_global <- function(samp_id, site_id,
                                      include_classes, gene_bins, sample_from_bin, N_REP = 100) {
  data_dir <- paste(here("data"), site_id, samp_id, sep = "/")
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

  spots_use <- metadata %>%
    filter(Sample == samp_id, EpiStroma %in% include_classes, !is.na(Subclone)) %>%
    pull(SpotID) %>%
    intersect(colnames(m))
  if (!length(spots_use)) return(NULL)
  
  m <- m[, spots_use, drop = FALSE]
  
  hpv_rows <- grep("^hpv", rownames(m))
  if (!length(hpv_rows)) return(NULL)
  
  # keep HPV genes with any counts in this sample
  hpv_present <- rowSums(m[hpv_rows, , drop = FALSE]) > 0
  hpv_genes_present <- rownames(m)[hpv_rows][hpv_present]
  if (!length(hpv_genes_present)) return(NULL)
  
  # observed TOTAL
  hpv_pos <- colSums(m[hpv_genes_present, , drop = FALSE]) > 0
  obs <- tibble(
    Sample   = samp_id,
    n_spots  = length(spots_use),
    n_HPVpos = sum(hpv_pos),
    obs_frac = mean(hpv_pos)
  )
  
  # expected TOTAL (global bins; mean across replicates)
  bins_this <- gene_bins[intersect(hpv_genes_present, names(gene_bins))]
  bins_this <- bins_this[!is.na(bins_this)]
  if (!length(bins_this)) return(NULL)
  
  rep_fracs <- replicate(N_REP, {
    ctrl_genes <- vapply(bins_this, sample_from_bin, FUN.VALUE = character(1))
    ctrl_genes <- ctrl_genes[!is.na(ctrl_genes)]
    ctrl_in    <- intersect(ctrl_genes, rownames(m))
    if (!length(ctrl_in)) return(0)
    mean(colSums(m[ctrl_in, , drop = FALSE]) > 0)
  })
  
  exp <- tibble(
    Sample   = samp_id,
    exp_frac = mean(rep_fracs, na.rm = TRUE)   # <- average only; SD not used for plotting
  )
  
  list(obs = obs, exp = exp)
}

res_list <- purrr::map2(
  hpvp$Sample, hpvp$Site,
  ~ get_expected_total_global(.x, .y, include_classes, gene_bins, sample_from_bin, N_REP)
)

obs_total_global <- bind_rows(purrr::map(res_list, "obs"))
exp_total_global <- bind_rows(purrr::map(res_list, "exp"))

# order samples by observed fraction (desc)
sample_order <- obs_total_global %>% arrange(desc(obs_frac)) %>% pull(Sample) %>% rev

# Observed Wilson CI
obs_plot_tot <- obs_total_global %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ci = list(binom::binom.confint(n_HPVpos, n_spots, methods = "wilson"))) %>%
  tidyr::unnest_wider(ci) %>%
  dplyr::rename(prop = obs_frac, prop_low = lower, prop_high = upper) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Sample_f = factor(Sample, levels = sample_order))

# Expected (mean only)
exp_plot_tot <- exp_total_global %>%
  dplyr::mutate(Sample_f = factor(Sample, levels = sample_order))

# Δ and its CI (from observed Wilson only)
delta_tot <- obs_plot_tot %>%
  dplyr::select(Sample, Sample_f, prop_low, prop_high, prop) %>%
  dplyr::left_join(exp_total_global, by = "Sample") %>%
  dplyr::transmute(
    Sample, Sample_f,
    delta    = prop - exp_frac,
    lo_delta = prop_low - exp_frac,
    hi_delta = prop_high - exp_frac
  )

# 3) Plots
# LEFT: observed (with Wilson CI) vs expected (points only)
p_left <- ggplot() +
  geom_errorbarh(
    data = obs_plot_tot,
    aes(y = Sample_f, xmin = prop_low, xmax = prop_high),
    height = 0.18, color = "black"
  ) +
  geom_point(
    data = obs_plot_tot,
    aes(x = prop, y = Sample_f, shape = "Observed", color = "Observed"),
    size = 2.2
  ) +
  geom_point(
    data = exp_plot_tot,
    aes(x = exp_frac, y = Sample_f, shape = "Expected", color = "Expected"),
    size = 2.2
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_shape_manual(values = c("Observed" = 17, "Expected" = 16), name = "Type") +
  scale_color_manual(values = c("Observed" = "#1b9e77", "Expected" = "#d95f02"), name = "Type") +
  labs(
    title = "HPV-positivity by sample",
    subtitle = "observed vs global bin-matched expected",
    x = "HPV+ spot proportion", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# RIGHT: Δ = observed − expected (CI from observed Wilson only)
lim <- max(0.03, ceiling(max(abs(c(delta_tot$lo_delta, delta_tot$hi_delta)), na.rm = TRUE) * 100) / 100)

p_right <- ggplot(delta_tot) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = 2, color = "grey55") +
  geom_vline(xintercept = 0, linetype = 1, color = "#0c2a50ff") +
  geom_errorbarh(
    aes(y = Sample_f, xmin = lo_delta, xmax = hi_delta, color = "Δ"),
    height = 0.18
  ) +
  geom_point(
    aes(x = delta, y = Sample_f, color = "Δ", shape = "Δ"),
    size = 2.2
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = c(-0.5, -0.3, -0.1, 0.1, 0.3, 0.5), limits = c(-lim, lim)) +
  scale_shape_manual(values = c("Δ" = 15), name = expression(Corrected~(Delta))) +
  scale_color_manual(values = c("Δ" = "#0072B2"), name = expression(Corrected~(Delta))) +
  labs(
    title = "Corrected difference",
    subtitle = expression(Delta == "Observed" - "Expected"),
    x = expression(Delta~"(percent)"), y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Combined figure (side by side)
p_inter <- (p_left | p_right) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggplot2::ggsave(filename = here("results/Paper_Figures/Revision_Figs/PANEL_HPV_Positivity_Rate_Across_Samples.png"), plot = p_inter, width = 12, height = 10, dpi = 300)
ggplot2::ggsave(filename = here("results/Paper_Figures/Revision_Figs/PANEL_HPV_Positivity_Rate_Across_Samples.pdf"), plot = p_inter, width = 12, height = 10, dpi = 300, device = "pdf")

p_inter_noleg <-(p_left + guides(colour = "none", shape = "none", fill = "none", linetype = "none")) |
  (p_right + guides(colour = "none", shape = "none", fill = "none", linetype = "none"))
p_intra_leg <- (p_left_intra | p_right_intra) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_both <- (p_inter_noleg / patchwork::plot_spacer() / p_intra_leg) +
  patchwork::plot_layout(heights = c(0.35, 0.005, 1))
ggplot2::ggsave(filename = here("results/Paper_Figures/Revision_Figs/PANEL_HPV_Positivity_Rate_Both_Inter_and_Intra_Tumor.pdf"), plot = p_both, width = 10, height = 10, dpi = 300, device = "pdf")



### ============ Diagnostic: HPV genes in global bins (per sample) ============
## 1) Global mean CPM + bins (named)
global_sum <- setNames(numeric(0), character(0))
global_n   <- setNames(integer(0), character(0))

for (i in seq_along(sample_parameters_hpvp$Sample)) {
  samp_id <- sample_parameters_hpvp$Sample[[i]]
  site_id <- sample_parameters_hpvp$Site[[i]]
  data_dir <- paste(here("data"), site_id, samp_id, sep = "/")
  m_i <- as.matrix(Matrix::readMM(paste(data_dir, grep("matrix(_[A-Za-z0-9]+)?\\.mtx$", list.files(data_dir), value = TRUE), sep = "/")))
  colnames(m_i) <- read.table(paste(data_dir, grep("barcodes(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 1]
  rownames(m_i) <- read.table(paste(data_dir, grep("features(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 2]
  if (sum(grep("GRCh38", rownames(m_i))) == 0) {
    m_i <- m_i
  } else {
    # Remove HPV genes and GRCh38 prefix
    m_i <- m_i[-grep("^hpv", rownames(m_i)), ]
    rownames(m_i) <- sub("GRCh38_", "", rownames(m_i))
  }

  spots_use <- metadata |>
    dplyr::filter(Sample == samp_id, EpiStroma %in% include_classes, !is.na(Subclone)) |>
    dplyr::pull(SpotID) |>
    intersect(colnames(m_i))
  if (!length(spots_use)) next
  
  cpm_i <- cpm_no_log(m_i[, spots_use, drop = FALSE])
  gsum  <- rowSums(cpm_i); genes <- names(gsum)
  
  new_genes <- setdiff(genes, names(global_sum))
  if (length(new_genes)) { global_sum[new_genes] <- 0; global_n[new_genes] <- 0L }
  global_sum[genes] <- global_sum[genes] + gsum
  global_n[genes]   <- global_n[genes]   + ncol(cpm_i)
}

global_mean_cpm <- global_sum / pmax(global_n, 1L)
stopifnot(length(global_mean_cpm) > 0L, !is.null(names(global_mean_cpm)))

pos_idx_g <- global_mean_cpm > 0
br_g <- safe_quantile_breaks(global_mean_cpm[pos_idx_g], N_BINS)
if (length(br_g) < 2) br_g <- range(global_mean_cpm[pos_idx_g])

bins_global <- rep(1L, length(global_mean_cpm)); names(bins_global) <- names(global_mean_cpm)
bins_global[pos_idx_g] <- cut(global_mean_cpm[pos_idx_g], breaks = br_g, include.lowest = TRUE, labels = FALSE)

## 2) Diagnostics – non-HPV bin sizes
all_genes      <- names(global_mean_cpm)
hpv_flag       <- grepl("^hpv", all_genes, perl = TRUE)
non_hpv_global <- all_genes[!hpv_flag]

bin_sizes <- tibble::tibble(bin = as.integer(bins_global[non_hpv_global])) |>
  dplyr::count(bin, name = "Freq") |>
  tidyr::complete(bin = seq_len(N_BINS), fill = list(Freq = 0L)) |>
  dplyr::mutate(frac = Freq / sum(Freq))

## 3) Map each sample’s HPV genes to global bins
get_prod <- function(g) sub("^hpv[^_]+__", "", g)

hpv_bins_by_sample <- purrr::map_dfr(seq_along(sample_parameters_hpvp$Sample), function(i) {
  samp_id <- sample_parameters_hpvp$Sample[[i]]
  site_id <- sample_parameters_hpvp$Site[[i]]
  data_dir <- paste(here("data"), site_id, samp_id, sep = "/")
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

  spots_use <- metadata |>
    dplyr::filter(Sample == samp_id, EpiStroma %in% include_classes) |>
    dplyr::pull(SpotID) |>
    intersect(colnames(m))
  if (!length(spots_use)) return(NULL)
  m <- m[, spots_use, drop = FALSE]
  
  hpv_rows <- grep("^hpv", rownames(m))
  if (!length(hpv_rows)) return(NULL)
  hpv_present <- rowSums(m[hpv_rows, , drop = FALSE]) > 0
  hpv_genes   <- rownames(m)[hpv_rows][hpv_present]
  if (!length(hpv_genes)) return(NULL)
  
  tibble::tibble(
    Sample  = samp_id,
    gene    = hpv_genes,
    product = factor(get_prod(hpv_genes), levels = c("E1","E2","E5","E6","E7","L1","L2")),
    bin     = as.integer(bins_global[hpv_genes])
  ) |>
    dplyr::filter(!is.na(bin))
})

if (nrow(hpv_bins_by_sample) == 0) {
  stop("No HPV genes mapped to global bins. Check filters or input matrices.")
}

# sample order (reuse earlier if you have it)
samp_levels <- hpv_bins_by_sample$Sample |> unique() |> sort() |> rev()
hpv_bins_by_sample <- hpv_bins_by_sample |>
  dplyr::mutate(Sample_f = factor(Sample, levels = samp_levels))

## 4) Plots
# Heatmap: #HPV genes per bin per sample
hpv_counts_per_bin <- hpv_bins_by_sample |>
  dplyr::count(Sample_f, bin, name = "n")

gg_hmap <- ggplot2::ggplot(hpv_counts_per_bin,
                           ggplot2::aes(x = bin, y = Sample_f, fill = n)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.2) +
  ggplot2::scale_fill_viridis_c(option = "C", direction = 1,
                                limits = c(0, max(hpv_counts_per_bin$n, 7L)),
                                name = "# HPV genes") +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10),
                              limits = c(1, N_BINS)) +
  ggplot2::labs(title = "HPV genes per global expression bin",
                subtitle = "heatmap of counts per sample",
                x = "Global expression bin", y = NULL) +
  ggplot2::theme_minimal(12) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 legend.position = "right")

# Dot panel: which HPV product sits in which bin (per sample)
gg_dots <- ggplot2::ggplot(hpv_bins_by_sample,
                           ggplot2::aes(x = bin, y = product, color = product)) +
  ggplot2::geom_point(size = 2.4, position = ggplot2::position_jitter(width = 0.15, height = 0)) +
  ggplot2::facet_grid(rows = ggplot2::vars(Sample_f), scales = "free_y", space = "free_y") +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10),
                              limits = c(1, N_BINS)) +
  ggsci::scale_color_jama(name = "HPV gene") +
  ggplot2::labs(title = "HPV products mapped to global expression bins",
                subtitle = "one dot per HPV gene present in the sample",
                x = "Global expression bin", y = "HPV gene") +
  ggplot2::theme_minimal(12) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 legend.position = "right")

# Optional: combine
diag_global_bins <- gg_hmap / gg_dots + patchwork::plot_layout(heights = c(0.45, 0.55))


# Are control bins large enough?
bins_used <- hpv_bins_by_sample$bin |> unique()
bin_sizes_used <- bin_sizes |> dplyr::filter(bin %in% bins_used)
summary(bin_sizes_used$Freq)
min_controls <- min(bin_sizes_used$Freq)
min_controls

# Are those controls present in each sample’s matrix?
check_controls_per_sample <- purrr::map_dfr(sample_parameters_hpvp$Sample, function(samp_id){
  site_id <- sample_parameters_hpvp$Site[sample_parameters_hpvp$Sample == samp_id][1]
  data_dir <- paste(here("data"), site_id, samp_id, sep = "/")
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
  in_mat <- intersect(names(bins_global), rownames(m))
  tibble::tibble(Sample=samp_id,
                 bin = as.integer(bins_global[in_mat])) |>
    dplyr::count(Sample, bin, name="n_controls_in_mat")
})
# Check bins used by HPV genes
check_controls_per_sample |>
  dplyr::filter(bin %in% bins_used) |>
  dplyr::group_by(bin) |>
  dplyr::summarise(min_in_mat = min(n_controls_in_mat), median_in_mat = median(n_controls_in_mat))




### ============ Differences in HPV-gene expression patterns between genetic subclones ============
get_strain <- function(g) sub("^((hpv[^_]+))__.*$", "\\1", g)
get_prod   <- function(g) sub("^hpv[^_]+__", "", g)
make_xpos <- function(ids, step, start = 1) {
  tibble::tibble(id = ids,
                 x  = (seq_along(ids) - 1) * step + start)
}

new_w  <- 0.50
gap    <- 0.08
step   <- new_w + gap 

# which samples to run
hpv_samples <- samples_metadata %>% filter(HPV == "HPVpos") %>% pull(Sample)
samples_to_run <- hpv_samples
clone_groups_map <- list(P5710 = list(E7 = c("1","4","5"), E2E6 = c("2","3")),
                         P5735 = list(L1L2 = c("2","3"), E2 = c("1")),
                         P5766 = list(L1_hi = c("1"), L1_lo = c("2"), Balance = c("3")))

results <- purrr::map(samples_to_run, function(target_sample) {
  site_id <- sample_parameters %>% filter(Sample == target_sample) %>% pull(Site) %>% unique() %>% .[[1]]
  
  # dense counts
  data_dir <- paste(here("data"), site_id, target_sample, sep = "/")
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

  # Epithelial/Mixed + non-NA subclone
  s_meta <- metadata %>%
    filter(Sample == target_sample, EpiStroma %in% c("Epithelial","Mixed"), !is.na(Subclone)) %>%
    dplyr::select(SpotID, Subclone, Sample, MPid)
  
  spots_use <- intersect(colnames(m), s_meta$SpotID)
  if (!length(spots_use)) return(list(sample = target_sample, p_hpv = NULL, p_states = NULL))
  m <- m[, spots_use, drop = FALSE]
  s_meta <- s_meta %>% filter(SpotID %in% spots_use)
  
  # HPV genes present; dominant strain
  hpv_rows <- grep("^hpv", rownames(m))
  if (!length(hpv_rows)) return(list(sample = target_sample, p_hpv = NULL, p_states = NULL))
  hpv_totals_all <- rowSums(m[hpv_rows, , drop = FALSE])
  hpv_present <- hpv_totals_all > 0
  hpv_genes_present <- rownames(m)[hpv_rows][hpv_present]
  
  hpv_gene_table <- tibble(
    gene      = hpv_genes_present,
    strain    = get_strain(hpv_genes_present),
    product   = get_prod(hpv_genes_present),
    total     = as.numeric(hpv_totals_all[hpv_present])
  )
  
  strain_totals <- hpv_gene_table %>% group_by(strain) %>% summarise(total = sum(total), .groups="drop") %>% arrange(desc(total))
  dominant_strain <- strain_totals$strain[1]
  dom_genes <- hpv_gene_table %>% filter(strain == dominant_strain) %>% pull(gene)
  if (!length(dom_genes)) return(list(sample = target_sample, p_hpv = NULL, p_states = NULL))
  
  # pseudobulk counts per subclone (genes x subclones)
  sub_levels <- s_meta$Subclone %>% as.character() %>% unique() %>% sort()
  pb_counts <- matrix(0, nrow = length(dom_genes), ncol = length(sub_levels),
                      dimnames = list(dom_genes, sub_levels))
  for (sc in sub_levels) {
    sc_spots <- s_meta %>% filter(Subclone == sc) %>% pull(SpotID)
    if (length(sc_spots)) pb_counts[, sc] <- rowSums(m[dom_genes, sc_spots, drop = FALSE])
  }
  
  pb_long <- as.data.frame(pb_counts) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "Subclone", values_to = "counts") %>%
    mutate(product = get_prod(gene))
  
  pb_sub_totals <- pb_long %>% group_by(Subclone) %>% summarise(total_hpvc = sum(counts), .groups="drop") %>% arrange(desc(total_hpvc))
  
  pb_comp <- pb_long %>% left_join(pb_sub_totals, by = "Subclone") %>%
    mutate(prop = if_else(total_hpvc > 0, counts / total_hpvc, NA_real_))
  
  prod_levels <- c("E1","E2","E5","E6","E7","L1","L2")
  pb_comp_plot <- pb_comp %>% filter(!is.na(prop)) %>% mutate(product = factor(product, levels = prod_levels))
  
  # Two-panel HPV plot (counts on top, composition below)
  sub_order <- pb_sub_totals$Subclone %>% as.character() %>% as.numeric() %>% sort() %>% as.character()
  x_pos <- make_xpos(sub_order, step) %>% dplyr::rename(Subclone = id)
  xmin <- min(x_pos$x) - new_w/2
  xmax <- max(x_pos$x) + new_w/2
  top_pad_comp <- 0.01
  
  pb_counts_df <- pb_sub_totals %>% mutate(Subclone = as.character(Subclone)) %>% left_join(x_pos, by = "Subclone")
  pb_comp_df   <- pb_comp_plot %>% mutate(Subclone = as.character(Subclone)) %>% left_join(x_pos, by = "Subclone")
  
  hpv_sum_dom <- colSums(m[dom_genes, spots_use, drop = FALSE])
  hpv_pos_dom <- hpv_sum_dom > 0
  
  spot_stats <- s_meta %>%
    dplyr::select(SpotID, Subclone) %>%
    dplyr::mutate(hpv_total = as.numeric(hpv_sum_dom[SpotID]),
                  hpv_pos   = hpv_pos_dom[SpotID])
  
  sub_stats <- spot_stats %>%
    dplyr::group_by(Subclone) %>%
    dplyr::summarise(
      n_spots    = dplyr::n(),
      n_HPVpos   = sum(hpv_pos, na.rm = TRUE),
      mean_hpvc  = mean(hpv_total[hpv_pos], na.rm = TRUE),   # mean across ALL spots (zeros included)
      .groups = "drop"
    ) %>%
    dplyr::mutate(Subclone = as.character(Subclone)) %>%
    dplyr::left_join(x_pos, by = "Subclone")
  
  # y-axis headroom
  ymax_pos  <- max(sub_stats$n_HPVpos,  na.rm = TRUE) * 1.03
  ymax_mean <- max(sub_stats$mean_hpvc, na.rm = TRUE) * 1.03
  
  # Top panel: # HPV+ spots per subclone
  posneg_df <- sub_stats %>%
    dplyr::mutate(n_HPVneg = pmax(n_spots - n_HPVpos, 0L)) %>%
    dplyr::select(Subclone, x, HPV_pos = n_HPVpos, HPV_neg = n_HPVneg) %>%
    tidyr::pivot_longer(cols = c(HPV_pos, HPV_neg),
                        names_to = "HPV_status", values_to = "n")
  
  # nice upper limit for total-height axis
  ymax_total <- max(sub_stats$n_spots, na.rm = TRUE)
  
  p_pos <- ggplot(posneg_df, aes(x = x, y = n, fill = HPV_status)) +
    geom_col(width = new_w) +
    scale_fill_manual(
      name   = NULL,
      values = c(HPV_pos = "gray35", HPV_neg = "gray80"),
      labels = c(HPV_pos = "HPV-positive", HPV_neg = "HPV-negative")
    ) +
    scale_x_continuous(limits = c(xmin, xmax),
                       breaks = x_pos$x, labels = x_pos$Subclone,
                       expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, ymax_total), expand = c(0, 0)) +
    labs(x = NULL, y = "Spots") +
    theme_minimal(12) +
    theme(
      panel.grid      = element_blank(),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      plot.margin     = margin(5.5, 5.5, 0, 5.5),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.margin = margin(b = 2),
      legend.text = element_text(size = 11)
    )
  
  # Middle panel: mean HPV counts per spot (across all spots)
  # If you want mean among HPV+ spots only, change to mean(hpv_total[hpv_pos], na.rm = TRUE)
  p_mean <- ggplot(sub_stats, aes(x = x, y = mean_hpvc)) +
    geom_col(width = new_w, fill = "#0c2a50ff") +
    scale_x_continuous(limits = c(xmin, xmax),
                       breaks = x_pos$x, labels = x_pos$Subclone,
                       expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, ymax_mean), expand = c(0, 0)) +
    labs(x = NULL, y = "Mean HPV\ncounts / spot") +
    theme_minimal(12) +
    theme(
      panel.grid   = element_blank(),
      legend.position = "none",
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin  = margin(0, 5.5, 0, 5.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
  
  p_comp <- ggplot(pb_comp_df, aes(x = x, y = prop, fill = product)) +
    geom_col(width = new_w) +
    scale_x_continuous(limits = c(xmin, xmax), breaks = x_pos$x, labels = x_pos$Subclone, expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, top_pad_comp)),
                       breaks = seq(0, 1, 0.25), labels = scales::percent_format()) +
    ggsci::scale_fill_jama(name = "HPV gene") +
    labs(x = "Subclone", y = "Within-subclone\nHPV composition") +
    theme_minimal(12) +
    theme(panel.grid = element_blank(),
          plot.margin = margin(0, 5.5, 5.5, 5.5),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8))
  
  p_hpv <- (p_pos / p_mean / patchwork::plot_spacer() / p_comp) +
    patchwork::plot_layout(heights = c(0.24, 0.24, 0.005, 1), guides = "keep")
  
  p_hpv <- cowplot::ggdraw(p_hpv) +
    theme(plot.background = element_rect(color = "black", fill = NA, linewidth = 1.4))
  ggplot2::ggsave(filename = paste0(here("results/Paper_Figures/Revision_Figs/"), target_sample, "_Per-Subclone_Relative_HPV_Gene_Expression.png"), plot = p_hpv, width = 8, height = 10, dpi = 300)
  ggplot2::ggsave(filename = paste0(here("results/Paper_Figures/Revision_Figs/"), target_sample, "_Per-Subclone_Relative_HPV_Gene_Expression.pdf"), plot = p_hpv, width = 7, height = 10, dpi = 300, device = "pdf")
  
  # Group similar HPV-gene-expressing groups and check their cell-state abundance
  # use the predefined subclone grouping list "clone_groups_map" to automatically group subclones
  grp_map <- clone_groups_map[[target_sample]]
  
  if (!is.null(grp_map) && length(grp_map) > 0) {
    # named vector: names = subclone IDs, values = group names
    sub2grp <- setNames(
      rep(names(grp_map), lengths(grp_map)),
      unlist(grp_map, use.names = FALSE)
    )
    
    s_meta <- s_meta %>%
      dplyr::mutate(
        Subclone   = as.character(Subclone),
        HPV_clones = unname(sub2grp[Subclone]),                    # map
        HPV_clones = ifelse(is.na(HPV_clones), "Other", HPV_clones) # fill
      )
    
    # drop "Other" if < 50 spots BEFORE factoring
    other_n <- sum(s_meta$HPV_clones == "Other", na.rm = TRUE)
    if (other_n < 50) {
      s_meta <- dplyr::filter(s_meta, HPV_clones != "Other")
      group_levels <- names(grp_map)
    } else {
      group_levels <- c(names(grp_map), "Other")
    }
    
    s_meta <- s_meta %>%
      dplyr::mutate(HPV_clones = factor(HPV_clones, levels = group_levels))
    
    # Use only the mapped groups for tests/plots
    groups_keep <- names(grp_map)
    
  } else {
    # No mapping provided for this sample
    s_meta <- s_meta %>%
      dplyr::mutate(HPV_clones = factor("Other", levels = "Other"))
    groups_keep <- character(0)
  }
  
  per_group_states <- s_meta %>% dplyr::select(HPV_clones, MPid) %>% table() %>% as.data.frame() %>% 
    group_by(HPV_clones) %>% mutate(Prop = Freq / sum(Freq)) %>% ungroup
  
  mp_plot_df <- per_group_states %>%
    dplyr::mutate(
      Group = ifelse(MPid %in% c("T_cell","B_cell","Macrophage","Fibroblast",
                                 "Endothelial","Complement","Skeletal_Muscle","Secretory_Norm"),
                     "Non-Malignant", "Malignant")
    ) %>%
    dplyr::group_by(MPid, HPV_clones) %>%
    dplyr::summarise(
      Malignant      = sum(Prop[Group == "Malignant"]),
      `Non-Malignant`= -sum(Prop[Group == "Non-Malignant"]),
      .groups = "drop"
    )
  
  # 2) Fix the group order (same as your Ord logic), then compute *numeric* x positions
  ord_df <- mp_plot_df %>%
    dplyr::group_by(HPV_clones) %>%
    dplyr::summarise(Ord = sum(`Non-Malignant`), .groups="drop") %>%
    dplyr::arrange(Ord)
  
  grp_order <- ord_df$HPV_clones
  x_pos_states <- make_xpos(grp_order, step) %>% dplyr::rename(HPV_clones = id)
  
  xmin_s <- min(x_pos_states$x) - new_w/2
  xmax_s <- max(x_pos_states$x) + new_w/2
  
  mp_plot_df <- mp_plot_df %>%
    dplyr::left_join(ord_df, by = "HPV_clones") %>%
    dplyr::left_join(x_pos_states, by = "HPV_clones")
  
  # 3) Plot with the SAME bar width (new_w) and SAME spacing (via x numeric scale)
  p_states <- ggplot(mp_plot_df, aes(x = x)) +
    # positive (Malignant) stacked upward; reverse to match your style
    geom_col(aes(y = Malignant, fill = MPid),
             width = new_w, position = position_stack(reverse = TRUE)) +
    # negative (Non-Malignant) stacked downward
    geom_col(aes(y = `Non-Malignant`, fill = MPid),
             width = new_w, position = "stack") +
    geom_hline(yintercept = 0) +
    scale_x_continuous(limits = c(xmin_s, xmax_s),
                       breaks = x_pos_states$x,
                       labels = x_pos_states$HPV_clones,
                       expand = c(0.01, 0.01)) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    scale_fill_manual(name = "States",
                      values = state_cols[levels(per_group_states$MPid)]) +
    theme_classic() +
    theme(
      axis.title   = element_blank(),
      axis.text.x  = element_text(size = 16),
      plot.title   = element_text(hjust = 0.5, size = 16),
      legend.text  = element_text(size = 14),
      legend.title = element_text(size = 16)
    ) +
    guides(fill = guide_legend(override.aes = list(size = 7), reverse = TRUE))
  ggplot2::ggsave(filename = paste0(here("results/Paper_Figures/Revision_Figs/"), target_sample, "__Per-Subclone_Relative_HPV_Gene_Expression_MP_Enrichment.png"), plot = p_states, width = 8, height = 10, dpi = 300)
  ggplot2::ggsave(filename = paste0(here("results/Paper_Figures/Revision_Figs/"), target_sample, "__Per-Subclone_Relative_HPV_Gene_Expression_MP_Enrichment.pdf"), plot = p_states, width = 7, height = 10, dpi = 300, device = "pdf")
})



### ============ Per-spot HPV gene expression patterns association with metaprogram identity ============
HPV_DET_THRESHOLD <- 1     # >=1 UMI → HPV+
PSEUDOCOUNT       <- 0.5   # for CLR
K_RANGE           <- 2:6   # candidate clusters
MALIGNANT_CLASSES <- c("Epithelial", "Mixed")  # malignant(+mixed)

# Helper functions
clr_transform <- function(counts_mat, pseudo = 0.5) {
  X <- counts_mat + pseudo
  S <- rowSums(X)
  P <- sweep(X, 1, S, "/")
  L <- log(P)
  sweep(L, 1, rowMeans(L), "-")
}
best_kmeans <- function(X, k_range = 2:6, seed = 1) {
  set.seed(seed)
  have_cluster <- requireNamespace("cluster", quietly = TRUE)
  best <- NULL; best_score <- -Inf
  for (k in k_range) {
    km <- kmeans(X, centers = k, nstart = 25, iter.max = 200)
    score <- if (have_cluster && k > 1) {
      sil <- cluster::silhouette(km$cluster, dist(X))
      mean(sil[, "sil_width"], na.rm = TRUE)
    } else km$betweenss / km$totss
    if (score > best_score) { best <- km; best_score <- score }
  }
  best
}
dominant_name <- function(mean_props_row, thresh = 0.4) {
  mx <- max(mean_props_row, na.rm = TRUE)
  if (!is.finite(mx) || mx < thresh) "Mixed" else names(which.max(mean_props_row))
}

# Collect all malignant(+mixed) spots across HPV+ samples
hpv_samples <- samples_metadata %>% filter(HPV == "HPVpos") %>% pull(Sample)
prods <- c("E1","E2","E5","E6","E7","L1","L2")

spot_prod_list <- map(hpv_samples, function(samp) {
  site_id <- sample_parameters %>% filter(Sample == samp) %>% pull(Site) %>% unique() %>% .[[1]]
  data_dir <- paste(here("data"), site_id, samp, sep = "/")
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

  s_meta <- metadata %>%
    filter(Sample == samp, EpiStroma %in% MALIGNANT_CLASSES, !is.na(Subclone)) %>%
    dplyr::select(SpotID, MPid, Sample)
  
  spots_use <- intersect(colnames(m), s_meta$SpotID)
  if (!length(spots_use)) return(NULL)
  m <- m[, spots_use, drop = FALSE]
  s_meta <- s_meta %>% filter(SpotID %in% spots_use)
  
  # build product counts per spot by summing across any strain (e.g., hpv16__E6, hpv33__E6 → "E6")
  prod_counts <- map_dfc(prods, function(p) {
    rows <- grep(paste0("^hpv[^_]+__", p, "$"), rownames(m))
    if (length(rows) == 0) rep(0, length(spots_use)) else colSums(m[rows, , drop = FALSE])
  })
  names(prod_counts) <- prods
  
  tibble(
    Sample = samp,
    SpotID = spots_use,
    MPid   = s_meta$MPid[match(spots_use, s_meta$SpotID)]
  ) %>% bind_cols(prod_counts)
})

spot_prod_df <- bind_rows(spot_prod_list)
stopifnot(nrow(spot_prod_df) > 0)

# Total HPV per spot and HPV+ flag
spot_prod_df <- spot_prod_df %>%
  mutate(tot_hpv = rowSums(across(all_of(prods))),
         HPV_pos = tot_hpv >= HPV_DET_THRESHOLD)

# Cluster HPV+ spots on composition
hpvpos_df <- spot_prod_df %>% filter(HPV_pos)

# Drop products that are zero across all HPV+ spots (rare but safe)
keep_cols <- prods[colSums(hpvpos_df[, prods, drop = FALSE]) > 0]
X_counts  <- as.matrix(hpvpos_df[, keep_cols, drop = FALSE])

# CLR (compositional) features + standardize
X_clr  <- clr_transform(X_counts, pseudo = PSEUDOCOUNT)
X_feat <- scale(X_clr)

km <- best_kmeans(X_feat, k_range = K_RANGE, seed = 1)
hpvpos_df$Cluster <- factor(km$cluster, levels = sort(unique(km$cluster)))

# cluster mean compositions (proportions)
cl_means <- as.data.frame(X_counts) %>%
  mutate(.cluster = km$cluster) %>%
  group_by(.cluster) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  mutate(total = rowSums(across(all_of(colnames(X_counts))))) %>%
  mutate(across(all_of(colnames(X_counts)), ~ .x / total)) %>%
  dplyr::select(-total)

# cutoff-based cluster naming: all genes with prop ≥ 0.3 are concatenated
prods_order <- c("E1","E2","E5","E6","E7","L1","L2")  # fixed ordering for names
genes_use   <- intersect(prods_order, colnames(cl_means))
cluster_labels <- apply(cl_means[, genes_use, drop = FALSE], 1, function(r) {
  hits <- genes_use[r[genes_use] >= 0.3]
  if (length(hits) == 0) "Mixed" else paste0(hits, collapse = "")
})
name_map <- setNames(cluster_labels, paste0("C", cl_means$.cluster))

# attach human-readable pattern
hpvpos_df <- hpvpos_df %>%
  mutate(HPV_pattern = name_map[paste0("C", Cluster)])

# label ALL malignant(+mixed) spots
cl_names <- colnames(spot_prod_df)
spot_labels_all <- spot_prod_df %>%
  left_join(hpvpos_df %>% dplyr::select(Sample, SpotID, HPV_pattern), by = c("Sample","SpotID")) %>%
  mutate(HPV_pattern = ifelse(is.na(HPV_pattern), "HPVneg", HPV_pattern),
         HPV_pattern = factor(HPV_pattern,
                              levels = c("HPVneg",
                                         unique(name_map[order(match(cl_names, c("E7","E6","E5","E2","E1","L1","L2","Mixed")))]))))

# Plots

# (1) cluster mean compositions
cl_means_long <- cl_means %>%
  mutate(cluster_name = name_map[paste0("C", .cluster)]) %>%
  dplyr::select(cluster_name, all_of(genes_use)) %>%
  pivot_longer(
    cols = all_of(genes_use),
    names_to = "gene", values_to = "prop"
  ) %>%
  mutate(gene = factor(gene, levels = prods_order))

p_cluster_means <- ggplot(cl_means_long, aes(x = cluster_name, y = prop, fill = gene)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  ggsci::scale_fill_jama(name = "HPV gene") +
  labs(title = "HPV+ cluster mean composition (pooled samples)",
       x = "HPV pattern", y = "Proportion") +
  theme_minimal(12) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggplot2::ggsave(filename = here("results/Paper_Figures/Revision_Figs/HPVpos_Spots_Cluster_Composition.png"), plot = p_cluster_means, width = 8, height = 10, dpi = 300)

# (2) MP composition by HPV pattern (pooled)
non_malig <- c("T_cell","B_cell","Macrophage","Fibroblast",
               "Endothelial","Complement","Skeletal_Muscle","Secretory_Norm")

# get the reference order for states (from your earlier plot); fall back to names(state_cols)
mp_levels <- if (exists("per_group_states")) levels(per_group_states$MPid) else names(state_cols)

per_pattern_states <- spot_labels_all %>%
  dplyr::select(HPV_pattern, MPid) %>%
  table() %>% as.data.frame() %>%
  dplyr::group_by(HPV_pattern) %>%
  dplyr::mutate(Prop = Freq / sum(Freq)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    MPid  = factor(MPid, levels = mp_levels),                   # <- force SAME order
    Group = ifelse(MPid %in% non_malig, "Non-Malignant", "Malignant")
  )

p_mp_by_pattern <- per_pattern_states %>%
  dplyr::group_by(MPid, HPV_pattern) %>%
  dplyr::summarise(
    Malignant       = sum(Prop[Group == "Malignant"]),
    `Non-Malignant` = -sum(Prop[Group == "Non-Malignant"]),
    .groups = "drop"
  ) %>%
  dplyr::group_by(HPV_pattern) %>%
  dplyr::mutate(Ord = sum(`Non-Malignant`)) %>%
  ggplot2::ggplot(ggplot2::aes(x = forcats::fct_reorder(HPV_pattern, Ord),
                               y = Malignant, fill = MPid)) +
  ggplot2::geom_col(position = ggplot2::position_stack(reverse = TRUE)) +  # same stacking
  ggplot2::geom_col(ggplot2::aes(y = `Non-Malignant`), position = "stack") +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
  ggplot2::scale_fill_manual(
    name   = "States",
    values = state_cols[mp_levels],
    breaks = mp_levels,
    limits = mp_levels,
    drop   = FALSE
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  ggplot2::labs(title = "MP composition by HPV pattern (pooled malignant+mixed spots)") +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.title   = ggplot2::element_blank(),
    axis.text.x  = ggplot2::element_text(size = 16),
    plot.title   = ggplot2::element_text(hjust = 0.5, size = 16),
    legend.text  = ggplot2::element_text(size = 14),
    legend.title = ggplot2::element_text(size = 16)
  ) +
  ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 7),
                                               reverse = TRUE))
ggplot2::ggsave(filename = here("results/Paper_Figures/Revision_Figs/MP_Composition_Per_HPVpos_Spots_Cluster.png"), plot = p_mp_by_pattern, width = 8, height = 10, dpi = 300)



# Tcell - Clonality Association (R1C3) -----------------------------------------

# Load per-sample parameter table
sample_parameters <- readRDS(file = here("metadata/per_sample_parameters_df.rds"))
samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))
# Load merged metadata tibble
metadata <- readRDS(file = here("metadata/merged_metadata_Extend.rds"))
# Load merged expression matrix
merged_mat <- readRDS(file = here("results/Generated_Data/merged_expr_mat.rds"))
# Load sample color palette
samp_cols <- setNames(samples_metadata$Color, samples_metadata$Sample)
# Signature for T-cell detection
t_chain_genes <- intersect(c("CD3D", "CD3E", "CD3G", "CD247"), rownames(merged_mat))
t_receptor_genes <- intersect(c("TRAC", "TRBC1", "TRBC2", "TRAT1"), rownames(merged_mat))

### ============ Percentage of T-cells containing spots in multi- vs uni-clonal samples (per-sample) ============
# Get the number of clones that each sample has
sample_clone_number <- as.data.frame(rowSums(as.matrix(table(metadata$Sample, metadata$Subclone)) > 0)) %>% 
  `colnames<-`("N_clones") %>% 
  tibble::rownames_to_column(var = "Sample") %>% 
  dplyr::mutate(HPV = samples_metadata$HPV[match(Sample, samples_metadata$Sample)],
                pEMT_Pattern = samples_metadata$pEMT_Pattern[match(Sample, samples_metadata$Sample)],
                Site = samples_metadata$Site[match(Sample, samples_metadata$Sample)]) %>% 
  dplyr::mutate(Group = dplyr::case_when(N_clones == 3 ~ NA,
                                         N_clones > 3 ~ "Multi-clonal",
                                         N_clones < 3 ~ "Uni-clonal"))

# Detect T-cells based on co-expression of a CD3 chains plus a TCR gene
has_chain <- colSums(merged_mat[t_chain_genes, , drop = FALSE] > 0, na.rm = TRUE) > 1
has_receptor <- colSums(merged_mat[t_receptor_genes, , drop = FALSE] > 0, na.rm = TRUE) > 1
spots_both <- colnames(merged_mat)[has_chain & has_receptor]

# T-cell containing spots abundance per sample
t_cells_meta <- metadata %>% dplyr::filter(Key %in% spots_both)
infiltrating_t_cells <- t_cells_meta %>% dplyr::filter(!Zone %in% c("Filtered_out", "Non_Epithelial")) %>% dplyr::group_by(Sample) %>% count()
total_spots_per_samp <- metadata %>% dplyr::group_by(Sample) %>% count()

plot_df <- full_join(as.data.frame(table(t_cells_meta$Sample)), sample_clone_number, by = c("Var1" = "Sample")) %>% 
  dplyr::mutate(Total_spots = total_spots_per_samp$n[match(Var1, total_spots_per_samp$Sample)],
                Inflit_T = infiltrating_t_cells$n[match(Var1, infiltrating_t_cells$Sample)],
                Inflit_T_prop = Inflit_T / Total_spots,
                T_prop = Freq / Total_spots) %>% 
  dplyr::arrange(desc(T_prop)) %>% 
  dplyr::rename(Sample = Var1)

p <- make_boxdot(plot_df[!is.na(plot_df$Group), ], yvar = "T_prop", x_var = "Group", ylab = " ", binwidth = 1/100, percent_y = TRUE, 
                 dot_palette = samp_cols, box_width = 0.6, dot_mm = 2, expend_mult_y = 0.015, tilt_x_text = TRUE)
p1 <- make_boxdot(plot_df, yvar = "T_prop", x_var = "Site", ylab = "", binwidth = 1/30, percent_y = TRUE, 
                  dot_palette = samp_cols, box_width = 0.6, dot_mm = 2, expend_mult_y = 0, tilt_x_text = TRUE, expand_top = 0.15)
p2 <- make_boxdot(plot_df, yvar = "T_prop", x_var = "HPV", ylab = "", binwidth = 1/50, percent_y = TRUE, 
                  dot_palette = samp_cols, box_width = 0.6, dot_mm = 2, expend_mult_y = 0.01, tilt_x_text = TRUE)
p3 <- make_boxdot(plot_df[!is.na(plot_df$Group), ], yvar = "Inflit_T_prop", x_var = "Group", ylab = "", binwidth = 1/100, percent_y = TRUE, 
                  dot_palette = samp_cols, box_width = 0.6, dot_mm = 2, expend_mult_y = 0.01, tilt_x_text = TRUE)     ## Only tumor infiltrating T-cells
p4 <- make_boxdot(plot_df, yvar = "Inflit_T_prop", x_var = "Site", ylab = "", binwidth = 1/80, percent_y = TRUE, 
                  dot_palette = samp_cols, box_width = 0.6, dot_mm = 2, expend_mult_y = 0.01, tilt_x_text = TRUE, expand_top = 0.15)     ## Only tumor infiltrating T-cells
p5 <- make_boxdot(plot_df, yvar = "Inflit_T_prop", x_var = "HPV", ylab = "", binwidth = 1/80, percent_y = TRUE, 
                  dot_palette = samp_cols, box_width = 0.6, dot_mm = 2, expend_mult_y = 0.008, tilt_x_text = TRUE)     ## Only tumor infiltrating T-cells

final_plot <- plot_grid(p, p1, p2, p3,p4, p5, ncol = 3, nrow = 2, align = "hv", axis = "tblr", rel_heights = rep(1, 2))
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C3_Per-sample_Grouped_T-cell_Abundance_Comparisons_Panel.pdf"), plot = final_plot, dpi = 300, width = 8, height = 6)
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C3_Multiclonal_vs_Uniclonal_Samples_T-cell_Abundance_Comparison.pdf"), plot = p, dpi = 300, width = 2, height = 3)


### ============ T-cell gene expression profile differences between Uni- and Muli-clonal samples ============
# t_cell_genes <- read_csv(file = here("Analysis/Revision_Figs/tcell_genes_expanded.csv"))$gene
t_cell_genes <- read_csv(file = here("aux_data/t_cell_markers_human_Tcells_v2.csv"))$genes
keep_samps <- with(plot_df, Sample[T_prop >= 0.05 & !is.na(Group)])

# 2) All T+ spots (chain + receptor), restricted to kept samples
has_chain    <- colSums(merged_mat[t_chain_genes, , drop = FALSE]    > 0, na.rm = TRUE) > 1
has_receptor <- colSums(merged_mat[t_receptor_genes, , drop = FALSE] > 0, na.rm = TRUE) > 1
t_spots_all  <- colnames(merged_mat)[has_chain & has_receptor]

spot_sample  <- sub("_.*$", "", t_spots_all)
t_spots_keep <- t_spots_all[spot_sample %in% keep_samps]

# 3) T-cell panel only
t_panel <- intersect(t_cell_genes, rownames(merged_mat))
X <- merged_mat[t_panel, t_spots_keep, drop = FALSE]

# 4) Pseudo-bulk per sample (mean of T-spots)
samp_vec <- sub("_.*$", "", colnames(X))
samps    <- unique(samp_vec)

PB <- sapply(samps, function(s) {
  idx <- samp_vec == s
  if (sum(idx) == 1L) X[, idx, drop = TRUE] else rowMeans(X[, idx, drop = FALSE], na.rm = TRUE)
}, simplify = "matrix")

# 5) Keep only Uni/Multi samples & set group
grp_by_sample <- setNames(sample_clone_number$Group, sample_clone_number$Sample)
grp_vec <- factor(grp_by_sample[colnames(PB)], levels = c("Uni-clonal","Multi-clonal"))
keep    <- !is.na(grp_vec)
PB      <- PB[, keep, drop = FALSE]
grp_vec <- droplevels(grp_vec[keep])

# 6) Compute Δ and mean (NO anchors)
mu_uni    <- rowMeans(PB[, grp_vec == "Uni-clonal",  drop = FALSE], na.rm = TRUE)
mu_multi  <- rowMeans(PB[, grp_vec == "Multi-clonal", drop = FALSE], na.rm = TRUE)
delta     <- mu_multi - mu_uni                                # y-axis
x_mean    <- rowMeans(PB, na.rm = TRUE)                       # x-axis (uncentered abundance)

res_gene <- tibble::tibble(
  gene  = rownames(PB),
  delta = delta,
  mean_across_samples = x_mean
) |> dplyr::arrange(dplyr::desc(delta))

# Add statistics
design <- stats::model.matrix(~ grp_vec)
fit    <- limma::lmFit(PB, design) |> limma::eBayes()
tt     <- limma::topTable(fit, coef = "grp_vecMulti-clonal",
                          number = Inf, sort.by = "none")
res_gene <- res_gene |>
  dplyr::left_join(
    tibble::tibble(
      gene   = rownames(tt),
      p_value = tt$P.Value,
      padj    = tt$adj.P.Val,
      logFC   = tt$logFC        # should closely match your `delta`
    ),
    by = "gene"
  )

# Δ vs mean plot (p-value free)
top_multi <- res_gene |>
  dplyr::filter(!is.na(delta)) |>
  dplyr::arrange(dplyr::desc(delta)) |>
  dplyr::slice_head(n = 10) |>
  dplyr::pull(gene)

top_uni <- res_gene |>
  dplyr::filter(!is.na(delta)) |>
  dplyr::arrange(delta) |>
  dplyr::slice_head(n = 10) |>
  dplyr::pull(gene)

res_gene <- res_gene |>
  dplyr::mutate(flag_top10 = gene %in% c(top_multi, top_uni))

p_delta_mean <- ggplot2::ggplot(
  dplyr::filter(res_gene, !is.na(mean_across_samples), !is.na(delta)),
  ggplot2::aes(x = mean_across_samples, y = delta)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::geom_point(ggplot2::aes(color = flag_top10), size = 1.9, alpha = 0.6) +
  ggplot2::scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "red"),
                              guide = "none") +
  ggrepel::geom_text_repel(
    data = dplyr::filter(res_gene, flag_top10),
    ggplot2::aes(label = gene),
    min.segment.length = 0, box.padding = 0.25, point.padding = 0.3,
    size = 3, max.overlaps = Inf
  ) +
  ggplot2::labs(
    x = "Mean Expression Across Samples (log2CPM)",
    y = "Mean Expression (Multi-clonal − Uni-clonal)"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  ggplot2::theme_classic()
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C3_Pseudobulk_Per-sample_Grouped_T-cell_Genes_Differntial_Expression_v2.pdf"), plot = p_delta_mean, dpi = 300, width = 6, height = 6)

# Another version of the plot with statistical significance encoded as point shape
res_gene2 <- res_gene %>%
  mutate(point_type = case_when(
    flag_top10 & padj < 0.10 ~ "Top-10 (FDR < 0.1)",
    flag_top10               ~ "Top-10",
    TRUE                     ~ "Other genes"
  )) %>%
  mutate(point_type = factor(point_type,
                             levels = c("Other genes","Top-10","Top-10 (FDR < 0.1)")
  ))

p_delta_mean <- ggplot(
  filter(res_gene2, !is.na(mean_across_samples), !is.na(delta)),
  aes(x = mean_across_samples, y = delta,
      color = point_type, shape = point_type, size = point_type)
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = 0.85) +
  ggrepel::geom_text_repel(
    data = filter(res_gene2, flag_top10),
    aes(label = gene),
    min.segment.length = 0, box.padding = 0.25, point.padding = 0.3,
    size = 3, max.overlaps = Inf, colour = "black", show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Other genes"        = "grey70",
    "Top-10"             = "#d6604d",
    "Top-10 (FDR < 0.1)" = "#d6604d"
  )) +
  scale_shape_manual(values = c(
    "Other genes"        = 16,  # circle
    "Top-10"             = 16,  # circle
    "Top-10 (FDR < 0.1)" = 8    # star
  )) +
  # bump sizes for the two highlighted groups
  scale_size_manual(values = c(
    "Other genes"        = 1.9,
    "Top-10"             = 2.6,
    "Top-10 (FDR < 0.1)" = 3.0
  ), guide = "legend") +
  labs(
    x = "Mean Expression Across Samples (log2CPM)",
    y = "Mean Expression (Multi-clonal − Uni-clonal)",
    color = NULL, shape = NULL, size = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  theme_classic() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.key = element_blank()
  ) +
  guides(
    # make legend symbols a touch larger & opaque
    color = guide_legend(override.aes = list(alpha = 1)),
    shape = guide_legend(override.aes = list(alpha = 1)),
    size  = guide_legend(override.aes = list(alpha = 1))
  )
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C3_Pseudobulk_Per-sample_Grouped_T-cell_Genes_Differntial_Expression_v3.pdf"), plot = p_delta_mean, dpi = 300, width = 6, height = 6)

# Plot T-cell subtype enrichment
counts_inter <- tcell_enrichment_counts(
  res_gene2, tcell_markers,
  top_n = 10, use_labels = TRUE,
  count_only_sig = FALSE,
  tag = "Inter-tumor"
)
scores_inter <- tcell_enrichment_scores(res_gene2, tcell_markers, top_n = 10, use_labels = TRUE)
scores_inter %>% group_by(direction) %>% arrange(desc(z), .by_group = TRUE) %>% slice_head(n = 5)

enriched_inter <- plot_tcell_enrichment(counts_inter, mode = "auto")
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C3_Pseudobulk_Per-sample_Grouped_T-cell_DEGs_Enrichment.pdf"), plot = enriched_inter, dpi = 300, width = 6, height = 6)
enriched_inter <- plot_tcell_enrichment(scores_inter, mode = "auto", pad_units = 0.1, fill_limits = c(0, 3), fill_breaks = 0:3)
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C3_Pseudobulk_Per-sample_Grouped_T-cell_DEGs_Enrichment_v2.pdf"), plot = enriched_inter, dpi = 300, width = 6, height = 6)



### ============ Intra-tumor multi- vs uni-clonal regions ============
# Subclones color palette
clone_cols <- setNames(c("black", "#0c2a50ff", "#BC80BD", "#88CCEE", "#FCCDE5", "#33A02C", "#FFED6F", "#FE9929", "#E31A1C", "#D9D9D9", "#B15928"),
                       c("Non_Malignant", "Mixed", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
# Load T cell signatures
tcell_markers <- list(
  General_T_identity = c(
    "CD3D","CD3E","CD3G","CD247","TRAC","TRBC1","TRBC2",
    "CD2","CD5","PTPRC","LCK","ZAP70","LAT","ITK","TRAT1","BCL11B","THEMIS"
  ),
  Lineage_CD4 = c("CD4"),
  Lineage_CD8 = c("CD8A","CD8B"),
  Naive_CentralMemory = c("CCR7","SELL","IL7R","TCF7","LEF1","BACH2","ID3"),
  Activation_Proliferation = c(
    "CD69","CD38","HLA-DRA","MKI67",
    "CD27","CD28","ICOS","CD40LG","TNFRSF9","TNFRSF4","TNFRSF18","TNFRSF25"
  ),
  Cytotoxic_Effector_CD8 = c(
    "PRF1","GZMB","GZMH","GZMA","GZMK","GNLY","NKG7","FASLG",
    "CCL5","CCL4","CCL3","XCL1","IFNG","TNF","KLRG1","B3GAT1"
  ),
  Treg = c("FOXP3","IL2RA","CTLA4","IKZF2","CCR8","TGFB1","IL10"),
  Th1 = c("TBX21","IFNG","CXCR3","CCR5"),
  Th2 = c("GATA3","IL4","IL5","IL13","CCR4"),
  Th17 = c("RORC","IL17A","IL17F","IL22","CCR6"),
  Tfh = c("BCL6","CXCR5","IL21","PDCD1","ICOS","CXCL13"),
  Tissue_Resident_TRM = c("ITGAE","ITGA1","CXCR6","ZNF683","ITGB7"),
  Exhaustion_Dysfunction = c(
    "PDCD1","LAG3","TIGIT","HAVCR2","CTLA4","BTLA","CD244","CD160","CD96","TOX","PRDM1","EOMES"
  )
)

meta_with_region <- list()
region_t_abund <- list()
de_t_genes <- list()

select_samp <- c("P5710", "P5624", "P5217", "P5684")

for(samp in select_samp) {
  message(paste0("Processing sample: ", samp))
  
  # Load expression matrix
  data_dir <- paste(here("data"), sample_parameters$Site[sample_parameters$Sample == samp], samp, sep = "/")
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
  # Load spatial image object
  spatial_image <- load_spatial_image(paste0(data_dir, "/spatial"))
  
  # Create Seurat object and manually annotate regions of interest
  samples <- paste(data_dir, "filtered_feature_bc_matrix.h5", sep = "/")
  imgs <- paste(data_dir, "spatial/tissue_hires_image.png", sep = "/")
  spotfiles <- paste(data_dir, "spatial/tissue_positions_list.csv", sep = "/")
  json <- paste(data_dir, "spatial/scalefactors_json.json", sep = "/")
  
  infoTable <- data.frame(samples, imgs, spotfiles, json)
  se <- STutility::InputFromTable(infoTable)
  se <- STutility::LoadImages(se, time.resolve = FALSE)
  # NOTE: ManualAnnotation() requires interactive input and cannot run unattended -- the result is saved below so this step only needs to be run once
  se <- STutility::ManualAnnotation(se, type = "raw", res = 1000)
  rownames(se@meta.data) <- scalop::substri(rownames(se@meta.data), sep = "_", pos = 1)

  samp_meta <- metadata[metadata$Sample == samp, ]
  samp_meta <- samp_meta %>% dplyr::mutate(Region = se@meta.data$labels[match(SpotID, rownames(se@meta.data))])
  ## OR ##
  # meta_with_region <- readRDS(file = here("results/Generated_Data/per-sample_clonal_regions_metadata.rds"))
  # samp_meta <- meta_with_region[[samp]]
  order_cells_full <- readRDS(file = paste0(here("results/Per_Sample_Results/"), samp, "/", samp, "_CNA_Subclones.rds"))
  samp_meta$Subclone <- ifelse(samp_meta$SpotID %in% Unlist(order_cells_full), names(Unlist(order_cells_full)[match(samp_meta$SpotID, Unlist(order_cells_full))]), NA)
  
  # Plot the selected region on top of the spatial distribution map of the subclones
  scale_used <- "lowres"
  coords <- .get_coords(spatial_image, scale = scale_used)
  df_all <- tibble::tibble(SpotID = rownames(coords),
                           y      = coords[,1],
                           x      = coords[,2]) %>% 
    dplyr::left_join(samp_meta %>% dplyr::select(SpotID, Region), by = "SpotID")
  df_sel <- df_all %>% dplyr::filter(Region %in% c("Multi-clonal", "Uni-clonal")) %>% droplevels()
  
  pts_sel <- sf::st_as_sf(df_sel, coords = c("x","y"), crs = NA)
  pts_all <- sf::st_as_sf(df_all, coords = c("x","y"), crs = NA)
  
  xy <- as.matrix(df_all[, c("x","y")])
  spacing <- tryCatch(FNN::get.knn(xy, k = 2)$nn.dist[,2] |> median(na.rm=TRUE), 
                      error = function(e) {
                        idx <- if (nrow(xy) > 3000) sample(nrow(xy), 3000) else seq_len(nrow(xy))
                        median(sapply(idx, function(i){
                          d <- sqrt((xy[i,1]-xy[,1])^2 + (xy[i,2]-xy[,2])^2)
                          sort(d[d>0], partial=1)[1]
                        }), na.rm=TRUE)
                      })
  r_mult <- 0.6
  r <- r_mult * spacing
  outline_off <- 0.18 * spacing
  
  regions_poly <- pts_sel %>%
    dplyr::group_by(Region) %>%
    dplyr::summarise(geometry = st_union(st_buffer(geometry, r)), .groups = "drop") %>%
    sf::st_make_valid()
  
  tissue_poly <- sf::st_union(sf::st_buffer(pts_all$geometry, r)) %>% sf::st_make_valid()
  tissue_interior <- sf::st_buffer(tissue_poly, -0.35*r)
  regions_interior <- sf::st_intersection(regions_poly, tissue_interior)
  regions_interior <- regions_interior[!st_is_empty(regions_interior), ]
  
  bnd_ml <- sf::st_cast(sf::st_boundary(regions_interior), "MULTILINESTRING")
  co <- as.data.frame(sf::st_coordinates(bnd_ml))
  bnd_df <- transform(co,
                      Region = regions_interior$Region[L1],
                      grp = interaction(L1, L2, drop = TRUE),
                      x = X, y = Y)
  yr <- range(df_all$y, na.rm = TRUE)
  bnd_df$y <- sum(yr) - bnd_df$y
  
  p <- plot_spatial_features(metadata = samp_meta, image_obj = spatial_image, color_by = "Subclone", cols = clone_cols, pt.size = 1.5, image.alpha = 0, rm_spots = NULL, na.value = "gray", stroke = 0, scale = scale_used) +
    guides(fill = guide_legend(override.aes = list(size = 5)))
  
  p <- p + ggplot2::geom_path(data = bnd_df, aes(x = x, y = y, group = grp), inherit.aes = FALSE, linewidth = 1, color = "red")
  ggplot2::ggsave(filename = paste0(here("results/Paper_Figures/Revision_Figs/"), samp, "_Marking_Multi_vs_Uni_Clonal_Regions.png"),
                  plot = p, width = 16, height = 14, units = "cm")
  
  # T-cell abundance and differental gene expression between the regions (multi vs uni clonal)
  t_cell_genes <- tcell_markers %>% unlist %>% unique
  t_chain_genes <- intersect(c("CD3D", "CD3E", "CD3G", "CD247"), rownames(m))
  t_receptor_genes <- intersect(c("TRAC", "TRBC1", "TRBC2", "TRAT1"), rownames(m))
  
  # Detect T-cells based on co-expression of a CD3 chains plus a TCR gene
  has_chain <- colSums(m[t_chain_genes, , drop = FALSE] > 0, na.rm = TRUE) > 0
  has_receptor <- colSums(m[t_receptor_genes, , drop = FALSE] > 0, na.rm = TRUE) > 0
  spots_both <- colnames(m)[has_chain & has_receptor]
  
  # T-cell containing spots abundance per sample
  region_meta <- samp_meta %>% 
    dplyr::filter(Region %in% c("Multi-clonal", "Uni-clonal") & !Zone %in% c("Filtered_out", "Non_Epithelial")) %>% 
    dplyr::group_by(Region) %>% dplyr::summarise(n_spots = n())
  t_cells_meta <- samp_meta %>% dplyr::filter(SpotID %in% spots_both) %>% 
    dplyr::filter(Region %in% c("Multi-clonal", "Uni-clonal") & !Zone %in% c("Filtered_out", "Non_Epithelial")) %>% 
    dplyr::group_by(Region) %>% dplyr::summarise(n_spots = n())
  regional_t_cell_props <- dplyr::left_join(region_meta, t_cells_meta, by = "Region") %>% 
    dplyr::mutate(T_prop = n_spots.y / n_spots.x)
  
  region_t_abund[[samp]] <- tibble(Sample = samp,
                                   Multi = regional_t_cell_props$T_prop[regional_t_cell_props$Region == "Multi-clonal"],
                                   Uni = regional_t_cell_props$T_prop[regional_t_cell_props$Region == "Uni-clonal"])
  meta_with_region[[samp]] <- samp_meta
  
  # T-cell genes differential expression
  filt_m <- filter_10x(m, complexity_thresh = NULL, percent_mt_thresh = NULL, gene_filt_method = "none", merge_dup_genes = FALSE,
                       centre = "none", log = TRUE, raw_to_CPM = TRUE)
  
  # meta_with_region <- readRDS(file = here("results/Generated_Data/per-sample_clonal_regions_metadata.rds"))
  samp_meta <- meta_with_region[[samp]]
  roi <- samp_meta %>%
    dplyr::filter(SpotID %in% spots_both,
                  Region %in% c("Multi-clonal","Multi-clonal","Uni-clonal")) %>%
    dplyr::select(SpotID, Region)
  spots_roi <- roi$SpotID
  X <- filt_m[, spots_roi, drop = FALSE]
  region <- factor(roi$Region, levels = c("Uni-clonal","Multi-clonal"))
  genes_interest <- intersect(t_cell_genes, rownames(X))
  # Filter genes with mean expression lower than that of the lowest (non-zero) mean expression of any of the general T-cell genes
  avg_general_t_expr <- rowMeans(X[tcell_markers$General_T_identity, ])
  filt_gene_cut <- min(avg_general_t_expr[avg_general_t_expr > 0])
  avg_t_expr <- rowMeans(X[genes_interest, ])
  genes_interest <- names(avg_t_expr[avg_t_expr > filt_gene_cut])
  
  t_def <- tcell_markers$General_T_identity
  t_score <- colMeans(X[t_def, , drop = FALSE])
  n_def <- length(t_def)
  sub_mat <- matrix(t_score, nrow = length(genes_interest), ncol = length(spots_roi), byrow = TRUE,
                    dimnames = list(genes_interest, spots_roi))
  
  for (g in intersect(genes_interest, t_def)) {
    sub_mat[g, ] <- (t_score * n_def - X[g, ])/(n_def - 1)
  }
  
  X_adj <- X[genes_interest, , drop = FALSE] - sub_mat
  idx <- split(seq_along(spots_roi), region)
  diff_adj <- rowMeans(X_adj[, idx[["Multi-clonal"]], drop = FALSE]) -
    rowMeans(X_adj[, idx[["Uni-clonal"]],  drop = FALSE])
  
  pvals <- sapply(genes_interest, function(g)
    wilcox.test(X_adj[g, idx[["Multi-clonal"]]],
                X_adj[g, idx[["Uni-clonal"]]], exact = FALSE)$p.value)
  
  res_tbl <- data.frame(
    gene = genes_interest,
    diff_adj = diff_adj[genes_interest],
    p = pvals,
    FDR = p.adjust(pvals, "BH"),
    row.names = NULL
  ) |> dplyr::arrange(FDR)
  
  de_t_genes[[samp]] <- res_tbl
}

saveRDS(region_t_abund, file = here("results/Generated_Data/per-sample_clonal_regions_T-cell_abundance.rds"))
saveRDS(meta_with_region, file = here("results/Generated_Data/per-sample_clonal_regions_metadata.rds"))
saveRDS(de_t_genes, file = here("results/Generated_Data/per-sample_clonal_regions_T-cell_DEGs.rds"))
test <- dplyr::bind_rows(region_t_abund)
plot_df <- test %>%
  tidyr::pivot_longer(
    cols = -Sample,              # or c(Multi, Uni)
    names_to  = "Group",
    values_to = "T_prop"
  ) %>%
  dplyr::mutate(Group = paste0(Group, "-clonal")) %>%
  dplyr::arrange(Sample, Group)
p <- make_boxdot(plot_df[!is.na(plot_df$Group), ], yvar = "T_prop", x_var = "Group", ylab = "Fraction of Spots Containing T-cells (%)", binwidth = 1/100, percent_y = TRUE, 
                 dot_palette = samp_cols, box_width = 0.6, dot_mm = 2, expend_mult_y = 0.015, tilt_x_text = TRUE)
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C3_Per-sample_Multiclonal_vs_Uniclonal_Region_T-cell_Abundance_Comparison.pdf"), plot = p, dpi = 300, width = 2, height = 3)

# Gene-level plot
de_t_genes <- readRDS(file = here("results/Generated_Data/per-sample_clonal_regions_T-cell_DEGs.rds"))
long <- imap_dfr(de_t_genes, ~{
  .x %>%
    transmute(gene = toupper(gene), diff_adj, Sample = .y)
}) %>%
  group_by(Sample, gene) %>%                    # just in case of duplicates
  summarise(diff_adj = mean(diff_adj, na.rm = TRUE), .groups = "drop")

tt_one <- function(x, mu = 0, conf.level = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- mean(x)
  s <- stats::sd(x)
  if (n < 2 || isTRUE(all.equal(s, 0))) {
    return(list(n = n, mean = m, sd = s, t = NA_real_, p = NA_real_,
                ci_lo = NA_real_, ci_hi = NA_real_, d = NA_real_, g = NA_real_))
  }
  tt <- t.test(x, mu = mu, conf.level = conf.level)
  d  <- m / s
  J  <- 1 - 3 / (4*(n - 1) - 1)                 # Hedges' small-sample corr
  g  <- J * d
  list(n = n, mean = m, sd = s,
       t = unname(tt$statistic), p = unname(tt$p.value),
       ci_lo = unname(tt$conf.int[1]), ci_hi = unname(tt$conf.int[2]),
       d = d, g = g)
}

# t-test per gene
res_t <- long %>%
  group_by(gene) %>%
  summarise(as_tibble(tt_one(diff_adj)), .groups = "drop") %>%
  mutate(FDR = p.adjust(p, "BH"))

# consistency diagnostics (on available samples per gene)
cons <- long %>%
  group_by(gene) %>%
  summarise(
    n_nonmiss = sum(is.finite(diff_adj)),
    n_pos     = sum(diff_adj > 0, na.rm = TRUE),
    n_neg     = sum(diff_adj < 0, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    n_consistent   = pmax(n_pos, n_neg),                 # # samples in majority direction
    frac_consistent = n_consistent / n_nonmiss
  )

de_meta_ttest <- res_t %>%
  left_join(cons, by = "gene") %>%
  arrange(FDR)

de_plot <- de_meta_ttest %>%
  transmute(
    Gene = gene,
    x    = mean,                        # meta log2FC (Multi − Uni)
    Pval = p,                           # raw p-value
    y    = -log10(pmax(Pval, 1e-300)),
    n_consistent, n_nonmiss
  ) %>%
  mutate(
    direction   = ifelse(x >= 0, "Up in multiclonal", "Up in uniclonal"),
    n_sig_dir_f = factor(pmin(n_consistent, 4L), levels = 0:4)
  ) %>%
  filter(!Gene %in% tcell_markers$General_T_identity) %>% 
  group_by(direction) %>%
  arrange(Pval, desc(abs(x)), .by_group = TRUE) %>%  # rank by p; tie-break by |LFC|
  mutate(
    rank_dir   = row_number(),
    ColorGroup = case_when(
      direction == "Up in multiclonal" & rank_dir <= 10 ~ "Up in multiclonal",
      direction == "Up in uniclonal"   & rank_dir <= 10 ~ "Up in uniclonal",
      TRUE ~ "Other"
    ),
    Label = ifelse(rank_dir <= 10, Gene, NA_character_)
  ) %>%
  ungroup()

sz_fac <- 0.65

# size map (numeric named vector), keep "0" for plotting but exclude it from legend
size_map_base   <- c("0"=0.15, "1"=1.6, "2"=3.0, "3"=5.2, "4"=8.0)
lvl_all         <- levels(de_plot$n_sig_dir_f)                 # c("0","1","2","3","4")
size_map_scaled <- (size_map_base * sz_fac)[lvl_all]

# star overlay scaled too
star_size   <- 3.2 * sz_fac
star_stroke <- 0.6 * sz_fac

# axes paddings (as you had)
xr   <- range(de_plot$x, na.rm = TRUE)
yr   <- range(de_plot$y, na.rm = TRUE)
xpad <- 0.04 * diff(xr)
ypad <- 0.01 * diff(yr)
upper <- yr[2] + 0.12 * diff(yr)

# label data (top-10 per direction)
labs_all <- subset(de_plot, rank_dir <= 10)

# --- key tweaks for clearer connector lines ---
pt_pad   <- grid::unit(0.5, "mm")   # very small: segment starts at the point
box_pad  <- grid::unit(3,   "mm")   # modest label spacing
nudgeY   <- 0.035 * diff(yr)        # a bit smaller upward nudge
seg_col  <- "black"
seg_size <- 0.4
seg_alpha<- 0.8

p <- ggplot(de_plot, aes(x = x, y = y)) +
  # base points (size encodes consistency)
  geom_point(aes(colour = ColorGroup, size = n_sig_dir_f), alpha = 0.95) +
  # star overlay (raw p < 0.05)
  geom_point(
    data = dplyr::filter(de_plot, Pval < 0.05),
    aes(shape = "p < 0.05"),
    size = star_size, stroke = star_stroke, colour = "black"
  ) +
  coord_cartesian(
    xlim = c(xr[1] - xpad, xr[2] + xpad),
    ylim = c(-ypad, upper),
    expand = FALSE, clip = "off"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_color_manual(
    name = "Direction (top 10 only)",
    values = c("Other" = "grey70",
               "Up in uniclonal" = "#4393c3",
               "Up in multiclonal" = "#d6604d")
  ) +
  # hide "0" from the legend by not listing it in breaks
  scale_size_manual(
    name   = "Consistency of\nthe direction",
    values = size_map_scaled,
    breaks = c("1","2","3","4")
  ) +
  scale_shape_manual(
    name   = "Per-gene p-value",
    values = c("p < 0.05" = 8)
  ) +
  labs(
    x = expression("meta log"[2]*"FC  (Multi − Uni)"),
    y = expression("-log"[10]*" p")
  ) +
  ggrepel::geom_text_repel(
    data  = labs_all,
    aes(label = Label),
    color = "black",
    size  = 3,                        # fixed label size as requested
    point.padding    = pt_pad,        # << anchor near the point center
    box.padding      = box_pad,
    nudge_y          = nudgeY,
    force            = 2.0,
    force_pull       = 0.2,
    min.segment.length = 0,
    segment.size     = seg_size,
    segment.alpha    = seg_alpha,
    segment.color    = seg_col,
    seed = 1, show.legend = FALSE
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position   = "right",
    legend.key.height = grid::unit(1.2, "lines")
  )
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C3_Per-sample_Multiclonal_vs_Uniclonal_Region_T-cell_Genes_Differntial_Expression.pdf"), plot = p, dpi = 300, width = 6, height = 6)

# Plot T-cell subtype enrichment
counts_intra <- tcell_enrichment_counts(
  de_plot, tcell_markers,
  top_n = 10, use_labels = TRUE,
  count_only_sig = FALSE,
  tag = "Intra-tumor"
)
scores_intra <- tcell_enrichment_scores(de_plot, tcell_markers, top_n = 10, use_labels = TRUE)

enriched_intra <- plot_tcell_enrichment(counts_intra, mode = "auto")
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C3_Per-sample_Multiclonal_vs_Uniclonal_Region_T-cell_DEGs_Enrichment.pdf"), plot = enriched_intra, dpi = 300, width = 6, height = 6)
enriched_intra <- plot_tcell_enrichment(scores_intra, mode = "auto", pad_units = 0.1, fill_limits = c(0, 3), fill_breaks = 0:3)
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C3_Per-sample_Multiclonal_vs_Uniclonal_Region_T-cell_DEGs_Enrichment_v2.pdf"), plot = enriched_intra, dpi = 300, width = 6, height = 6)


# HPV Status Affect MP Abundance (R1C7) ----------------------------------------

# Load state colors 
state_cols <- read.delim(file = here("aux_data/state_col_tab_extend.tsv"), sep = "\t", header = TRUE)
state_cols <- setNames(state_cols$V2, state_cols$V1)
# Load per-sample parameter table
sample_parameters <- readRDS(file = here("metadata/per_sample_parameters_df.rds"))
# Load sample level metadata
samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))
# Load merged metadata tibble
metadata <- readRDS(file = here("metadata/merged_metadata_Extend.rds"))

# Create the same figure as Fig.3c but controlled for anatomical subsite variablility ===> compare HPV+ and HPV- samples derived from the base of tongue.
sum_hpv_status_MPs_prop <- metadata %>% dplyr::mutate(HPVstat = as.factor(samples_metadata$HPV[match(Sample, samples_metadata$Sample)])) %>% 
  # dplyr::filter(!is.na(MPid) & MPid != "Skeletal_Muscle" & Sample %in% c("P5903", "P5767", "P5764", "P5735", "P5692", "P5756", "P5751")) %>% droplevels %>%
  dplyr::filter(!is.na(MPid) & MPid != "Skeletal_Muscle" & Sample %in% samples_metadata$Sample[samples_metadata$Site == "Oropharynx"]) %>% droplevels %>%
  dplyr::select(MPid, HPVstat) %>% table() %>% as.data.frame() %>% 
  group_by(HPVstat) %>% mutate(Prop = Freq / sum(Freq)) %>% ungroup
p <- sum_hpv_status_MPs_prop %>% 
  dplyr::mutate(Group = ifelse(MPid %in% c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement", "Secretory_Norm"), "Non-Malignant", "Malignant")) %>% 
  dplyr::group_by(MPid, HPVstat) %>% 
  summarize(Malignant = sum(Prop[Group == "Malignant"]), `Non-Malignant` = -sum(Prop[Group == "Non-Malignant"])) %>%
  dplyr::group_by(HPVstat) %>% 
  dplyr::mutate(Ord = sum(`Non-Malignant`)) %>% 
  ggplot(aes(x = fct_reorder(HPVstat, Ord), y = Malignant, fill = MPid)) + 
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_col(aes(y = `Non-Malignant`), position = "stack") +
  geom_hline(aes(yintercept = 0)) + 
  scale_fill_manual(name = "States", values = state_cols[levels(sum_hpv_status_MPs_prop$MPid)]) + scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) + 
  scale_x_discrete(expand = c(0.3, 0.2)) +
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 16),
                          legend.text = element_text(size = 14), legend.title = element_text(size = 16), aspect.ratio = 1.4) + 
  guides(fill = guide_legend(override.aes = list(size = 7), reverse = TRUE)) 
ggplot2::ggsave(filename = here("results/Paper_Figures/Revision_Figs/HPVpos_vs_HPVneg_MP_Proportions_ONLY_OPSCC.png"), plot = p, width = 10, height = 8)
ggplot2::ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C7_HPVpos_vs_HPVneg_MP_Proportions_ONLY_OPSCC.pdf"), plot = p, device = "pdf", dpi = 300, width = 10, height = 8)

# Draw statistics
MPs_per_hpv_stats <- metadata %>% dplyr::mutate(HPVstat = as.factor(samples_metadata$HPV[match(Sample, samples_metadata$Sample)])) %>% 
  # dplyr::filter(!is.na(MPid) & MPid != "Skeletal_Muscle" & Sample %in% c("P5903", "P5767", "P5764", "P5735", "P5692", "P5756", "P5751")) %>% droplevels
  dplyr::filter(!is.na(MPid) & MPid != "Skeletal_Muscle" & Sample %in% samples_metadata$Sample[samples_metadata$Site == "Oropharynx"]) %>% droplevels
prop_stat_tab <- test_state_prop_diff(MPs_per_hpv_stats$MPid, MPs_per_hpv_stats$Sample, MPs_per_hpv_stats$HPVstat, transform = NULL)
write.csv(prop_stat_tab, file = here("results/Generated_Data/State_Proportion_BOT-HPVneg_vs_HPVpos_Ttest_Results_Extend.csv"))




# pEMT-core vs HPVneg-pEMT-edge to Control for HPV Status (R1C8) ---------------

metadata <- readRDS(file = here("metadata/merged_metadata_Extend.rds"))

# Load Node attributes tables
node_attr_tab <- read.csv(file = here("results/Generated_Data/Node_Abundance_Neighborhood_Network_corrected.csv"))
node_attr_tab <- node_attr_tab %>% dplyr::select(State, pEMTcore, HPVneg_pEMTmargin, Coherence_Core, Coherence_HPVnegMargin) %>% 
  dplyr::mutate(State = stringr::str_replace_all(State, "TNF_Signaling", "Inflammatory"))

# Compute differential abundace and coherence and plot
plot_hm <- node_attr_tab %>% dplyr::reframe(State = State,
                                            Abundance = pEMTcore - HPVneg_pEMTmargin,
                                            Coherence = Coherence_Core - Coherence_HPVnegMargin) %>% 
  dplyr::mutate(Abundance = Abundance * 100) %>% 
  reshape2::melt(.) %>% dplyr::mutate(State = factor(State, levels = c("T_cell", "B_cell", "Macrophage", "Complement", "Fibroblast", "Endothelial", "Cell_Cycle", "pEMT", "Epithelial", "Hypoxia", "Inflammatory", "Senescence")))
p1 <- gmap(plot_hm[plot_hm$variable == "Abundance", ], x = State, y = variable, fill = value, limits = c(-7, 7), midpoint = c(-2, 2), geom = "tile", tile.size = 0,
           y.name = "", legend.height = 0.6, legend.width = 8, x.num = F, y.num = F, minimal = TRUE, legend.title = "Percent Difference") +
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 0.5, axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 14), title = element_text(size = 14), legend.direction = "horizontal") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  facet_grid(. ~ State, scales = "free") + theme(strip.text.x = element_text(size = 14, face = "bold"))
ggplot2::ggsave(filename = here("results/Paper_Figures/Revision_Figs/pEMT-core_vs_HPVneg-pEMT-edge_Differential_Abundance_Heattmap.pdf"), plot = p1, device = "pdf", dpi = 300, width = 20, height = 8)

p2 <- gmap(plot_hm[plot_hm$variable == "Coherence", ], x = State, y = variable, fill = value, limits = c(-0.3, 0.3), midpoint = c(-0.085, 0.085), geom = "tile", tile.size = 0,
           y.name = "", legend.height = 0.6, legend.width = 8, x.num = F, y.num = F, minimal = TRUE, legend.title = "Score Difference") +
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 0.5, axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 14), title = element_text(size = 14), legend.direction = "horizontal") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  facet_grid(. ~ State, scales = "free") + theme(strip.text.x = element_text(size = 14, face = "bold"))
ggplot2::ggsave(filename = here("results/Paper_Figures/Revision_Figs/pEMT-core_vs_HPVneg-pEMT-edge_Differential_Coherence_Heattmap.pdf"), plot = p2, device = "pdf", dpi = 300, width = 20, height = 8)


# Create a pEMT-core vs HPVneg-pEMT-edge network
network_ls <- list(pEMT_core = read.csv(file = here("results/Generated_Data/pEMT_CORE_Neighborhood_Network_corrected.csv")) %>% 
                     dplyr::mutate(Alpha = abs(Effect_Size)) %>% dplyr::rename_with(~paste(., "pEMTcore", sep = "."), 3:ncol(.)),
                   pEMT_margin = read.csv(file = here("results/Generated_Data/HPV_Negative_pEMT_MARGIN_Neighborhood_Network_corrected.csv")) %>% 
                     dplyr::mutate(Alpha = abs(Effect_Size)) %>% dplyr::rename_with(~paste(., "pEMTmargin", sep = "."), 3:ncol(.)))
merged_network <- network_ls %>% reduce(full_join, by = c("Ref_State" = "Ref_State", "Neighb_State" = "Neighb_State")) %>% 
  dplyr::filter(!Ref_State %in% c("LowQ", "Skeletal_Muscle", "Secretory"), !Neighb_State %in% c("LowQ", "Skeletal_Muscle", "Secretory"))
merged_network <- merged_network %>% dplyr::mutate(IDnames = paste0(Ref_State, " (interacts with) ", Neighb_State))

merged_network <- merged_network %>% dplyr::mutate(Group = dplyr::case_when(!is.na(Interaction_Type.pEMTcore) & !is.na(Interaction_Type.pEMTmargin) & Interaction_Type.pEMTcore == Interaction_Type.pEMTmargin ~ "Shared",
                                                                            !is.na(Interaction_Type.pEMTcore) & !is.na(Interaction_Type.pEMTmargin) & Interaction_Type.pEMTcore != Interaction_Type.pEMTmargin & Interaction_Type.pEMTcore == "Drawn" ~ "pEMT_core",
                                                                            !is.na(Interaction_Type.pEMTcore) & !is.na(Interaction_Type.pEMTmargin) & Interaction_Type.pEMTcore != Interaction_Type.pEMTmargin & Interaction_Type.pEMTmargin == "Drawn" ~ "pEMT_margin",
                                                                            !is.na(Interaction_Type.pEMTcore) & is.na(Interaction_Type.pEMTmargin) ~ "pEMT_core",
                                                                            !is.na(Interaction_Type.pEMTmargin) & is.na(Interaction_Type.pEMTcore) ~ "pEMT_margin")) %>% 
  dplyr::mutate(Interaction_Type = dplyr::case_when(Group == "Shared" ~ Interaction_Type.pEMTcore,
                                                    Group == "pEMT_core" ~ Interaction_Type.pEMTcore,
                                                    Group == "pEMT_margin" ~ Interaction_Type.pEMTmargin),
                Filt_n = dplyr::case_when(Group == "Shared" ~ Filt_n.pEMTcore + Filt_n.pEMTmargin,
                                          Group == "pEMT_core" ~ Filt_n.pEMTcore,
                                          Group == "pEMT_margin" ~ Filt_n.pEMTmargin),
                Effect_Size = dplyr::case_when(Group == "Shared" ~ (Effect_Size.pEMTcore + Effect_Size.pEMTmargin) / 2,
                                               Group == "pEMT_core" ~ Effect_Size.pEMTcore,
                                               Group == "pEMT_margin" ~ Effect_Size.pEMTmargin)) %>% 
  dplyr::mutate(Frac_Filt_n = dplyr::case_when(Group == "Shared" ~ Filt_n / 11,
                                               Group == "pEMT_core" ~ Filt_n.pEMTcore / 6,
                                               Group == "pEMT_margin" ~ Filt_n.pEMTmargin / 5),
                Alpha = ifelse(Effect_Size > 0, Effect_Size, 0)) %>% 
  dplyr::select(Ref_State, Neighb_State, IDnames, Group, Interaction_Type, Filt_n, Effect_Size, Frac_Filt_n, Alpha)
write.csv(merged_network, file = here("results/Generated_Data/pEMT-core_vs_HPVneg-pEMT-edge_Neighborhood_Network.csv"), row.names = FALSE, na = '', quote = FALSE)


group_props_df <- as.data.frame(table(metadata$MPid[metadata$Sample %in% c("P5219", "P5624", "P5709", "P5418", "P5659", "P5909",
                                                                           "P5448", "P5707", "P5666", "P5684", "P5756")])) %>% 
  dplyr::mutate(State_Prop = Freq / sum(Freq)) %>% dplyr::select(-Freq) %>% 
  dplyr::rename(State = Var1) %>% 
  dplyr::filter(!State %in% c("LowQ", "Skeletal_Muscle", "Secretory_Malig", "Secretory_Norm"))
# Rescale the proportions to 1 after the removal of "LowQ", "Skeletal_Muscle" and "Secretory" states.
group_props_df <- group_props_df %>% dplyr::mutate(across(where(is.numeric), ~ .x / sum(.x))) %>% dplyr::mutate(State = stringr::str_replace_all(State, "Inflammatory", "TNF_Signaling"))

coherence_scores <- readRDS(file = here("results/Generated_Data/all_samples_MPs_coherence_scores.rds"))
mp_coher_df <- do.call(cbind.data.frame, coherence_scores[c("P5219", "P5624", "P5709", "P5418", "P5659", "P5909",
                                                            "P5448", "P5707", "P5666", "P5684", "P5756")]) %>%
  dplyr::reframe(Coherence = unname(rowMeans(., na.rm = TRUE))) %>% 
  dplyr::mutate(State = sort(names(coherence_scores$P5219)), .before = 1)
group_props_df <- group_props_df %>% dplyr::left_join(mp_coher_df, by = "State")
write.csv(group_props_df, file = here("results/Generated_Data/pEMT-core_vs_HPVneg-pEMT-edge_Node_Abundance_Neighborhood_Network.csv"), row.names = FALSE, na = '', quote = FALSE)



# OSM-OSMR Spatial Analysis (R1C9) ---------------------------------------------

source(here("analysis_scripts/functions/Connectivity_Functions.R"))
# Load per-sample parameter table
sample_parameters <- readRDS(file = here("metadata/per_sample_parameters_df.rds"))
# Load sample level metadata
samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))
# Load metadata
metadata <- readRDS(file = here("metadata/merged_metadata_Extend.rds"))
# Load merged expression matrix
merged_mat <- readRDS(file = here("results/Generated_Data/merged_expr_mat.rds"))
# Load state colors 
state_cols <- read.delim(file = here("aux_data/state_col_tab_extend.tsv"), sep = "\t", header = TRUE)
state_cols <- setNames(state_cols$V2, state_cols$V1)
# Sample color palette
samp_cols <- setNames(samples_metadata$Color, samples_metadata$Sample)

samp_info_tab <- samples_metadata %>% dplyr::filter(pEMT_Pattern != "none")


### ============ Panel A - Fraction of pEMT Spots Having Inf-Resp Neighbor ============
samples <- samp_info_tab %>% dplyr::filter(pEMT_Pattern == "Margin_Enriched") %>% dplyr::pull(Sample)
results <- map_dfr(samples, function(samp) {
  samp_meta  <- filter(metadata, Sample == samp)
  adj_mat <- .neighbors_table_func(samp_meta, spot_class = "MPid", spot_name = "Key")
  pemt_spots <- samp_meta %>% 
    dplyr::filter(!is.na(MPid) & MPid == "pEMT") %>% 
    dplyr::pull(Key)
  pemt_adj_mat <- adj_mat[pemt_spots, , drop = FALSE]
  drop_homotypic <- apply(pemt_adj_mat, 1, function(x) all(is.na(x) | x == "pEMT"))
  pemt_adj_mat <- pemt_adj_mat[!drop_homotypic, , drop = FALSE]
  has_inf_neigh <- apply(pemt_adj_mat, 1, function(nei_labels) {
    any(nei_labels == "Inflammatory", na.rm = TRUE)
  })
  frac_pEMT_with_inf <- mean(has_inf_neigh)
  count_with_inf <- sum(has_inf_neigh)
  total_pEMT <- nrow(pemt_adj_mat)
  
  tibble(Sample = samp,
         total_pEMT_spots = total_pEMT,
         pEMT_with_inf_neighbors = count_with_inf,
         frac_pEMT_with_inf = frac_pEMT_with_inf)
})

samples <- samp_info_tab %>% dplyr::filter(pEMT_Pattern == "Core_Enriched") %>% dplyr::pull(Sample)
results2 <- map_dfr(samples, function(samp) {
  samp_meta  <- filter(metadata, Sample == samp)
  adj_mat <- .neighbors_table_func(samp_meta, spot_class = "MPid", spot_name = "Key")
  pemt_spots <- samp_meta %>% 
    dplyr::filter(!is.na(MPid) & MPid == "pEMT") %>% 
    dplyr::pull(Key)
  pemt_adj_mat <- adj_mat[pemt_spots, , drop = FALSE]
  drop_homotypic <- apply(pemt_adj_mat, 1, function(x) all(is.na(x) | x == "pEMT"))
  pemt_adj_mat <- pemt_adj_mat[!drop_homotypic, , drop = FALSE]
  has_inf_neigh <- apply(pemt_adj_mat, 1, function(nei_labels) {
    any(nei_labels == "Inflammatory", na.rm = TRUE)
  })
  frac_pEMT_with_inf <- mean(has_inf_neigh)
  count_with_inf <- sum(has_inf_neigh)
  total_pEMT <- nrow(pemt_adj_mat)
  
  tibble(Sample = samp,
         total_pEMT_spots = total_pEMT,
         pEMT_with_inf_neighbors = count_with_inf,
         frac_pEMT_with_inf = frac_pEMT_with_inf)
})


res1 <- results  %>% mutate(Set = "Margin_Enriched")
res2 <- results2 %>% mutate(Set = "Core_Enriched")
sample_order <- c(unique(res2$Sample[order(res2$frac_pEMT_with_inf, decreasing = T)]), unique(res1$Sample[order(res1$frac_pEMT_with_inf, decreasing = T)]))

combined <- bind_rows(res1, res2) %>%
  dplyr::mutate(Set = factor(Set, levels = c("Core_Enriched", "Margin_Enriched")),
                Sample = factor(Sample, levels = sample_order))

pp1 <- make_boxdot(combined, "frac_pEMT_with_inf", "Fraction of pEMT Spots Having Inf-Resp Neighbor (%)", percent_y = TRUE,
                   dot_palette = samp_cols, binwidth = 1/100, box_width = 0.6, dot_mm = 2, expend_mult_y = 0.01)
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C9_Core_vs_Edge_Samples_INF-Resp_Proximity_to_pEMT_vs_Non-pEMT_v3.pdf"), 
       plot = pp1, dpi = 300, width = 2, height = 3)


### ============ Panel B - Fraction of OSM expressing spots per metaprogram in a given sample ============
osm_pos <- merged_mat["OSM", ] %>% .[. > 1] %>% names()

per_mp_osm_expr_spots <- metadata %>% 
  dplyr::mutate(osm_pos = as.integer(ifelse(Key %in% osm_pos, 1, 0)),
                pEMT_Pattern = samples_metadata$pEMT_Pattern[match(Sample, samples_metadata$Sample)]) %>% 
  dplyr::filter(!is.na(MPid) & MPid != "Skeletal_Muscle") %>%      # remove MPs that are anatomical size specific
  dplyr::group_by(MPid) %>% 
  dplyr::summarize(Freq = mean(osm_pos)) %>% 
  dplyr::arrange(desc(Freq)) %>% 
  dplyr::mutate(MPid = factor(MPid, levels = rev(MPid)),
                pct = percent(Freq, accuracy = 0.1))

p <- ggplot(per_mp_osm_expr_spots, aes(x = Freq, y = MPid)) +
  geom_segment(aes(x = 0, xend = Freq, y = MPid, yend = MPid, colour = MPid),
               linewidth = 1) +
  geom_point(aes(fill = MPid), size = 5, colour = "black", shape = 21) +
  geom_text(aes(label = pct), nudge_x = 0.013, size = 5) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = state_cols) +
  scale_colour_manual(values = state_cols) +
  labs(x = "Fraction of Spots Expressing OSM per Metaprogram", y = NULL) +
  theme_classic() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16), legend.position = "none", aspect.ratio = 0.7)
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C9_Per-MP_Fraction_of_OSM_Expressing_Spots.pdf"), 
       plot = p, dpi = 300, width = 8, height = 5)


### ============ Panel C - Fraction of INF-Resp spots out of all spots in a given sample ============
inf_resp_fracs <- metadata %>% 
  dplyr::filter(!is.na(MPid)) %>% 
  dplyr::group_by(Sample) %>% 
  dplyr::summarise(frac_inflammatory = mean(MPid == "Inflammatory"),
                   n_total = n(),
                   n_inflammatory = sum(MPid == "Inflammatory")) %>%
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(frac_inflammatory)) %>% 
  dplyr::left_join(samples_metadata[, c("Sample", "pEMT_Pattern")], by = "Sample") %>% 
  dplyr::filter(pEMT_Pattern != "none") %>% 
  dplyr::mutate(pEMT_Pattern = factor(pEMT_Pattern, levels = c("Core_Enriched", "Margin_Enriched")))

pp3 <- make_boxdot(inf_resp_fracs, "frac_inflammatory", "Fraction of Inf-Resp Spots (%)", 
                   dot_palette = samp_cols, binwidth = 1/120, percent_y = TRUE, box_width = 0.6, dot_mm = 2, expend_mult_y = 0)
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C9_Core_vs_Edge_Samples_Fraction_of_Inf-Resp_Spots.pdf"), 
       plot = pp3, dpi = 300, width = 2, height = 3)


### ============ Panel D - Fraction of OSM-producing Inf-Resp spots out of all spots in a given sample ============
inf_meta <- metadata %>%
  dplyr::mutate(pEMT_Pattern = samples_metadata$pEMT_Pattern[match(Sample, samples_metadata$Sample)]) %>%
  dplyr::filter(!is.na(MPid) & MPid == "Inflammatory" & pEMT_Pattern != "none") %>%
  dplyr::select(Key, SpotID, Sample, MPid, pEMT_Pattern) %>%
  dplyr::mutate(osm_expr = as.numeric(merged_mat["OSM", Key]),
                osm_pos  = as.integer(osm_expr > 1),
                Sample = factor(Sample),
                pEMT_Pattern = factor(pEMT_Pattern, levels = c("Core_Enriched", "Margin_Enriched")))

# Logistic regression with sample as random effect - tests whether Inf-Resp spots in core samples have a higher probability of OSM expression than those in edge samples, adjusting for sample-to-sample variability.
fit <- lme4::glmer(osm_pos ~ pEMT_Pattern + (1 | Sample),
                   data = inf_meta,
                   family = binomial)
broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)

inf_meta <- inf_meta %>%
  mutate(pred_resp_conditional = predict(fit, type = "response"))
pred_sample <- inf_meta %>%
  group_by(Sample, pEMT_Pattern) %>%
  summarise(pred_prob = mean(pred_resp_conditional),
            n_spots   = n(),
            .groups = "drop")

p <- make_boxdot(pred_sample, "pred_prob", "Predicted Probability of OSM Expression in Inf-Resp spots", 
                 dot_palette = samp_cols, binwidth = 1/40, box_width = 0.6, dot_mm = 2, expend_mult_y = 0.005)
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C9_Core_vs_Edge_Samples_Predicted_OSM_Expressing_Inf-Resp_Spots.pdf"), 
       plot = p, dpi = 300, width = 2, height = 3)


### ============ Panel E - Inf-Resp Spots Spatial Coherence ============
coherence_scores <- readRDS(file = here("results/Generated_Data/all_samples_MPs_coherence_scores.rds"))
inf_resp_coher <- map_dbl(coherence_scores, "TNF_Signaling")
inf_resp_coher_df <- data.frame(Sample = names(inf_resp_coher),
                                Coherence = unname(inf_resp_coher),
                                pEMT_Pattern = samples_metadata$pEMT_Pattern[match(names(inf_resp_coher), samples_metadata$Sample)]) %>% 
  dplyr::filter(pEMT_Pattern != "none")

pp2 <- make_boxdot(inf_resp_coher_df, "Coherence", "Inf-Resp Spots Spatial Coherence Score", 
                   dot_palette = samp_cols, binwidth = 1/30, box_width = 0.6, dot_mm = 2, expend_mult_y = 0)
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C9_Core_vs_Edge_Samples_Spatial_Coherence_of_Inf-Resp_Spots_v3.pdf"), 
       plot = pp2, dpi = 300, width = 2, height = 3)


### ============ Panel F - Association between OSM expression and spatial coherence in Inf-Resp spots ============
osm_pos <- merged_mat["OSM", ] %>% .[. > 1] %>% names()
osm_pos_inf_spots <- metadata$Key[!is.na(metadata$MPid) & metadata$MPid == "Inflammatory" & metadata$Key %in% osm_pos]

test2 <- purrr::map(samples_metadata$Sample, function(samp) {
  samp_meta  <- filter(metadata, Sample == samp)
  adj_mat <- .neighbors_table_func(samp_meta, spot_class = "MPid", spot_name = "Key")
  inf_spots <- samp_meta %>%
    dplyr::filter(!is.na(MPid) & MPid == "Inflammatory") %>% 
    dplyr::pull(Key)
  inf_adj_mat <- adj_mat[inf_spots, , drop = FALSE]
  is_coherent <- apply(inf_adj_mat, 1, function(x) sum(is.na(x) | x == "Inflammatory") >= 4)
  
  list(coher_inf_spots = names(is_coherent[is_coherent]),
       not_coher_inf_spots = names(is_coherent[!is_coherent]))
})
names(test2) <- samples_metadata$Sample

long <- imap_dfr(test2, ~{
  tibble(
    Sample = .y,
    spot   = c(.x$coher_inf_spots, .x$not_coher_inf_spots),
    coherent = c(rep(TRUE,  length(.x$coher_inf_spots)),
                 rep(FALSE, length(.x$not_coher_inf_spots)))
  )
}) %>%
  mutate(OSM_pos = spot %in% osm_pos_inf_spots,
         pEMT_Pattern = samples_metadata$pEMT_Pattern[match(Sample, samples_metadata$Sample)])

# Build pooled 2x2 table
tab_pool <- long %>%
  mutate(coherent = coherent %in% c(TRUE, "TRUE", 1),
         OSM_pos = OSM_pos %in% c(TRUE, "TRUE", 1)) %>%
  count(coherent, OSM_pos) %>%
  pivot_wider(names_from = OSM_pos, values_from = n, values_fill = 0)

m_pool <- matrix(
  c(tab_pool$`TRUE` [tab_pool$coherent == TRUE],
    tab_pool$`FALSE`[tab_pool$coherent == TRUE],
    tab_pool$`TRUE` [tab_pool$coherent == FALSE],
    tab_pool$`FALSE`[tab_pool$coherent == FALSE]),
  nrow = 2, byrow = TRUE,
  dimnames = list(coherent = c("TRUE","FALSE"), OSM = c("TRUE","FALSE"))
)

# Fisher exact test
ft <- fisher.test(m_pool)
# Proportions + CIs - prop.test (Wilson/score CI)
k <- m_pool[, "TRUE"]
n <- rowSums(m_pool)
df <- tibble(coherent = c(TRUE, FALSE), k = k, n = n) %>%
  mutate(p = k / n)
cis <- t(mapply(function(ki, ni) {
  ci <- prop.test(ki, ni, correct = FALSE)$conf.int
  c(lo = ci[1], hi = ci[2])
}, df$k, df$n))

df$lo <- cis[, "lo"]
df$hi <- cis[, "hi"]
df$grp <- factor(df$coherent, levels = c(TRUE, FALSE),
                 labels = c("Coherent", "Non\ncoherent"))

# Plot
lw <- 0.9
bar_w <- 0.6
ann_text <- sprintf("OR = %.2f", unname(ft$estimate))

p <- ggplot(df, aes(x = grp, y = p, fill = grp)) +
  geom_col(width = bar_w, colour = "black", linewidth = lw, show.legend = FALSE) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = bar_w * 0.2, linewidth = lw) +
  scale_fill_manual(values = c("Coherent" = "#d94841", "Non-coherent" = "grey60")) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.00, 0.08))) +
  labs(x = NULL, y = "OSM Expression Fraction (%)") +
  ggpubr::theme_classic2(base_size = 10) +
  theme(
    axis.text  = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = lw),
    axis.line  = element_line(colour = "black", linewidth = lw),
    plot.margin = margin(4, 4, 2, 4, "pt")
  ) +
  annotate("text", x = 1.5, y = min(1, max(df$hi, na.rm = TRUE) * 1.08), label = ann_text, size = 4) +
  ggh4x::force_panelsizes(rows = unit(60, "mm"), cols = unit(30, "mm"), respect = TRUE)
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C9_OSM_Expression_in_Coherent_vs_Non-coherent_Inf-Resp_Spots.pdf"),
       plot = p, dpi = 300, width = 2, height = 3)


### ============ Panel G - OSM producing Inf-Resp colocalizing and/or neighboring pEMT spots ============
osm_pos_inf_spots <- metadata$Key[!is.na(metadata$MPid) & metadata$MPid == "Inflammatory" & metadata$Key %in% osm_pos]
cell_cols <- setdiff(names(dplyr::select(metadata, Fibroblast:Secretory_Norm)), "Inflammatory")
osm_spots_colocal <- metadata %>%
  dplyr::filter(Key %in% osm_pos_inf_spots) %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(across(all_of(cell_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  dplyr::mutate(pEMT_Pattern = samples_metadata$pEMT_Pattern[match(Sample, samples_metadata$Sample)]) %>%
  dplyr::filter(pEMT_Pattern != "none") %>%
  dplyr::mutate(total = rowSums(across(all_of(cell_cols))),
                across(all_of(cell_cols), ~ if_else(total > 0, .x / total, NA_real_))) %>%
  dplyr::select(Sample, all_of(cell_cols), pEMT_Pattern)

p <- make_boxdot(osm_spots_colocal, "pEMT", "Colocalization of OSM Expressing Inf-Resp and pEMT Spots", percent_y = TRUE, 
                 dot_palette = samp_cols, binwidth = 1/100, box_width = 0.6, dot_mm = 2, expend_mult_y = 0.01)
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C9_Core_vs_Edge_Samples_Colocalization_of_OSM_Expressing_Inf-Resp_and_pEMT_Spots.pdf"), 
       plot = p, dpi = 300, width = 2, height = 3)


### ============ Panel H - Fraction of pEMT spots out of all spots in a given sample ============
malignant_states <- c("Senescence","Secretory_Malig","Epithelial", "Cell_Cycle","Hypoxia", "pEMT", "LowQ")
pemt_fracs <- metadata %>% 
  dplyr::filter(!is.na(MPid) & MPid %in% malignant_states) %>% 
  dplyr::group_by(Sample) %>% 
  dplyr::summarise(frac_pEMT = mean(MPid == "pEMT"),
                   n_total = n(),
                   n_pEMT = sum(MPid == "pEMT")) %>%
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(frac_pEMT)) %>% 
  dplyr::left_join(samples_metadata[, c("Sample", "pEMT_Pattern")], by = "Sample") %>% 
  dplyr::filter(pEMT_Pattern != "none") %>% 
  dplyr::mutate(pEMT_Pattern = factor(pEMT_Pattern, levels = c("Core_Enriched", "Margin_Enriched")))

pp5 <- make_boxdot(pemt_fracs, "frac_pEMT", "Fraction of pEMT Spots (%)", binwidth = 1/100, percent_y = TRUE, 
                   dot_palette = samp_cols, box_width = 0.6, dot_mm = 2, expend_mult_y = 0.04)
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C9_Core_vs_Edge_Samples_Fraction_of_pEMT_Spots.pdf"), 
       plot = pp5, dpi = 300, width = 2, height = 3)


### ============ Panel I - Fraction of OSMR-producing pEMT spots out of all spots in a given sample ============
pEMT_meta <- metadata %>%
  dplyr::mutate(pEMT_Pattern = samples_metadata$pEMT_Pattern[match(Sample, samples_metadata$Sample)]) %>%
  dplyr::filter(!is.na(MPid) & MPid == "pEMT" & pEMT_Pattern != "none") %>%
  dplyr::select(Key, SpotID, Sample, MPid, pEMT_Pattern) %>%
  dplyr::mutate(osmr_expr = as.numeric(merged_mat["OSMR", Key]),
                osmr_pos  = as.integer(osmr_expr > 1),
                Sample = factor(Sample),
                pEMT_Pattern = factor(pEMT_Pattern, levels = c("Core_Enriched", "Margin_Enriched")))

# Logistic regression with sample as random effect - tests whether Inf-Resp spots in core samples have a higher probability of OSM expression than those in edge samples, adjusting for sample-to-sample variability.
fit <- lme4::glmer(osmr_pos ~ pEMT_Pattern + (1 | Sample),
                   data = pEMT_meta,
                   family = binomial)
broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)

pEMT_meta <- pEMT_meta %>%
  mutate(pred_resp_conditional = predict(fit, type = "response"))
pred_sample <- pEMT_meta %>%
  group_by(Sample, pEMT_Pattern) %>%
  summarise(pred_prob = mean(pred_resp_conditional),
            n_spots   = n(),
            .groups = "drop")
p <- make_boxdot(pred_sample, "pred_prob", "Predicted Probability of OSMR Expression in pEMT spots", 
                 dot_palette = samp_cols, binwidth = 1/30, box_width = 0.6, dot_mm = 2, expend_mult_y = 0.005)
ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C9_Core_vs_Edge_Samples_Predicted_OSMR_Expressing_pEMT_Spots.pdf"), 
       plot = p, dpi = 300, width = 2, height = 3)


# HPV On/Off/Neg State Fractions (R1C6) ----------------------------------------

# Load per-sample parameter table
sample_parameters <- readRDS(file = here("metadata/per_sample_parameters_df.rds"))
# Load sample level metadata
samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))
# Load merged metadata tibble
metadata <- readRDS(file = here("metadata/merged_metadata_Extend.rds"))
# Load state colors 
state_cols <- read.delim(file = here("aux_data/state_col_tab_extend.tsv"), sep = "\t", header = TRUE)
state_cols <- setNames(state_cols$V2, state_cols$V1)


## 1) Label HPVon / HPVoff per HPV-positive sample
# helper: return tibble(Key, SpotID, Sample, HPV_spot) for one sample
hpv_spot_calls <- function(samp_id, site_id, metadata, base_dir = here::here()) {
  data_dir <- paste(here("data"), site_id, samp_id, sep = "/")
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

  samp_meta <- metadata %>% filter(Sample == samp_id) %>% dplyr::select(Key, SpotID, Sample)
  
  # keep only columns present in metadata (order matters)
  m <- m[, samp_meta$SpotID, drop = FALSE]
  
  hpv_rows <- grep("^hpv", rownames(m), ignore.case = TRUE)
  hpv_on   <- if (length(hpv_rows) == 0L) {
    rep(FALSE, ncol(m))
  } else {
    colSums(m[hpv_rows, , drop = FALSE]) > 0   # >0 total HPV counts
  }
  
  tibble(
    Key      = paste(samp_id, colnames(m), sep = "_"),
    SpotID   = colnames(m),
    Sample   = samp_id,
    HPV_spot = ifelse(hpv_on, "HPVon", "HPVoff")
  )
}

# build calls for all HPVpos samples in your table
hpvp <- samples_metadata %>%
  dplyr::select(Sample, Site, HPV) %>%
  dplyr::filter(HPV == "HPVpos")

hpv_calls <- purrr::map2_dfr(hpvp$Sample, "Oropharynx_CA",
                             ~ hpv_spot_calls(.x, .y, metadata))


## 2) Per-spot group (HPVneg / HPVoff / HPVon) + MPid
malignant_states <- c("Senescence", "Cell_Cycle", "Hypoxia", "pEMT", "Epithelial", "Secretory_Malig")
# join with sample-level HPV status and MPid per spot
spot_groups <- metadata %>%
  dplyr::filter(Sample %in% samples_metadata$Sample[samples_metadata$Site == "Oropharynx"]) %>% 
  dplyr::select(Key, SpotID, Sample, MPid) %>%
  dplyr::left_join(samples_metadata %>% dplyr::select(Sample, HPV), by = "Sample") %>%
  # add HPVon/off where available (only HPVpos samples get a value)
  dplyr::left_join(hpv_calls %>% dplyr::select(Key, HPV_spot), by = "Key") %>%
  dplyr::mutate(
    HPV_group = dplyr::case_when(
      HPV == "HPVneg"              ~ "HPVneg",
      HPV == "HPVpos" & HPV_spot=="HPVon"  ~ "HPVon",
      HPV == "HPVpos" & HPV_spot=="HPVoff" ~ "HPVoff",
      TRUE ~ NA_character_
    ),
    HPV_group = factor(HPV_group, levels = c("HPVneg","HPVoff","HPVon"))
  ) %>%
  dplyr::filter(!is.na(HPV_group), !is.na(MPid), MPid %in% malignant_states) %>% 
  droplevels()


## 3) Fractions of cellular states per group
cell_fracs <- spot_groups %>%
  count(MPid, HPV_group, name = "n") %>%
  complete(MPid, HPV_group, fill = list(n = 0)) %>% 
  group_by(MPid) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

order_levels <- cell_fracs %>%
  dplyr::filter(HPV_group == "HPVneg") %>%
  dplyr::arrange(dplyr::desc(frac)) %>%
  dplyr::pull(MPid)

cell_fracs$MPid <- factor(cell_fracs$MPid, levels = rev(order_levels))

overall <- spot_groups %>%
  count(HPV_group, name = "n") %>%
  mutate(MPid = "All cells") %>%
  group_by(MPid) %>% mutate(frac = n / sum(n)) %>% ungroup()

cell_fracs_all <- bind_rows(
  overall %>% dplyr::select(HPV_group, MPid, n, frac),
  cell_fracs
)


## 4) Plot (stacked bars, fractions within each MPid)
pal <- c(HPVneg = "#9BD489", HPVoff = "#0F6A2B", HPVon = "#E3362D")  # tweak to taste

p <- ggplot(cell_fracs, aes(x = MPid, y = frac, fill = HPV_group)) +
  geom_col(width = 0.85, color = "white", size = 0.3, position = "stack") +
  scale_fill_manual(values = pal, name = NULL,
                    breaks = c("HPVneg","HPVoff","HPVon"),
                    labels = c("HPVneg cells","HPVoff cells","HPVon cells")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Fraction of cells") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p1 <- ggplot(overall, aes(x = MPid, y = frac, fill = HPV_group)) +
  geom_col(width = 0.85, color = "white", size = 0.3) +
  scale_fill_manual(values = pal, name = NULL,
                    labels = c("HPVneg cells","HPVoff cells","HPVon cells")) +
  scale_y_continuous(labels = scales::percent_format(1), limits = c(0,1), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Fraction of cells") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot2::ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C6_HPV_on_off_neg_MP_Fractions.pdf"), plot = p, width = 6, height = 4, dpi = 300, device = "pdf")
ggplot2::ggsave(filename = here("results/Paper_Figures/Revision_Figs/R1_C6_OPSCC_HPV_Groups_Spot_Fractions.pdf"), plot = p1, width = 3, height = 4, dpi = 300, device = "pdf")
