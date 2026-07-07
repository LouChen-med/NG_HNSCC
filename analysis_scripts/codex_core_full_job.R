#!/usr/bin/env Rscript

## ---------- robust arg parser (keeps everything after first '=') ----------
parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  if (!length(a)) return(list())
  # allow both "--key=val" and "key=val"
  a <- sub("^--", "", a)
  keys <- sub("=.*$", "", a)
  vals <- sub("^[^=]*=", "", a)
  as.list(stats::setNames(vals, keys))
}
`%||%` <- function(a,b) if (is.null(a) || !nzchar(a)) b else a

opt <- parse_args()

## ---------- required args ----------
req <- c("core_id","cells_csv","mask_dir","bound_dir","tools_r",
         "out_dir_res","out_dir_out","out_dir_plot","px2um")
miss <- setdiff(req, names(opt))
if (length(miss)) stop("Missing args: ", paste(miss, collapse=", "))

core_id     <- opt$core_id
cells_csv   <- opt$cells_csv
mask_dir    <- opt$mask_dir
bound_dir   <- opt$bound_dir
tools_r     <- opt$tools_r
out_dir_res <- opt$out_dir_res
out_dir_out <- opt$out_dir_out
out_dir_plot<- opt$out_dir_plot
flip_map_csv<- opt$flip_map_csv %||% ""
px2um       <- as.numeric(opt$px2um)
ncores_dist <- as.integer(opt$ncores %||% "4")

# distance options
source_label <- opt$source_label %||% "pEMT"
boundary_flag<- tolower(opt$boundary %||% "true") %in% c("true","t","1","yes","y")
force_rebuild<- tolower(opt$force_rebuild %||% "false") %in% c("true","t","1","yes","y")

## ---------- target string parser ----------
# Accepts "Name=A|B|C;Other=X"  (semicolon-separated groups; values pipe-separated)
parse_targets <- function(x) {
  if (is.null(x) || !nzchar(x)) return(list())
  # normalize group separators to ';'
  parts <- strsplit(x, ";", fixed = TRUE)[[1]]
  out <- list()
  for (p in parts) {
    p <- trimws(p)
    if (!nzchar(p)) next
    kv <- strsplit(p, "=", fixed = TRUE)[[1]]
    if (length(kv) < 2) next
    name <- trimws(kv[1])
    # keep everything after the FIRST '=' to allow '=' in values (rare but safe)
    vals <- sub("^[^=]*=", "", p)
    toks <- strsplit(vals, "|", fixed = TRUE)[[1]]
    toks <- unique(trimws(toks[nzchar(toks)]))
    out[[name]] <- toks
  }
  out
}

targets_raw <- parse_targets(opt$targets %||% "")

suppressPackageStartupMessages({
  library(sf); library(readr); library(dplyr)
})
sf::sf_use_s2(FALSE)

## ---------- idempotent skip for *both* outputs ----------
res_rds <- file.path(out_dir_res, sprintf("Core_%s.rds", core_id))
out_rds <- file.path(out_dir_out, sprintf("Core_%s.rds", core_id))
if (file.exists(res_rds) && file.exists(out_rds) && !force_rebuild) {
  message("Both outputs exist for ", core_id, " — skipping.")
  quit(save="no", status=0)
}

## ---------- tools ----------
if (!file.exists(tools_r)) stop("tools_r not found: ", tools_r)
source(tools_r)

## ---------- cache: load res if present unless force_rebuild ----------
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)
res <- NULL
if (!force_rebuild && file.exists(res_rds)) {
  message("Using cached processed core for ", core_id, " → ", res_rds)
  res <- tryCatch(readRDS(res_rds), error = function(e) {message("Failed to read cached res: ", e$message); NULL})
}

## ---------- build res if cache missing ----------
if (is.null(res)) {
  if (!file.exists(cells_csv)) stop("cells_csv not found: ", cells_csv)
  mask_path  <- file.path(mask_dir,  sprintf("Core_%s_mask.geojson",  core_id))
  bound_path <- file.path(bound_dir, sprintf("Core_%s_tumor.geojson", core_id))
  if (!file.exists(mask_path))  stop("Mask missing: ",  mask_path)
  if (!file.exists(bound_path)) stop("Boundary missing: ", bound_path)

  cells    <- readr::read_csv(cells_csv, show_col_types = FALSE)
  cells_sf <- sf::read_sf(mask_path)

  flip_map <- NULL
  if (nzchar(flip_map_csv)) {
    if (!file.exists(flip_map_csv)) stop("flip_map_csv not found: ", flip_map_csv)
    fm <- readr::read_csv(flip_map_csv, show_col_types = FALSE)
    if (!all(c("core","flip") %in% names(fm))) stop("flip_map_csv must have columns: core, flip")
    flip_map <- stats::setNames(as.logical(fm$flip), fm$core)
  }

  res <- codex_process_core_viz(
    cells_sf_core    = cells_sf,
    cells_tbl_core   = cells,
    boundary_geojson = bound_path,
    core_id          = core_id,
    flip_map         = flip_map,
    px2um            = px2um
  )

  # light meta for sanity
  res$meta <- list(core_id = core_id, px2um = px2um, built_at = Sys.time())
  saveRDS(res, res_rds)
  message("Saved processed core → ", res_rds)
}

stopifnot(is.list(res), !is.null(res$prep), !is.null(res$nest_bound))

## ---------- log labels & targets (raw vs filtered-to-present) ----------
allowed_labels <- sort(unique(res$prep$cells_sf$assignment))

cat("Available labels in core:", paste(allowed_labels, collapse = ", "), "\n")
cat("Parsed targets (raw):\n")
if (length(targets_raw)) {
  for (nm in names(targets_raw)) {
    cat(sprintf("  %s={%s}\n", nm, paste(targets_raw[[nm]], collapse="|")))
  }
} else {
  cat("  (none)\n")
}

# filter to labels actually present in this core
targets_list <- lapply(targets_raw, function(v) intersect(v, allowed_labels))
targets_list <- targets_list[lengths(targets_list) > 0]

cat("Parsed targets (filtered to present labels):\n")
if (length(targets_list)) {
  for (nm in names(targets_list)) {
    cat(sprintf("  %s={%s}\n", nm, paste(targets_list[[nm]], collapse="|")))
  }
} else {
  cat("  (none)\n")
}
cat("Boundary flag:", boundary_flag, "\n")

## ---------- 1) plot and save core with boundary ----------
my_pal <- c("Fibroblast" = "#88CCEE", "Endothelial" = "#DDCC77", "B-cell" = "#AA4499", "Tctx-cell" = "sandybrown", 
            "Neutrophil" = "#DD1C77", "M2-like_Macrophage" = "#FF2400", "M1-like_Macrophage" = "pink2", "Monocyte" = "#1b9e77",
            "Myeloid" = "#ADDD8E", "Epithelial" = "#999933", "pEMT" = "#E4CAB5", "Senescence" = "#CC6677", "Cell-Cycle" = "#117733", "Hypoxia" = "#888888", "Undefined" = "grey85")
p <- codex_plot_core(res$prep$cells_sf,
  cell_palette = my_pal,
  boundaries = res$nest_bound,
  boundary_linewidths = 0.8,
  boundary_colours = "white")
dir.create(out_dir_plot, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(out_dir_plot, sprintf("Core_%s_with_Boundary.png", core_id)), p, width = 14, height = 12, dpi = 300)

## ---------- 2) infiltration (myeloid example) ----------
myeloid_levels <- c("Myeloid","Neutrophil","M2-like_Macrophage")
myel_infil     <- codex_run_infiltration(res, myeloid_levels)

## ---------- 3) pEMT minimal distances ----------
dist_tab <- codex_min_dist_summary(
  prep_obj     = res,
  source       = source_label,
  targets_list = targets_list,   # may be empty
  boundary     = boundary_flag,  # compute boundary distances if TRUE
  px2um        = px2um,
  n_cores      = ncores_dist
)

## ---------- 4) save analysis ----------
dir.create(out_dir_out, recursive = TRUE, showWarnings = FALSE)
out <- list(
  core_id        = core_id,
  Myeloid_infilt = myel_infil[c("by_nest","summary")],
  Distance_table = dist_tab
)
saveRDS(out, out_rds)
message("Saved: ", out_rds)