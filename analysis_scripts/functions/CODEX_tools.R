# codex_tumor_tools.R
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Toolkit for CODEX TMA spatial analysis:
#
#   1) Prepare core + align tumour boundary (external or de-novo)
#   2) Infiltration metrics
#   3) Minimal distance summaries (cell→cell, cell→boundary) incl. multi-target
#   4) Pairwise proximity (fixed-radius permutations)
#   5) Tumour zones (edge / mid / core) from boundary offsets
#   6) Density (global + local-in-mask)
#   7) Plotting
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# IMPORTANT ASSUMPTIONS ABOUT INPUT DATA
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. cells_sf_core  : sf MULTIPOLYGON layer of cell segmentation for ONE core.
#                     Must contain columns:
#                       id          – unique cell id (character)
#                       assignment  – cell type label (factor or chr)
#                     Geometry coordinates are all in the SAME planar system
#                     (pixels or µm) and already oriented correctly (flip_y done).
#                     CRS may be NA; that's fine (preferred) for image coords.
# 2. cells_tbl_core : plain tibble/data.frame with at least:
#                       cell_id     – same IDs as cells_sf_core$id
#                       assignment  – cell type (redundant but convenient)
#                       centroid_x / centroid_y (OPTIONAL; if missing we derive
#                         centroids from polygons)
#                     You may include intensity, markers, etc.; not required.
#
# If you do **NOT** have a separate metadata table, you can pass `cells_tbl_core = NULL`
# and everything will be derived from `cells_sf_core`.


# Load Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(dbscan)
  library(parallel)
  library(rlang)
  library(tibble)
  library(RColorBrewer)
  library(viridis)
})


# Global Defaults (edit to suit the project) -------------------------------
codex_defaults <- list(
  malignant_levels = c("Epithelial","pEMT"),
  myeloid_levels   = c("Myeloid","M1-like_Macrophage","M2-like_Macrophage","Monocyte","Neutrophil","Cycling_Myeloid"),
  eps      = 30,      # DBSCAN eps (in your coordinate units)
  minPts   = 10,      # DBSCAN minPts
  size_min = 50,      # drop malignant clusters smaller than this
  simplify_tol = 1,   # pre-union simplify tolerance
  r_merge      = 10,  # soft merge radius (gap-closing)
  hole_min_cells       = 200,  # carve interior pockets with >=200 cells
  hole_max_malign_frac = 0.20  # carve pockets with <=20% malignant
)


# Utility Helpers ------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a
.codex_use_planar <- function() { if (sf::sf_use_s2()) sf::sf_use_s2(FALSE); invisible(TRUE) }
codex_parse_core_id <- function(x) sub("^fov0_([A-Za-z0-9]+)_.*$", "\\1", x)

# Flip Y wrt reference bbox (used for mask alignment)
flip_y_to_ref <- function(sf_obj, ref_sf) {
  bb   <- sf::st_bbox(ref_sf)
  ymin <- as.numeric(bb["ymin"]); ymax <- as.numeric(bb["ymax"])
  flip_coords <- function(m) { m[,2] <- ymin + ymax - m[,2]; m }
  g <- sf::st_geometry(sf_obj)
  g_flip <- lapply(g, function(geom) {
    cls <- class(geom)[2]
    switch(cls,
           "POINT"           = sf::st_point(flip_coords(matrix(sf::st_coordinates(geom), ncol = 2))[1,], dim = "XY"),
           "MULTIPOINT"      = sf::st_multipoint(flip_coords(sf::st_coordinates(geom)[,1:2, drop = FALSE])),
           "LINESTRING"      = sf::st_linestring(flip_coords(sf::st_coordinates(geom)[,1:2, drop = FALSE])),
           "MULTILINESTRING" = { parts <- lapply(geom, function(ls) sf::st_linestring(flip_coords(ls))); sf::st_multilinestring(parts) },
           "POLYGON"         = { rings <- lapply(geom, flip_coords); sf::st_polygon(rings) },
           "MULTIPOLYGON"    = { polys <- lapply(geom, function(pg) lapply(pg, flip_coords)); sf::st_multipolygon(polys) },
           stop("Unsupported geometry type: ", cls)
    )
  })
  sf::st_geometry(sf_obj) <- sf::st_sfc(g_flip, crs = sf::st_crs(sf_obj))
  sf_obj
}


# Prepare Core & Align Tumour Boundary ---------------------------------------

# Prepare: enforce planar, add assignment, create centroids; returns list
codex_prepare_core <- function(cells_sf_core,
                               cells_tbl_core = NULL,
                               flip_y_axis = TRUE,
                               na.rm = TRUE,
                               set_crs = FALSE,
                               px2um = NULL,
                               malignant_levels = codex_defaults$malignant_levels) {
  .codex_use_planar()
  stopifnot(inherits(cells_sf_core, "sf"))
  if (!"id" %in% names(cells_sf_core))
    stop("cells_sf_core must contain column 'id' (unique cell ids).")
  if (!"assignment" %in% names(cells_sf_core)) {
    if (is.null(cells_tbl_core))
      stop("No 'assignment' in cells_sf_core and no cells_tbl_core provided.")
    cells_sf_core <- dplyr::left_join(
      cells_sf_core,
      dplyr::select(cells_tbl_core, id = cell_id, assignment),
      by = "id"
    )
  }
  # flip & planar
  cells_sf_core <- if (flip_y_axis) flip_y_to_ref(cells_sf_core, cells_sf_core) else cells_sf_core
  cells_sf_core <- sf::st_set_crs(cells_sf_core, NA)
  if (set_crs && !is.null(px2um)) {
    sf::st_geometry(cells_sf_core) <- sf::st_geometry(cells_sf_core) * px2um
  }
  if (na.rm) cells_sf_core <- cells_sf_core %>% dplyr::filter(!is.na(assignment))

  cent <- sf::st_centroid(cells_sf_core) %>% dplyr::select(cell_id = id, assignment)
  coords <- sf::st_coordinates(cent)
  cells_tbl_out <- tibble::tibble(
    cell_id    = cent$cell_id,
    assignment = cent$assignment,
    centroid_x = coords[,1],
    centroid_y = coords[,2]
  )
  cent_mal <- cent %>% dplyr::filter(.data$assignment %in% malignant_levels)
  list(cells_sf = cells_sf_core, cells_tbl = cells_tbl_out, cent = cent, cent_mal = cent_mal)
}

# External boundary loader (geojson) → MULTIPOLYGON + MULTILINESTRING
codex_load_external_boundary <- function(path_geojson, ref_sf, id_col = "nest_id") {
  pol <- sf::read_sf(path_geojson) %>%
    sf::st_cast("MULTIPOLYGON", warn = FALSE) %>%
    sf::st_make_valid() %>%
    sf::st_set_crs(sf::st_crs(ref_sf))
  if (!id_col %in% names(pol)) pol[[id_col]] <- seq_len(nrow(pol))
  list(poly = pol, bound = sf::st_cast(pol, "MULTILINESTRING"))
}

# Quick mask score: fraction of malignant centroids inside poly
.codex_score_mask_fast <- function(poly, cells_sf, mal_classes = c("Epithelial","pEMT")) {
  mal <- cells_sf[cells_sf$assignment %in% mal_classes, ]
  if (!nrow(mal)) return(0)
  poly_u <- sf::st_union(sf::st_geometry(poly))
  pts    <- sf::st_centroid(mal)
  inside <- lengths(sf::st_within(pts, poly_u)) > 0
  mean(inside)
}

# Auto align mask (flip if needed to match cell coordinate frame)
codex_align_boundary <- function(poly, cells_sf,
                                 flip_fun = flip_y_to_ref,
                                 mal_classes = c("Epithelial","pEMT"),
                                 force_flip_mask = NA) {
  poly <- sf::st_set_crs(poly, sf::st_crs(cells_sf))
  if (isTRUE(force_flip_mask)) {
    poly1 <- flip_fun(poly, cells_sf)
    return(list(poly = poly1, bound = sf::st_cast(poly1, "MULTILINESTRING"), flipped = TRUE,  score = NA))
  }
  if (identical(force_flip_mask, FALSE)) {
    return(list(poly = poly,  bound = sf::st_cast(poly,  "MULTILINESTRING"), flipped = FALSE, score = NA))
  }
  s0 <- .codex_score_mask_fast(poly, cells_sf, mal_classes)
  poly1 <- flip_fun(poly, cells_sf)
  s1 <- .codex_score_mask_fast(poly1, cells_sf, mal_classes)
  if (s1 > s0) list(poly = poly1, bound = sf::st_cast(poly1,"MULTILINESTRING"), flipped = TRUE,  score = s1)
  else         list(poly = poly,  bound = sf::st_cast(poly, "MULTILINESTRING"), flipped = FALSE, score = s0)
}

# DBSCAN cluster malignant centroids
codex_cluster_malignant <- function(cent_mal,
                                    eps      = codex_defaults$eps,
                                    minPts   = codex_defaults$minPts,
                                    size_min = codex_defaults$size_min) {
  .codex_use_planar()
  mat <- sf::st_coordinates(cent_mal)
  db  <- dbscan::dbscan(mat, eps = eps, minPts = minPts)
  cent_mal$nest_id <- db$cluster
  cent_mal_lab <- cent_mal %>%
    dplyr::filter(.data$nest_id > 0) %>%
    dplyr::add_count(.data$nest_id, name = "cluster_size") %>%
    dplyr::filter(.data$cluster_size >= size_min) %>%
    dplyr::select(-.data$cluster_size)
  cent_mal_lab
}

# Nest boundary from malignant polygons: (tight union) → soften → valid
codex_nest_union_tight <- function(cells_sf_core, cent_mal_lab,
                                   simplify_tol = codex_defaults$simplify_tol,
                                   r_gap_close  = 0) {
  .codex_use_planar()
  mal_cells <- cells_sf_core %>%
    dplyr::left_join(cent_mal_lab %>% sf::st_drop_geometry() %>% dplyr::select(id = cell_id, nest_id), by = "id") %>%
    dplyr::filter(!is.na(.data$nest_id))
  if (simplify_tol > 0)
    mal_cells <- sf::st_simplify(mal_cells, dTolerance = simplify_tol, preserveTopology = TRUE)
  nest_poly <- mal_cells %>%
    dplyr::group_by(.data$nest_id) %>%
    dplyr::summarise(do_union = TRUE) %>%
    sf::st_make_valid()
  if (r_gap_close > 0) {
    nest_poly <- nest_poly %>%
      dplyr::mutate(geometry = sf::st_buffer(.data$geometry,  r_gap_close)) %>%
      dplyr::mutate(geometry = sf::st_buffer(.data$geometry, -r_gap_close)) %>%
      sf::st_make_valid()
  }
  nest_poly
}

codex_nest_soften <- function(nest_poly, r_merge = codex_defaults$r_merge) {
  .codex_use_planar()
  if (r_merge <= 0) return(nest_poly)
  nest_poly %>%
    dplyr::mutate(geometry = sf::st_buffer(.data$geometry,  r_merge)) %>%
    dplyr::mutate(geometry = sf::st_buffer(.data$geometry, -r_merge)) %>%
    sf::st_make_valid()
}

# Optional: carve interior low-malignant “lakes”
codex_nest_carve_tme <- function(cells_sf_core, nest_soft,
                                 malignant_levels      = codex_defaults$malignant_levels,
                                 hole_min_cells        = codex_defaults$hole_min_cells,
                                 hole_max_malign_frac  = codex_defaults$hole_max_malign_frac,
                                 r_cell = NULL) {
  .codex_use_planar()
  if (is.null(r_cell)) {
    med_area <- stats::median(as.numeric(sf::st_area(cells_sf_core)))
    r_cell   <- sqrt(med_area / pi)
  }
  cells_pts <- sf::st_centroid(cells_sf_core)

  refine_one <- function(nid) {
    outer  <- nest_soft %>% dplyr::filter(.data$nest_id == nid)
    inside <- sf::st_within(cells_pts, outer, sparse = FALSE)[,1]
    pts_in <- cells_pts[inside,]
    pts_in$is_malign <- pts_in$assignment %in% malignant_levels
    eps_sub <- r_cell * 2
    db <- dbscan::dbscan(sf::st_coordinates(pts_in), eps = eps_sub, minPts = 5)
    pts_in$sub_id <- db$cluster
    sub_sum <- pts_in %>% sf::st_drop_geometry() %>%
      dplyr::group_by(.data$sub_id) %>%
      dplyr::summarise(n_tot = dplyr::n(), n_mal = sum(.data$is_malign), frac_mal = .data$n_mal/.data$n_tot, .groups = "drop")
    cut_ids <- sub_sum %>%
      dplyr::filter(.data$sub_id > 0, .data$n_tot >= hole_min_cells, .data$frac_mal <= hole_max_malign_frac) %>%
      dplyr::pull(.data$sub_id)
    if (!length(cut_ids)) return(outer)
    hole_geom <- pts_in %>% dplyr::filter(.data$sub_id %in% cut_ids) %>%
      sf::st_buffer(r_cell) %>%
      dplyr::summarise(do_union = TRUE) %>% sf::st_make_valid()
    touch <- sf::st_touches(hole_geom, outer, sparse = FALSE)[,1]
    if (any(touch)) hole_geom <- hole_geom[!touch,]
    if (!nrow(hole_geom)) return(outer)
    sf::st_sf(nest_id = nid, geometry = sf::st_difference(sf::st_geometry(outer), sf::st_union(hole_geom)))
  }
  purrr::map_dfr(nest_soft$nest_id, refine_one) %>% sf::st_make_valid()
}

codex_nest_boundary <- function(nest_poly) sf::st_cast(nest_poly, "MULTILINESTRING")

# Wrapper: process core for downstream analysis (no infiltration here)
codex_process_core_viz <- function(cells_sf_core,
                                   cells_tbl_core   = NULL,
                                   core_id          = NULL,
                                   flip_map         = NULL,
                                   flip_y_axis      = TRUE,          # fallback if core not in flip_map
                                   boundary_geojson = NULL,          # optional external boundary
                                   # explicit args for codex_prepare_core
                                   malignant_levels = codex_defaults$malignant_levels,
                                   na.rm            = TRUE,
                                   set_crs          = FALSE,
                                   px2um            = NULL,
                                   # de-novo boundary params
                                   eps              = codex_defaults$eps,
                                   minPts           = codex_defaults$minPts,
                                   size_min         = codex_defaults$size_min,
                                   simplify_tol     = codex_defaults$simplify_tol,
                                   r_gap_close      = codex_defaults$r_merge/2,
                                   r_merge          = codex_defaults$r_merge,
                                   carve_tme        = FALSE) {
  # Decide flip per core (flip_map overrides flip_y_axis when available)
  flip_cells <- if (!is.null(core_id) && !is.null(flip_map) && core_id %in% names(flip_map)) {
    isTRUE(flip_map[[core_id]])
  } else {
    isTRUE(flip_y_axis)
  }

  # 1) Prepare core (only pass args codex_prepare_core actually takes)
  prep <- codex_prepare_core(
    cells_sf_core      = cells_sf_core,
    cells_tbl_core     = cells_tbl_core,
    malignant_levels   = malignant_levels,
    flip_y_axis        = flip_cells,
    na.rm              = na.rm,
    set_crs            = set_crs,
    px2um              = px2um
  )

  # 2) Boundary (external vs de-novo)
  if (!is.null(boundary_geojson)) {
    ext <- codex_load_external_boundary(boundary_geojson, ref_sf = prep$cells_sf)
    alg <- codex_align_boundary(ext$poly, prep$cells_sf)  # auto-choose flip for mask if needed
    nest_poly  <- alg$poly
    nest_bound <- alg$bound
    aligned    <- list(flipped = alg$flipped, score = alg$score)
  } else {
    cent_mal_lab <- codex_cluster_malignant(prep$cent_mal, eps = eps, minPts = minPts, size_min = size_min)
    nest_poly <- codex_nest_union_tight(prep$cells_sf, cent_mal_lab, simplify_tol = simplify_tol, r_gap_close = r_gap_close)
    nest_poly <- codex_nest_soften(nest_poly, r_merge = r_merge)
    if (isTRUE(carve_tme))
      nest_poly <- codex_nest_carve_tme(prep$cells_sf, nest_poly)
    nest_bound <- codex_nest_boundary(nest_poly)
    aligned <- list(flipped = FALSE, score = NA_real_)
  }

  list(core_id = core_id, flip_used = flip_cells, prep = prep, nest_poly = nest_poly, nest_bound = nest_bound, aligned = aligned)
}


# Infiltration Computation ---------------------------------------------------

codex_infiltration <- function(cells_sf_core, nest_poly, target_types) {
  .codex_use_planar()
  infil_pts <- sf::st_centroid(cells_sf_core) %>% dplyr::filter(.data$assignment %in% target_types)
  infil_in  <- sf::st_join(infil_pts, nest_poly, join = sf::st_within, left = FALSE)
  by_nest   <- infil_in %>% sf::st_drop_geometry() %>% dplyr::count(.data$nest_id, name = "n_infil")
  tumour_area <- sum(sf::st_area(nest_poly))
  summary <- tibble::tibble(
    total_infil = nrow(infil_in),
    tumour_area_um2 = as.numeric(tumour_area),
    density = nrow(infil_in) / as.numeric(tumour_area)
  )
  list(by_cell = infil_in, by_nest = by_nest, summary = summary)
}

# Wrapper
codex_run_infiltration <- function(prep_obj, target_types) {
  codex_infiltration(prep_obj$prep$cells_sf, prep_obj$nest_poly, target_types)
}


# Minimal Distances -----------------------------------------------------------

# Distances from a source cell-type to (a) one/many target groups and/or (b) boundary
codex_min_dist_summary <- function(prep_obj,
                                   source,                      # e.g. "pEMT" or c("pEMT","Epithelial")
                                   targets_list = NULL,         # char vec or named list of char vecs; ignored if boundary-only
                                   boundary = FALSE,            # include source→boundary distances?
                                   px2um = NULL,                # µm per pixel; else set pix_per_um
                                   pix_per_um = NULL,           # pixels per µm (alt to px2um)
                                   n_cores = 1,                 # parallel across target groups
                                   return_per_cell = FALSE) {   # also return per-cell distances
  # ---- checks / units --------------------------------------------------------
  stopifnot(is.list(prep_obj),
            !is.null(prep_obj$prep$cells_sf),
            inherits(prep_obj$prep$cells_sf, "sf"),
            "assignment" %in% names(prep_obj$prep$cells_sf))

  # derive px→µm
  if (is.null(px2um)) {
    if (!is.null(prep_obj$prep$px2um)) {
      px2um <- prep_obj$prep$px2um
    } else if (!is.null(pix_per_um)) {
      px2um <- 1 / pix_per_um
    } else {
      stop("Provide `px2um` or `pix_per_um`, or store scale as `prep_obj$prep$px2um`.")
    }
  }

  cells_sf <- prep_obj$prep$cells_sf
  bnd      <- prep_obj$nest_bound

  # Ensure point geometry & a usable ID
  pts <- if (inherits(sf::st_geometry(cells_sf), "sfc_POINT")) cells_sf else sf::st_centroid(cells_sf)
  if (!"cell_id" %in% names(pts)) {
    pts$cell_id <- if ("id" %in% names(pts)) pts$id else seq_len(nrow(pts))
  }

  # Source subset
  src <- dplyr::filter(pts, .data$assignment %in% !!source)
  if (!nrow(src)) {
    out <- tibble::tibble(
      target    = character(0),
      n_source  = 0L,
      mean_um   = numeric(0),
      sd_um     = numeric(0),
      median_um = numeric(0),
      p25_um    = numeric(0),
      p75_um    = numeric(0)
    )
    return(if (return_per_cell) list(summary = out,
                                     per_cell = tibble::tibble(cell_id = character(),
                                                               dist_um = numeric(),
                                                               target = character()))
           else out)
  }

  # ---- helper to summarize a numeric vector ----------------------------------
  .summ <- function(v_um) tibble::tibble(
    n_source  = length(v_um),
    mean_um   = mean(v_um, na.rm = TRUE),
    sd_um     = stats::sd(v_um, na.rm = TRUE),
    median_um = stats::median(v_um, na.rm = TRUE),
    p25_um    = stats::quantile(v_um, 0.25, na.rm = TRUE),
    p75_um    = stats::quantile(v_um, 0.75, na.rm = TRUE)
  )

  # ---- boundary-only mode ----------------------------------------------------
  do_boundary <- function() {
    if (is.null(bnd) || !inherits(bnd, "sf") || !nrow(bnd))
      stop("`boundary=TRUE` but `prep_obj$nest_bound` is missing or empty.")
    d_px <- apply(as.matrix(sf::st_distance(src, bnd)), 1, min)
    d_um <- as.numeric(d_px) * px2um
    s <- dplyr::bind_cols(tibble::tibble(target = "Boundary"), .summ(d_um))
    if (!return_per_cell) return(list(summary = s, per_cell = NULL))
    pc <- tibble::tibble(cell_id = src$cell_id, dist_um = d_um, target = "Boundary")
    list(summary = s, per_cell = pc)
  }

  # ---- normalize targets_list ------------------------------------------------
  # cases:
  #   * NULL (OK only if boundary=TRUE)
  #   * character vector -> one row per label (own group)
  #   * named list -> each element is a (possibly multi-label) group
  if (is.null(targets_list)) {
    if (!boundary) stop("Provide `targets_list` or set `boundary=TRUE`.")
    # boundary-only
    return(do_boundary())
  }

  if (is.character(targets_list)) {
    # treat each label as its own named group
    nm <- unique(targets_list)
    targets_list <- stats::setNames(as.list(nm), nm)
  } else if (is.list(targets_list)) {
    # keep as-is (must be character vectors)
    ok <- all(vapply(targets_list, function(x) is.character(x) && length(x) >= 1, logical(1)))
    if (!ok) stop("Each element of `targets_list` must be a non-empty character vector of labels.")
    if (is.null(names(targets_list)) || any(names(targets_list) == ""))
      stop("`targets_list` must be a *named* list; names label the target groups.")
  } else {
    stop("`targets_list` must be a character vector or a named list of character vectors.")
  }

  # ---- compute a single target group ----------------------------------------
  do_target <- function(gname, labvec) {
    tgt <- dplyr::filter(pts, .data$assignment %in% labvec)
    if (!nrow(tgt)) {
      d_um <- rep(NA_real_, nrow(src))
    } else {
      dmat <- sf::st_distance(src, tgt)  # |src| x |tgt|
      # avoid zero self-match if same cell ends up in both sets
      if ("cell_id" %in% names(tgt)) {
        col_match <- match(src$cell_id, tgt$cell_id, nomatch = 0)
        bad <- which(col_match > 0)
        if (length(bad)) dmat[cbind(bad, col_match[bad])] <- Inf
      }
      d_px <- apply(as.matrix(dmat), 1, min)
      d_um <- as.numeric(d_px) * px2um
    }
    s <- dplyr::bind_cols(tibble::tibble(target = gname), .summ(d_um))
    pc <- tibble::tibble(cell_id = src$cell_id, dist_um = d_um, target = gname)
    list(summary = s, per_cell = pc)
  }

  # ---- run targets (parallel if requested) -----------------------------------
  group_names <- names(targets_list)
  res_list <- if (n_cores > 1) {
    parallel::mclapply(group_names, function(nm) do_target(nm, targets_list[[nm]]),
                       mc.cores = n_cores)
  } else {
    lapply(group_names, function(nm) do_target(nm, targets_list[[nm]]))
  }

  sum_tgt <- dplyr::bind_rows(lapply(res_list, `[[`, "summary"))
  pc_tgt  <- dplyr::bind_rows(lapply(res_list, `[[`, "per_cell"))

  # ---- optionally add boundary in same run -----------------------------------
  if (isTRUE(boundary)) {
    b <- do_boundary()
    sum_all <- dplyr::bind_rows(sum_tgt, b$summary)
    pc_all  <- if (return_per_cell) dplyr::bind_rows(pc_tgt, b$per_cell) else NULL
  } else {
    sum_all <- sum_tgt
    pc_all  <- if (return_per_cell) pc_tgt else NULL
  }

  if (!return_per_cell) return(sum_all)
  list(summary = sum_all, per_cell = pc_all)
}


# Proximity (Fixed-Radius, Permutation Test) ---------------------------------

.nn_radius <- function(coords, r_px) {
  id <- dbscan::frNN(coords, eps = r_px, sort = FALSE)$id
  lapply(seq_along(id), function(i) setdiff(id[[i]], i))
}
.count_matrix <- function(labels, nn, levels) {
  idx <- match(labels, levels); K <- length(levels)
  M <- matrix(0L, K, K, dimnames = list(levels, levels))
  for (i in seq_along(nn)) {
    j <- nn[[i]]; if (!length(j)) next
    tab <- table(factor(idx[j], levels = 1:K))
    M[, idx[i]] <- M[, idx[i]] + as.integer(tab)
  }
  M
}

codex_proximity_simple <- function(centroids_df,   # x,y,label plus optional STRAT
                                   radii_um   = c(20,50,100),
                                   pix_per_um = 1/0.5067,
                                   stratum    = c("none","nest","zone"),
                                   stratum_col= NULL,
                                   n_perm     = 2000,
                                   n_cores    = 1) {
  .codex_use_planar()
  stopifnot(all(c("x","y","label") %in% names(centroids_df)))
  stratum <- match.arg(stratum)
  if (stratum != "none") {
    stopifnot(!is.null(stratum_col), stratum_col %in% names(centroids_df))
    centroids_df$STRAT <- factor(centroids_df[[stratum_col]])
  } else centroids_df$STRAT <- factor("all")

  lv <- sort(unique(centroids_df$label))
  results <- vector("list", length(radii_um)); names(results) <- paste0(radii_um,"um")

  for (r_um in radii_um) {
    r_px <- r_um * pix_per_um
    per_strat <- lapply(split(centroids_df, centroids_df$STRAT), function(D) {
      coords <- as.matrix(D[,c("x","y")]); labs <- as.character(D$label)
      nn  <- .nn_radius(coords, r_px)
      C   <- .count_matrix(labs, nn, lv); diag(C) <- 0L
      P   <- sweep(C, 2, pmax(colSums(C),1L), "/")
      if (n_perm > 0) {
        seeds <- sample.int(.Machine$integer.max, n_perm)
        one <- function(seed) {
          set.seed(seed); lp <- sample(labs)
          C2 <- .count_matrix(lp, nn, lv); diag(C2) <- 0L
          sweep(C2, 2, pmax(colSums(C2),1L), "/")
        }
        permL <- if (n_cores>1) parallel::mclapply(seeds, one, mc.cores = n_cores) else lapply(seeds, one)
        permA <- simplify2array(permL)  # KxKxB
      } else permA <- array(NA_real_, dim = c(length(lv),length(lv),0L))
      list(Pobs = P, perm = permA)
    })
    Pobs_sum <- Reduce(`+`, lapply(per_strat, `[[`, "Pobs"))
    perm_sum <- Reduce(function(a,b){
      if (length(dim(a))==0) return(b)
      if (length(dim(b))==0) return(a)
      stopifnot(identical(dim(a), dim(b))); a+b
    }, lapply(per_strat, `[[`, "perm"))

    K <- length(lv); obs_vec <- as.vector(Pobs_sum)
    rows <- rep(lv, times = K); cols <- rep(lv, each = K)

    if (length(dim(perm_sum))==0 || dim(perm_sum)[3]==0L) {
      results[[paste0(r_um,"um")]] <- tibble::tibble(
        radius_um=r_um, Ref_State=cols, Neighb_State=rows, Prop_Neighb=obs_vec,
        Z_Score=NA_real_, P_val=NA_real_)
      next
    }
    perm_flat <- apply(perm_sum, 3, as.vector)  # (K*K) x B
    mu <- rowMeans(perm_flat, na.rm = TRUE)
    sd <- apply(perm_flat, 1, stats::sd, na.rm = TRUE)
    ge <- sweep(perm_flat, 1, obs_vec, `>=`)
    gt <- sweep(perm_flat, 1, obs_vec, `>`)
    R <- (rowSums(ge, na.rm = TRUE)+1)/(ncol(perm_flat)+1)
    L <- (ncol(perm_flat)-rowSums(gt, na.rm = TRUE)+1)/(ncol(perm_flat)+1)
    p <- pmin(R, L)
    z <- ifelse(sd>0, (obs_vec - mu)/sd, 0)
    col_has <- colSums(Pobs_sum) > 0
    mask <- rep(col_has, each = K); z[!mask] <- NA_real_; p[!mask] <- NA_real_

    results[[paste0(r_um,"um")]] <- tibble::tibble(
      radius_um=r_um, Ref_State=cols, Neighb_State=rows,
      Prop_Neighb=obs_vec, Z_Score=z, P_val=p)
  }
  dplyr::bind_rows(results)
}


# Tumour Zones (Edge / Mid / Core) & Zone Heatmaps ---------------------------

.codex_as_sfc <- function(x) {
  if (inherits(x, "sf")) sf::st_geometry(x) else x
}

.codex_make_poly <- function(g) {
  g <- .codex_as_sfc(g)
  g <- sf::st_make_valid(g)

  # If geometry collections exist, extract polygons from them
  if (any(sf::st_geometry_type(g) %in% "GEOMETRYCOLLECTION")) {
    g2 <- sf::st_collection_extract(g, "POLYGON")
    if (length(g2) > 0) g <- g2
  }

  # Cast polygons consistently (POLYGON -> MULTIPOLYGON; MULTIPOLYGON stays)
  suppressWarnings({
    if (length(g) > 0) g <- sf::st_cast(g, "MULTIPOLYGON")
  })

  g
}


# Offsets outward/inward from boundary (positive dist grows outward; we need inward only)
codex_split_tumour_zones <- function(nest_poly, pix_per_um, edge_um = 100, mid_um = 200) {
  .codex_use_planar()

  nest <- .codex_make_poly(nest_poly)

  edge_px <- edge_um * pix_per_um
  mid_px  <- mid_um  * pix_per_um

  # inward buffers (negative distance)
  inner1 <- sf::st_buffer(nest, dist = -edge_px) %>% sf::st_make_valid()
  inner2 <- sf::st_buffer(nest, dist = -mid_px)  %>% sf::st_make_valid()

  # Make sure these remain polygonal types (and safe if empty)
  inner1 <- .codex_make_poly(inner1)
  inner2 <- .codex_make_poly(inner2)

  # Element-wise zone construction
  edge <- sf::st_difference(nest, inner1) %>% sf::st_make_valid()
  mid  <- sf::st_difference(inner1, inner2) %>% sf::st_make_valid()
  core <- inner2

  zones <- list(edge = edge, mid = mid, core = core)
  attr(zones, "nest_poly") <- nest
  attr(zones, "inner1")    <- inner1
  attr(zones, "inner2")    <- inner2
  zones
}

# Boundaries for plotting
codex_zone_boundaries <- function(zones) {
  .codex_use_planar()

  nest  <- attr(zones, "nest_poly")
  inner1 <- attr(zones, "inner1")
  inner2 <- attr(zones, "inner2")

  if (is.null(nest) || is.null(inner1) || is.null(inner2))
    stop("`zones` must come from codex_split_tumour_zones() (attributes missing).")

  out <- list(
    edge = sf::st_cast(sf::st_boundary(nest),  "MULTILINESTRING"),
    mid  = sf::st_cast(sf::st_boundary(inner1),"MULTILINESTRING"),
    core = sf::st_cast(sf::st_boundary(inner2),"MULTILINESTRING")
  )

  # Drop empty boundaries (e.g. if inner buffers collapse for small nests)
  out <- out[!vapply(out, function(x) length(x) == 0 || all(sf::st_is_empty(x)), logical(1))]
  out
}


codex_add_zone_boundaries_to_plot <- function(
    p,
    zone_bounds,
    boundary_colours    = c(edge = "white", mid = "orange", core = "red"),
    boundary_linewidths = c(edge = 1.2,  mid = 0.8,   core = 0.8),
    boundary_linetypes  = c(edge = "solid", mid = "solid", core = "solid")
) {
  for (nm in names(zone_bounds)) {
    p <- p + ggplot2::geom_sf(
      data = sf::st_sf(boundary = nm, geometry = zone_bounds[[nm]]),
      inherit.aes = FALSE,
      colour = boundary_colours[[nm]],
      linewidth = boundary_linewidths[[nm]],
      linetype = boundary_linetypes[[nm]]
    )
  }
  p
}


# Assign each cell to: Non-malignant / Tumor-edge / Mid-tumor / Tumor-core
codex_assign_cells_to_tumour_zones <- function(cells_tbl_core, zones, use_intersects = TRUE) {
  .codex_use_planar()

  stopifnot(all(c("cell_id", "centroid_x", "centroid_y") %in% names(cells_tbl_core)))

  nest_poly <- attr(zones, "nest_poly")
  if (is.null(nest_poly)) {
    stop("`zones` must come from codex_split_tumour_zones() (nest_poly attr missing).")
  }

  # Keep CRS consistent (even if NA)
  crs_use <- sf::st_crs(nest_poly)

  pts <- sf::st_as_sf(
    cells_tbl_core,
    coords = c("centroid_x", "centroid_y"),
    remove = FALSE,
    crs = crs_use
  )

  make_valid_nonempty <- function(g) {
    g <- .codex_as_sfc(g)
    g <- sf::st_make_valid(g)
    if (length(g) == 0) return(g)
    g[!sf::st_is_empty(g)]
  }

  nest_poly <- make_valid_nonempty(nest_poly)
  z_edge    <- make_valid_nonempty(zones$edge)
  z_mid     <- make_valid_nonempty(zones$mid)
  z_core    <- make_valid_nonempty(zones$core)

  rel <- if (isTRUE(use_intersects)) sf::st_intersects else sf::st_within

  in_any <- function(poly) {
    if (length(poly) == 0) return(rep(FALSE, nrow(pts)))
    lengths(rel(pts, poly)) > 0
  }

  in_nest <- in_any(nest_poly)
  in_core <- in_any(z_core)
  in_mid  <- in_any(z_mid)
  in_edge <- in_any(z_edge)

  Zone <- dplyr::case_when(
    !in_nest ~ "Non-malignant",
    in_core  ~ "Tumor-core",
    in_mid   ~ "Mid-tumor",
    in_edge  ~ "Tumor-edge",
    TRUE     ~ "Tumor-edge"   # inside nest but missed by zone polygons (rare)
  )

  tibble::tibble(
    cell_id = cells_tbl_core$cell_id,
    Zone    = Zone
  )
}


codex_prepare_zone_heatmap_allcells <- function(
    rmat,
    core_obj,
    zones,
    cells_meta,
    markers,
    zone_levels = c("Non-malignant", "Tumor-edge", "Mid-tumor", "Tumor-core")
) {
  stopifnot(is.matrix(rmat))
  stopifnot(all(markers %in% rownames(rmat)))
  stopifnot(is.data.frame(cells_meta))
  stopifnot(all(c("cell_id", "cell_name") %in% colnames(cells_meta)))
  stopifnot(!is.null(core_obj$prep$cells_tbl))

  # ---- 1) Assign EVERY core cell to a tumour zone (expensive part) ----
  zone_df <- codex_assign_cells_to_tumour_zones(
    cells_tbl_core = core_obj$prep$cells_tbl,
    zones          = zones,
    use_intersects = TRUE
  )

  if (!all(c("cell_id", "Zone") %in% colnames(zone_df))) {
    stop("`codex_assign_cells_to_tumour_zones()` must return columns: cell_id, Zone")
  }

  # Normalize Zone labels to canonical ones
  zone_df <- zone_df %>%
    dplyr::mutate(
      Zone = dplyr::case_when(
        Zone %in% c("non_malignant", "non-malignant", "NonMalignant") ~ "Non-malignant",
        Zone %in% c("edge", "tumor_edge", "Tumor-edge", "tumour_edge") ~ "Tumor-edge",
        Zone %in% c("mid", "mid_tumor", "Mid-tumor", "mid-tumour") ~ "Mid-tumor",
        Zone %in% c("core", "tumor_core", "Tumor-core", "tumour_core") ~ "Tumor-core",
        TRUE ~ as.character(Zone)
      ),
      Zone = factor(Zone, levels = zone_levels)
    )

  # ---- 2) Add cell_name (match to rmat colnames) ----
  id_map <- cells_meta %>%
    dplyr::select(cell_id, cell_name) %>%
    dplyr::distinct()

  zone_df <- zone_df %>%
    dplyr::left_join(id_map, by = "cell_id")

  # ---- 3) Filter to cells that exist in rmat (NOT downsampling) ----
  zone_df <- zone_df %>%
    dplyr::filter(!is.na(cell_name), cell_name %in% colnames(rmat))

  if (nrow(zone_df) == 0) {
    stop("After mapping to cell_name and filtering to rmat, no cells remain. Check cell_id↔cell_name mapping.")
  }

  # ---- 4) Build marker-by-cell matrix (selected markers only, all cells) ----
  cell_order_all <- zone_df$cell_name
  mat_all <- rmat[markers, cell_order_all, drop = FALSE]  # markers × cells

  # ---- 5) Prepare annotations (Zone + assignment) ----
  core_cells <- core_obj$prep$cells_tbl %>%
    dplyr::select(cell_id, assignment) %>%
    dplyr::distinct()

  ann_df <- zone_df %>%
    dplyr::left_join(core_cells, by = "cell_id") %>%
    dplyr::mutate(assignment = ifelse(is.na(assignment), "Undefined", assignment)) %>%
    dplyr::select(cell_name, Zone, assignment)

  # ---- 6) Long df for ggplot (still all cells) ----
  heat_long_all <- as.data.frame(t(mat_all), stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column("cell_name") %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(markers),
      names_to = "Marker",
      values_to = "Score"
    ) %>%
    dplyr::left_join(ann_df, by = "cell_name") %>%
    dplyr::mutate(
      Marker = factor(Marker, levels = markers),
      Zone   = factor(Zone, levels = zone_levels)
    )

  list(
    core_id        = core_obj$core_id,
    markers        = markers,
    zone_levels    = zone_levels,
    zone_df_all    = zone_df,
    ann_df_all     = ann_df,
    cell_order_all = cell_order_all,
    mat_all        = mat_all,
    heat_long_all  = heat_long_all
  )
}


codex_downsample_normalize_plot_zone_heatmap <- function(
    prep,                       # output of codex_prepare_zone_heatmap_allcells()
    max_cells_per_zone = 2000,  # start smaller for visual debugging; increase later
    prioritize_informative = TRUE,
    informativeness_fun = c("max", "mean", "l2"),
    normalization = c("marker_robust_z", "marker_zscore", "marker_minmax", "none"),
    clip_z = 3,
    seed = 1,
    order_within_zone = c("assignment", "informativeness", "none"),
    assignment_order = NULL,
    show_assignment_bar = TRUE,     # <--- NEW
    hide_x_text = TRUE
) {
  informativeness_fun <- match.arg(informativeness_fun)
  normalization <- match.arg(normalization)
  order_within_zone <- match.arg(order_within_zone)

  stopifnot(is.list(prep), !is.null(prep$mat_all), !is.null(prep$ann_df_all))
  stopifnot(is.matrix(prep$mat_all))

  mat_all <- prep$mat_all
  ann_all <- prep$ann_df_all
  markers <- prep$markers
  zone_levels <- prep$zone_levels

  # ---- make joins robust ----
  ann_all <- ann_all %>%
    dplyr::mutate(
      cell_name  = as.character(cell_name),
      Zone       = factor(as.character(Zone), levels = zone_levels),
      assignment = as.character(assignment)
    ) %>%
    dplyr::filter(!is.na(cell_name), !is.na(Zone))

  # ---- 1) informativeness per cell ----
  inf_score <- switch(
    informativeness_fun,
    max  = apply(mat_all, 2, max,  na.rm = TRUE),
    mean = apply(mat_all, 2, mean, na.rm = TRUE),
    l2   = sqrt(colSums(mat_all^2, na.rm = TRUE))
  )

  inf_df <- tibble::tibble(
    cell_name = colnames(mat_all),
    informativeness = as.numeric(inf_score)
  ) %>%
    dplyr::left_join(ann_all, by = "cell_name") %>%
    dplyr::filter(!is.na(Zone))

  # ---- 2) choose cells per zone ----
  set.seed(seed)

  pick_zone <- function(dfz) {
    if (nrow(dfz) == 0) return(character(0))
    if (is.finite(max_cells_per_zone) && nrow(dfz) > max_cells_per_zone) {
      if (prioritize_informative) {
        dfz <- dfz %>%
          dplyr::arrange(dplyr::desc(informativeness)) %>%
          dplyr::slice_head(n = max_cells_per_zone)
      } else {
        dfz <- dfz %>% dplyr::slice_sample(n = max_cells_per_zone)
      }
    }
    dfz$cell_name
  }

  chosen_cells <- unlist(lapply(zone_levels, function(z) {
    pick_zone(inf_df %>% dplyr::filter(Zone == z))
  }), use.names = FALSE)

  chosen_cells <- unique(chosen_cells)
  if (!length(chosen_cells)) stop("No cells selected after downsampling.")

  mat_sub <- mat_all[, chosen_cells, drop = FALSE]  # markers × cells

  # ---- 3) normalize (per marker across ALL chosen cells) ----
  mat_norm <- mat_sub

  if (normalization == "marker_zscore") {
    mu  <- rowMeans(mat_sub, na.rm = TRUE)
    sdv <- apply(mat_sub, 1, sd, na.rm = TRUE)
    sdv[sdv == 0 | is.na(sdv)] <- 1
    mat_norm <- sweep(mat_sub, 1, mu, "-")
    mat_norm <- sweep(mat_norm, 1, sdv, "/")
    mat_norm[mat_norm >  clip_z] <-  clip_z
    mat_norm[mat_norm < -clip_z] <- -clip_z

  } else if (normalization == "marker_robust_z") {
    med  <- apply(mat_sub, 1, stats::median, na.rm = TRUE)
    madv <- apply(mat_sub, 1, stats::mad, na.rm = TRUE)
    madv[madv == 0 | is.na(madv)] <- 1
    mat_norm <- sweep(mat_sub, 1, med, "-")
    mat_norm <- sweep(mat_norm, 1, madv, "/")
    mat_norm[mat_norm >  clip_z] <-  clip_z
    mat_norm[mat_norm < -clip_z] <- -clip_z

  } else if (normalization == "marker_minmax") {
    mn  <- apply(mat_sub, 1, min, na.rm = TRUE)
    mx  <- apply(mat_sub, 1, max, na.rm = TRUE)
    rng <- mx - mn
    rng[rng == 0 | is.na(rng)] <- 1
    mat_norm <- sweep(mat_sub, 1, mn, "-")
    mat_norm <- sweep(mat_norm, 1, rng, "/")
  }

  # ---- 4) build ordering table (THIS is what makes structure visible) ----
  ann_sub <- inf_df %>%
    dplyr::filter(cell_name %in% chosen_cells) %>%
    dplyr::mutate(
      Zone = factor(as.character(Zone), levels = zone_levels)
    )

  if (!is.null(assignment_order)) {
    ann_sub <- ann_sub %>%
      dplyr::mutate(assignment = factor(assignment, levels = assignment_order))
  } else {
    # keep deterministic order if user didn't specify
    ann_sub <- ann_sub %>%
      dplyr::mutate(assignment = factor(assignment, levels = sort(unique(assignment))))
  }

  if (order_within_zone == "assignment") {
    ann_sub <- ann_sub %>%
      dplyr::group_by(Zone) %>%
      dplyr::arrange(assignment, dplyr::desc(informativeness), .by_group = TRUE) %>%
      dplyr::mutate(x_in_zone = dplyr::row_number()) %>%
      dplyr::ungroup()

  } else if (order_within_zone == "informativeness") {
    ann_sub <- ann_sub %>%
      dplyr::group_by(Zone) %>%
      dplyr::arrange(dplyr::desc(informativeness), .by_group = TRUE) %>%
      dplyr::mutate(x_in_zone = dplyr::row_number()) %>%
      dplyr::ungroup()

  } else {
    ann_sub <- ann_sub %>%
      dplyr::group_by(Zone) %>%
      dplyr::mutate(x_in_zone = dplyr::row_number()) %>%
      dplyr::ungroup()
  }

  # Ensure matrix columns match final ordering (zone-wise)
  cell_order <- ann_sub %>%
    dplyr::arrange(Zone, x_in_zone) %>%
    dplyr::pull(cell_name)

  mat_norm <- mat_norm[, cell_order, drop = FALSE]

  # ---- 5) long df for heatmap with TRUE faceting by Zone ----
  df_heat <- as.data.frame(t(mat_norm), stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column("cell_name") %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(markers),
      names_to = "Marker",
      values_to = "Score"
    ) %>%
    dplyr::left_join(
      ann_sub %>% dplyr::select(cell_name, Zone, assignment, x_in_zone),
      by = "cell_name"
    ) %>%
    dplyr::mutate(
      Marker = factor(Marker, levels = markers),
      Zone   = factor(Zone, levels = zone_levels)
    )

  # ---- 6) plotting scales ----
  if (normalization %in% c("marker_zscore", "marker_robust_z")) {
    fill_scale <- ggplot2::scale_fill_gradient2(
      low = "#053061", mid = "white", high = "#67001f",
      midpoint = 0, name = "Score\n(z)"
    )
  } else {
    fill_scale <- ggplot2::scale_fill_gradientn(
      colours = c("#f7f7f7","#fddbc7","#f4a582","#d6604d","#b2182b","#67001f"),
      name = "Score"
    )
  }

  p_heat <- ggplot2::ggplot(df_heat, ggplot2::aes(x = x_in_zone, y = Marker, fill = Score)) +
    ggplot2::geom_tile(width = 1, height = 1) +
    ggplot2::facet_grid(. ~ Zone, scales = "free_x", space = "free_x") +
    fill_scale +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(face = "bold", size = 13),
      axis.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = if (hide_x_text) ggplot2::element_blank() else ggplot2::element_text(size = 6),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = paste0(prep$core_id, " — tumour zonation heatmap"))

  # ---- optional assignment bar ----
  p <- p_heat
  p_assign <- NULL

  if (isTRUE(show_assignment_bar)) {
    df_bar <- ann_sub %>%
      dplyr::select(Zone, x_in_zone, assignment)

    p_assign <- ggplot2::ggplot(df_bar, ggplot2::aes(x = x_in_zone, y = 1, fill = assignment)) +
      ggplot2::geom_tile(width = 1, height = 1) +
      ggplot2::facet_grid(. ~ Zone, scales = "free_x", space = "free_x") +
      ggplot2::theme_void(base_size = 12) +
      ggplot2::theme(
        strip.text.x = ggplot2::element_blank(),
        legend.position = "right"
      ) +
      ggplot2::labs(fill = "Assignment")

    if (requireNamespace("patchwork", quietly = TRUE)) {
      p <- p_assign / p_heat + patchwork::plot_layout(heights = c(0.35, 5))
    } else {
      message("Install `patchwork` to stack the assignment bar above the heatmap. Returning heatmap only.")
      p <- p_heat
    }
  }

  list(
    core_id = prep$core_id,
    chosen_cells = chosen_cells,
    ann_sub = ann_sub,
    cell_order = cell_order,
    mat_norm = mat_norm,
    df_heat = df_heat,
    plot_heat = p_heat,
    plot_assign = p_assign,
    plot = p
  )
}


codex_plot_celltype_zone_enrichment <- function(
    prep,  # output of codex_prepare_zone_heatmap_allcells()
    cell_types = c(
      "Epithelial","Senescence","Cell-Cycle","Hypoxia","pEMT",
      "Fibroblast","Endothelial","Neutrophil","Monocyte",
      "M1-like_Macrophage","M2-like_Macrophage","Tctx-cell","B-cell"
    ),
    zone_levels = NULL,  # default: prep$zone_levels
    # core metric before row-normalization:
    metric = c("log2oe", "frac_in_zone"),
    pseudocount = 0.5,   # for log2(O/E); 0.5 is a good default for counts
    # row-wise normalization to make it cell-type–centric:
    row_normalize = c("zscore", "center", "none"),
    clip = 3,            # clip after row_normalize (keeps palette readable)
    order_rows_by = c("max_zone", "none"),
    show_best_zone = TRUE,
    show_values = FALSE,
    value_digits = 2
) {
  metric <- match.arg(metric)
  row_normalize <- match.arg(row_normalize)
  order_rows_by <- match.arg(order_rows_by)

  stopifnot(is.list(prep), !is.null(prep$ann_df_all))
  ann <- prep$ann_df_all

  if (is.null(zone_levels)) zone_levels <- prep$zone_levels
  stopifnot(!is.null(zone_levels), length(zone_levels) > 0)

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
  })

  # ---- 1) counts table (cell_type x zone) ----
  df0 <- ann %>%
    dplyr::transmute(
      Zone = factor(as.character(Zone), levels = zone_levels),
      cell_type = as.character(assignment)
    ) %>%
    dplyr::filter(!is.na(Zone), !is.na(cell_type)) %>%
    dplyr::filter(cell_type %in% cell_types)

  if (nrow(df0) == 0) stop("No rows after filtering to requested cell_types. Check naming/levels.")

  cnt <- df0 %>%
    dplyr::count(cell_type, Zone, name = "obs") %>%
    tidyr::complete(
      cell_type = cell_types,
      Zone = factor(zone_levels, levels = zone_levels),
      fill = list(obs = 0)
    )

  # totals
  N <- sum(cnt$obs)
  ct_tot <- cnt %>% dplyr::group_by(cell_type) %>% dplyr::summarise(ct_n = sum(obs), .groups = "drop")
  z_tot  <- cnt %>% dplyr::group_by(Zone)      %>% dplyr::summarise(z_n  = sum(obs), .groups = "drop")

  cnt <- cnt %>%
    dplyr::left_join(ct_tot, by = "cell_type") %>%
    dplyr::left_join(z_tot,  by = "Zone") %>%
    dplyr::mutate(
      exp = ifelse(N > 0, (ct_n * z_n) / N, 0)
    )

  # ---- 2) base metric ----
  if (metric == "log2oe") {
    cnt <- cnt %>%
      dplyr::mutate(value_raw = log2((obs + pseudocount) / (exp + pseudocount)),
                    value_label = "log2(O/E)")
  } else { # "frac_in_zone"
    # fraction of this zone that is this cell_type (zone-centric if used raw)
    # we only use it if you still want it, but typically log2oe is better.
    cnt <- cnt %>%
      dplyr::mutate(value_raw = ifelse(z_n > 0, obs / z_n, 0),
                    value_label = "Fraction in zone")
  }

  # ---- 3) row-wise normalization to make it cell-type–centric ----
  cnt <- cnt %>%
    dplyr::group_by(cell_type) %>%
    dplyr::mutate(
      value = dplyr::case_when(
        row_normalize == "none"   ~ value_raw,
        row_normalize == "center" ~ value_raw - mean(value_raw, na.rm = TRUE),
        row_normalize == "zscore" ~ {
          s <- stats::sd(value_raw, na.rm = TRUE)
          if (is.na(s) || s == 0) 0 else (value_raw - mean(value_raw, na.rm = TRUE)) / s
        }
      ),
      value_label2 = dplyr::case_when(
        row_normalize == "none"   ~ value_label,
        row_normalize == "center" ~ paste0(value_label, " (row-centered)"),
        row_normalize == "zscore" ~ paste0(value_label, " (row z)")
      )
    ) %>%
    dplyr::ungroup()

  # clip after normalization
  if (!is.null(clip) && is.finite(clip)) {
    cnt <- cnt %>% dplyr::mutate(value = pmax(pmin(value, clip), -clip))
  }

  # ---- 4) optionally order rows by which zone is "best" ----
  best <- NULL
  if (order_rows_by == "max_zone") {
    best <- cnt %>%
      dplyr::group_by(cell_type) %>%
      dplyr::slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(cell_type, best_zone = Zone, best_value = value)

    # order by best_zone then best_value, so you can read patterns like:
    # "these cell types peak in Non-malignant", then "... peak in Tumor-core", etc.
    row_levels <- best %>%
      dplyr::arrange(best_zone, dplyr::desc(best_value)) %>%
      dplyr::pull(cell_type)

    cnt <- cnt %>%
      dplyr::mutate(cell_type = factor(cell_type, levels = rev(row_levels)))
  } else {
    cnt <- cnt %>%
      dplyr::mutate(cell_type = factor(cell_type, levels = rev(cell_types)))
  }

  cnt <- cnt %>%
    dplyr::mutate(Zone = factor(Zone, levels = zone_levels))

  # ---- 5) plot ----
  fill_scale <- ggplot2::scale_fill_gradient2(
    low = "#053061", mid = "white", high = "#67001f",
    midpoint = 0,
    name = cnt$value_label2[1]
  )

  p <- ggplot2::ggplot(cnt, ggplot2::aes(x = Zone, y = cell_type, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    fill_scale +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold", size = 12),
      axis.text.y = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    ) +
    ggplot2::labs(title = paste0(prep$core_id, " — cell-type enrichment across zones"))

  if (isTRUE(show_values)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = round(value, value_digits)),
      size = 3
    )
  }

  # optionally mark each row's best zone
  if (isTRUE(show_best_zone)) {
    best_pts <- cnt %>%
      dplyr::group_by(cell_type) %>%
      dplyr::slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()

    p <- p + ggplot2::geom_point(
      data = best_pts,
      ggplot2::aes(x = Zone, y = cell_type),
      inherit.aes = FALSE,
      shape = 21, fill = "black", color = "black", size = 2.2
    )
  }

  list(df = cnt, best = best, plot = p)
}


codex_plot_celltype_zone_enrichment_dual <- function(
    prep,
    cell_types,
    zone_levels = prep$zone_levels,
    drop_assignments = c("Undefined", "Myeloid"),   # adjust if you want to keep Myeloid
    eps = 1e-4,                                    # small stabilizer in proportion space
    smooth_alpha = 0.5,                             # Dirichlet smoothing on proportions
    fill_mode = c("row_z", "raw"),                  # what to color by
    clip_z = 2.5,
    zone_for_row_norm = c("all", "tumour_only"),    # row normalization across zones
    tumour_zones = c("Tumor-edge", "Mid-tumor", "Tumor-core"),
    mark_row_max = TRUE,
    mark_col_topk = 1,                              # 0 disables
    order_cell_types = TRUE,                        # order by zone of max enrichment
    base_size = 18
) {
  fill_mode <- match.arg(fill_mode)
  zone_for_row_norm <- match.arg(zone_for_row_norm)

  stopifnot(is.list(prep), !is.null(prep$ann_df_all))
  ann <- prep$ann_df_all

  # ---- 1) counts per (cell_type, zone) ----
  df <- ann %>%
    dplyr::transmute(
      Zone = as.character(Zone),
      cell_type = assignment
    ) %>%
    dplyr::filter(!is.na(Zone), Zone %in% zone_levels) %>%
    dplyr::filter(cell_type %in% cell_types)

  if (length(drop_assignments)) {
    df <- df %>% dplyr::filter(!cell_type %in% drop_assignments)
  }

  # complete grid so missing combos are 0
  counts <- df %>%
    dplyr::count(cell_type, Zone, name = "n") %>%
    tidyr::complete(
      cell_type = cell_types[cell_types %in% unique(df$cell_type)],
      Zone = zone_levels,
      fill = list(n = 0)
    )

  # zone sizes (for sanity/debug)
  zone_sizes <- counts %>%
    dplyr::group_by(Zone) %>%
    dplyr::summarise(N_zone = sum(n), .groups = "drop")

  # ---- 2) proportions with smoothing (Dirichlet-ish) ----
  # p(z,ct) = (n + alpha) / (N_zone + alpha*K)
  K <- length(unique(counts$cell_type))
  counts <- counts %>%
    dplyr::left_join(zone_sizes, by = "Zone") %>%
    dplyr::mutate(
      p_obs = (n + smooth_alpha) / (N_zone + smooth_alpha * K)
    )

  # ---- 3) zone-balanced expected: mean_z p_obs (unweighted across zones) ----
  exp_bal <- counts %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(p_exp_bal = mean(p_obs), .groups = "drop")

  counts <- counts %>%
    dplyr::left_join(exp_bal, by = "cell_type") %>%
    dplyr::mutate(
      log2OE_bal = log2((p_obs + eps) / (p_exp_bal + eps))
    )

  # ---- 4) row-wise normalization (cell-type centric) ----
  zones_for_row <- if (zone_for_row_norm == "tumour_only") {
    intersect(zone_levels, tumour_zones)
  } else {
    zone_levels
  }

  row_stats <- counts %>%
    dplyr::filter(Zone %in% zones_for_row) %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(mu = mean(log2OE_bal), sd = stats::sd(log2OE_bal), .groups = "drop") %>%
    dplyr::mutate(sd = dplyr::if_else(is.na(sd) | sd == 0, 1, sd))

  counts <- counts %>%
    dplyr::left_join(row_stats, by = "cell_type") %>%
    dplyr::mutate(
      row_z = (log2OE_bal - mu) / sd,
      row_z = pmax(pmin(row_z, clip_z), -clip_z)
    )

  # ---- 5) column-wise normalization (zone centric) ----
  col_stats <- counts %>%
    dplyr::group_by(Zone) %>%
    dplyr::summarise(mu_c = mean(log2OE_bal), sd_c = stats::sd(log2OE_bal), .groups = "drop") %>%
    dplyr::mutate(sd_c = dplyr::if_else(is.na(sd_c) | sd_c == 0, 1, sd_c))

  counts <- counts %>%
    dplyr::left_join(col_stats, by = "Zone") %>%
    dplyr::mutate(
      col_z = (log2OE_bal - mu_c) / sd_c,
      col_z = pmax(pmin(col_z, clip_z), -clip_z)
    )

  # ---- 6) pick row-max (best zone per cell type) ----
  row_max <- counts %>%
    dplyr::filter(Zone %in% zones_for_row) %>%
    dplyr::group_by(cell_type) %>%
    dplyr::slice_max(order_by = row_z, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  # ---- 7) pick col-topk (best cell types per zone) ----
  col_top <- NULL
  if (mark_col_topk > 0) {
    col_top <- counts %>%
      dplyr::group_by(Zone) %>%
      dplyr::slice_max(order_by = col_z, n = mark_col_topk, with_ties = FALSE) %>%
      dplyr::ungroup()
  }

  # ---- 8) ordering (optional) ----
  counts <- counts %>%
    dplyr::mutate(
      Zone = factor(Zone, levels = zone_levels)
    )

  if (order_cell_types) {
    # order rows by zone of row-max, then by strength
    ord <- row_max %>%
      dplyr::mutate(Zone = factor(Zone, levels = zone_levels)) %>%
      dplyr::arrange(Zone, dplyr::desc(row_z)) %>%
      dplyr::pull(cell_type) %>%
      unique()
    counts <- counts %>% dplyr::mutate(cell_type = factor(cell_type, levels = rev(ord)))
    row_max <- row_max %>% dplyr::mutate(cell_type = factor(cell_type, levels = rev(ord)))
    if (!is.null(col_top)) col_top <- col_top %>% dplyr::mutate(cell_type = factor(cell_type, levels = rev(ord)))
  } else {
    counts <- counts %>% dplyr::mutate(cell_type = factor(cell_type, levels = rev(cell_types)))
    row_max <- row_max %>% dplyr::mutate(cell_type = factor(cell_type, levels = rev(cell_types)))
    if (!is.null(col_top)) col_top <- col_top %>% dplyr::mutate(cell_type = factor(cell_type, levels = rev(cell_types)))
  }

  fill_var <- if (fill_mode == "row_z") "row_z" else "log2OE_bal"

  # ---- 9) plot ----
  p <- ggplot2::ggplot(
    counts,
    ggplot2::aes(x = Zone, y = cell_type, fill = .data[[fill_var]])
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::scale_fill_gradient2(
      low = "#053061", mid = "white", high = "#67001f",
      midpoint = 0,
      name = if (fill_mode == "row_z") "log2(O/E)\n(row z)" else "log2(O/E)\n(zone-balanced)"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, face = "bold"),
      axis.text.y = ggplot2::element_text(face = "plain"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    ) +
    ggplot2::labs(title = paste0(prep$core_id, " — cell-type enrichment across zones"))

  if (isTRUE(mark_row_max)) {
    p <- p + ggplot2::geom_point(
      data = row_max,
      ggplot2::aes(x = Zone, y = cell_type),
      inherit.aes = FALSE,
      shape = 16, size = 3, color = "black"
    )
  }

  if (!is.null(col_top) && nrow(col_top) > 0) {
    p <- p + ggplot2::geom_point(
      data = col_top,
      ggplot2::aes(x = Zone, y = cell_type),
      inherit.aes = FALSE,
      shape = 17, size = 3.2, color = "black"
    )
  }

  list(
    core_id = prep$core_id,
    zone_sizes = zone_sizes,
    df = counts,
    row_max = row_max,
    col_top = col_top,
    plot = p
  )
}


# codex_core_cellid_to_cellname <- function(core_obj, cells_meta) {
#   core_ids <- core_obj$prep$cells_tbl$cell_id
#
#   cells_meta |>
#     dplyr::filter(.data$cell_id %in% core_ids) |>
#     dplyr::distinct(.data$cell_id, .data$cell_name)
# }
#
#
#
# codex_build_zone_heatmap_inputs <- function(
#     rmat,
#     core_obj,
#     zones,
#     cells_meta,
#     markers,
#     zone_levels = c("Non-malignant","Tumor-edge","Mid-tumor","Tumor-core"),
#     max_cells_per_zone = Inf,
#     seed = 1,
#     order_within_zone = c("assignment", "none")
# ) {
#   order_within_zone <- match.arg(order_within_zone)
#
#   # 1) zone per cell_id
#   zone_tbl <- codex_assign_cells_to_tumour_zones(
#     prep_obj = core_obj, zones = zones, zone_levels = zone_levels
#   )
#
#   # 2) add cell_name for rmat
#   id_map <- codex_core_cellid_to_cellname(core_obj, cells_meta)
#
#   zone_tbl <- zone_tbl |>
#     dplyr::left_join(id_map, by = "cell_id") |>
#     dplyr::filter(!is.na(.data$cell_name))
#
#   # 3) keep only cell_names present in rmat
#   zone_tbl <- zone_tbl |>
#     dplyr::filter(.data$cell_name %in% colnames(rmat))
#
#   # 4) optional downsample per zone (display only)
#   if (is.finite(max_cells_per_zone)) {
#     set.seed(seed)
#     zone_tbl <- zone_tbl |>
#       dplyr::group_by(.data$zone) |>
#       dplyr::mutate(.rand = stats::runif(dplyr::n())) |>
#       dplyr::arrange(.data$.rand, .by_group = TRUE) |>
#       dplyr::slice_head(n = max_cells_per_zone) |>
#       dplyr::ungroup() |>
#       dplyr::select(-.data$.rand)
#   }
#
#   # 5) order columns
#   if (order_within_zone == "assignment") {
#     zone_tbl <- zone_tbl |>
#       dplyr::arrange(.data$zone, .data$assignment, .data$cell_name)
#   } else {
#     zone_tbl <- zone_tbl |>
#       dplyr::arrange(.data$zone, .data$cell_name)
#   }
#
#   # 6) subset rmat (markers × selected cells)
#   miss_markers <- setdiff(markers, rownames(rmat))
#   if (length(miss_markers)) stop("Markers not in rmat: ", paste(miss_markers, collapse = ", "))
#
#   cell_names <- zone_tbl$cell_name
#   idx <- match(cell_names, colnames(rmat))
#   keep <- !is.na(idx)
#   cell_names <- cell_names[keep]
#   idx <- idx[keep]
#
#   mat <- rmat[markers, idx, drop = FALSE]
#   colnames(mat) <- cell_names
#
#   # column split (zone facet)
#   zone_split <- zone_tbl$zone
#   names(zone_split) <- zone_tbl$cell_name
#
#   list(
#     mat = mat,
#     zone_tbl = zone_tbl,
#     zone_split = zone_split
#   )
# }
#
#
#
# codex_plot_zone_marker_heatmap <- function(
#     mat,
#     zone_split,
#     zone_tbl,
#     title = NULL,
#     scale_rows = c("none", "z"),
#     show_celltype_bar = TRUE
# ) {
#   scale_rows <- match.arg(scale_rows)
#
#   if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
#       !requireNamespace("circlize", quietly = TRUE)) {
#     stop("Please install ComplexHeatmap + circlize.")
#   }
#
#   M <- mat
#   if (scale_rows == "z") {
#     M <- t(scale(t(M)))
#     M[is.na(M)] <- 0
#   }
#
#   # Color mapping: marker scores are ~0..1, but z-scores can be negative
#   if (scale_rows == "none") {
#     col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("white", "#fddbc7", "#b2182b"))
#   } else {
#     q <- stats::quantile(M, c(0.02, 0.5, 0.98), na.rm = TRUE)
#     col_fun <- circlize::colorRamp2(as.numeric(q), c("#2166ac", "white", "#b2182b"))
#   }
#
#   top_anno <- NULL
#   if (show_celltype_bar) {
#     # a simple annotation bar for your assignment (cell type/state)
#     assign_vec <- zone_tbl$assignment
#     names(assign_vec) <- zone_tbl$cell_name
#
#     top_anno <- ComplexHeatmap::HeatmapAnnotation(
#       Zone = zone_split,
#       CellType = assign_vec,
#       show_annotation_name = TRUE,
#       annotation_name_side = "left"
#     )
#   } else {
#     top_anno <- ComplexHeatmap::HeatmapAnnotation(Zone = zone_split)
#   }
#
#   ComplexHeatmap::Heatmap(
#     M,
#     name = "Score",
#     col = col_fun,
#     column_split = zone_split,
#     cluster_columns = FALSE,
#     cluster_rows = FALSE,
#     show_column_names = FALSE,
#     top_annotation = top_anno,
#     row_names_side = "left",
#     column_title = title,
#     use_raster = TRUE
#   )
# }


# Cellular Density ------------------------------------------------------------

# Global density (ALL cells inside mask)
codex_global_tumour_density_allcells <- function(prep, tumour_poly, pix_per_um, core_id = NA_character_) {
  .codex_use_planar()
  pts <- prep$cent
  inside <- lengths(sf::st_within(pts, sf::st_union(tumour_poly))) > 0
  n_in   <- sum(inside)
  area_um2 <- as.numeric(sum(sf::st_area(tumour_poly)) / (pix_per_um^2))
  tibble::tibble(core = as.character(core_id),
         tumour_area_um2 = area_um2,
         n_cells_inmask = n_in,
         cells_per_mm2 = n_in / (area_um2 / 1e6),
         details = list(tibble::tibble(label = c("inside","outside"),
                               n = c(n_in, nrow(pts)-n_in))))
}

# Local density (count neighbours or occupied area inside a disk; hybrid/fast)
codex_local_density_inmask <- function(prep, tumour_poly,
                                       ref_states = c("pEMT","Epithelial"),
                                       r_um = 50, pix_per_um = 1/0.5067,
                                       mode = c("hybrid","exact"),
                                       min_overlap = 0.75,
                                       exclude_self = TRUE) {
  .codex_use_planar()
  mode <- match.arg(mode)
  r_px <- r_um * pix_per_um
  maskU <- sf::st_union(tumour_poly)
  centroids_sf <- prep$cent
  need <- lengths(sf::st_within(centroids_sf, maskU)) > 0
  centroids_sf <- centroids_sf[need,]
  coords <- sf::st_coordinates(centroids_sf)
  labs   <- centroids_sf$assignment
  idx_ref <- which(labs %in% ref_states)

  if (mode == "hybrid") {
    # neighbour count within radius, constrained to inside-mask
    id <- dbscan::frNN(coords, eps = r_px, sort = FALSE)$id
    out <- lapply(idx_ref, function(i) {
      neigh <- setdiff(id[[i]], if (exclude_self) i else integer(0))
      if (!length(neigh)) return(c(cell_id = centroids_sf$cell_id[i], label = labs[i], n_in_disk = 0, frac_area = NA))
      # drop neighbours whose centroids are outside mask (already filtered)
      n_in <- length(neigh)
      c(cell_id = centroids_sf$cell_id[i], label = labs[i], n_in_disk = n_in, frac_area = NA)
    })
    tibble::as_tibble(do.call(rbind, out)) %>%
      dplyr::mutate(n_in_disk = as.integer(n_in_disk))
  } else {
    # exact: clip disk with mask & cell polygons to compute occupied fraction
    # expensive — only use for spot checks
    stop("mode='exact' not shown here to keep file compact; reuse your validated exact variant when needed.")
  }
}


# Visium-Like Spot Simulation (Hex Grid) --------------------------------------
# Build a "core ROI" polygon from the cells (fast + robust)
# Default: convex hull of all centroids (cheap; good for round-ish TMA cores)
codex_core_roi <- function(prep_obj, method = c("convex_centroids", "bbox"),
                           buffer_px = 0) {
  .codex_use_planar()
  method <- match.arg(method)

  stopifnot(is.list(prep_obj), !is.null(prep_obj$prep), inherits(prep_obj$prep$cent, "sf"))
  pts <- prep_obj$prep$cent

  roi <- switch(
    method,
    convex_centroids = {
      hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(pts)))
      sf::st_sfc(hull, crs = sf::st_crs(pts)) %>% sf::st_sf(geometry = .)
    },
    bbox = {
      bb <- sf::st_bbox(pts)
      p <- sf::st_polygon(list(rbind(
        c(bb["xmin"], bb["ymin"]),
        c(bb["xmax"], bb["ymin"]),
        c(bb["xmax"], bb["ymax"]),
        c(bb["xmin"], bb["ymax"]),
        c(bb["xmin"], bb["ymin"])
      )))
      sf::st_sf(geometry = sf::st_sfc(p, crs = sf::st_crs(pts)))
    }
  )

  if (!is.null(buffer_px) && buffer_px != 0) {
    roi <- roi %>% dplyr::mutate(geometry = sf::st_buffer(.data$geometry, buffer_px)) %>% sf::st_make_valid()
  }
  roi
}

# Generate a hex grid (point centers) covering bbox; then optionally rotate
.codex_hex_grid_points <- function(bbox, dx, dy) {
  xs <- seq(as.numeric(bbox["xmin"]) - dx, as.numeric(bbox["xmax"]) + dx, by = dx)
  ys <- seq(as.numeric(bbox["ymin"]) - dy, as.numeric(bbox["ymax"]) + dy, by = dy)

  out <- lapply(seq_along(ys), function(i) {
    y <- ys[i]
    xoff <- if ((i %% 2) == 0) dx / 2 else 0
    data.frame(x = xs + xoff, y = rep(y, length(xs)))
  })
  dplyr::bind_rows(out)
}

.codex_rotate_xy <- function(x, y, theta_rad, cx, cy) {
  xr <- (x - cx) * cos(theta_rad) - (y - cy) * sin(theta_rad) + cx
  yr <- (x - cx) * sin(theta_rad) + (y - cy) * cos(theta_rad) + cy
  list(x = xr, y = yr)
}

# Core simulation: choose rotation+offset that yields MAX number of spots fully inside ROI
codex_simulate_visium_grid_core <- function(prep_obj,
                                            pix_per_um = 1/0.5067,
                                            spot_diameter_um = 55,
                                            spot_pitch_um    = 100,
                                            roi_method = c("convex_centroids", "bbox"),
                                            roi_buffer_um = 0,
                                            rotations_deg = seq(0, 30, by = 5),
                                            offset_steps = 6,
                                            require_full_spot_inside = TRUE) {
  .codex_use_planar()
  roi_method <- match.arg(roi_method)

  stopifnot(is.list(prep_obj), !is.null(prep_obj$core_id), !is.null(prep_obj$prep$cells_tbl))
  cells_tbl <- prep_obj$prep$cells_tbl

  # units → pixels
  r_px  <- (spot_diameter_um / 2) * pix_per_um
  dx_px <- spot_pitch_um * pix_per_um
  dy_px <- dx_px * sqrt(3) / 2

  roi_buffer_px <- roi_buffer_um * pix_per_um
  roi <- codex_core_roi(prep_obj, method = roi_method, buffer_px = roi_buffer_px)
  roi_u <- sf::st_union(sf::st_geometry(roi))

  # If we require the entire spot circle to fit, shrink ROI by radius
  roi_inner <- if (isTRUE(require_full_spot_inside)) {
    sf::st_buffer(roi_u, dist = -r_px) %>% sf::st_make_valid()
  } else roi_u

  # Centroid of ROI for rotation center
  cc <- sf::st_coordinates(sf::st_centroid(sf::st_sfc(roi_u)))[1, ]
  cx <- as.numeric(cc[1]); cy <- as.numeric(cc[2])

  bb <- sf::st_bbox(roi_u)

  best <- list(n = -Inf, rot = NA_real_, offx = NA_real_, offy = NA_real_, pts = NULL)

  # Local offset search in [0, dx) × [0, dy)
  offx_grid <- seq(0, dx_px, length.out = offset_steps + 1)[- (offset_steps + 1)]
  offy_grid <- seq(0, dy_px, length.out = offset_steps + 1)[- (offset_steps + 1)]

  base_pts <- .codex_hex_grid_points(bb, dx = dx_px, dy = dy_px)

  for (rot in rotations_deg) {
    th <- rot * pi / 180

    for (ox in offx_grid) for (oy in offy_grid) {
      tmp <- base_pts
      tmp$x <- tmp$x + ox
      tmp$y <- tmp$y + oy

      rr <- .codex_rotate_xy(tmp$x, tmp$y, th, cx, cy)
      tmp$x <- rr$x; tmp$y <- rr$y

      pts_sf <- sf::st_as_sf(tmp, coords = c("x", "y"), crs = sf::st_crs(roi))  # POINT sf
      inside <- lengths(sf::st_within(pts_sf, roi_inner)) > 0
      n_in <- sum(inside)

      if (n_in > best$n) {
        best$n <- n_in
        best$rot <- rot
        best$offx <- ox
        best$offy <- oy
        best$pts <- tmp[inside, , drop = FALSE]
      }
    }
  }

  # Assemble spots table
  spots_tbl <- best$pts
  if (is.null(spots_tbl) || !nrow(spots_tbl)) {
    return(list(
      core_id = prep_obj$core_id,
      params = list(r_px = r_px, dx_px = dx_px, dy_px = dy_px, best = best),
      spots_tbl = tibble::tibble(),
      spots_sf  = sf::st_sf()
    ))
  }

  spots_tbl <- spots_tbl %>%
    dplyr::mutate(
      core_id = prep_obj$core_id,
      spot_id = paste0(prep_obj$core_id, "_spot_", dplyr::row_number()),
      radius_px = r_px,
      diameter_px = 2 * r_px,
      diameter_um = spot_diameter_um,
      pitch_um = spot_pitch_um,
      rotation_deg = best$rot,
      offset_x_px = best$offx,
      offset_y_px = best$offy
    ) %>%
    dplyr::relocate("core_id", "spot_id")

  spots_sf <- sf::st_as_sf(spots_tbl, coords = c("x", "y"), crs = sf::st_crs(roi))

  list(
    core_id = prep_obj$core_id,
    params  = list(r_px = r_px, dx_px = dx_px, dy_px = dy_px, roi_method = roi_method,
                   roi_buffer_um = roi_buffer_um, rotations_deg = rotations_deg,
                   offset_steps = offset_steps, require_full_spot_inside = require_full_spot_inside,
                   best = list(n = best$n, rot = best$rot, offx = best$offx, offy = best$offy)),
    roi     = roi,
    spots_tbl = dplyr::as_tibble(spots_tbl),
    spots_sf  = spots_sf
  )
}

# Assign each cell to at most 1 spot: nearest spot, if within radius
codex_assign_cells_to_spots <- function(cells_tbl, spots_tbl) {
  stopifnot(all(c("centroid_x","centroid_y","cell_id","assignment") %in% names(cells_tbl)))
  stopifnot(all(c("spot_id","x","y","radius_px") %in% names(spots_tbl)))

  if (!nrow(spots_tbl)) {
    return(tibble::tibble(cell_id = cells_tbl$cell_id, spot_id = NA_character_, dist_px = NA_real_))
  }

  # Fast nearest-neighbour: each cell → nearest spot
  if (!requireNamespace("RANN", quietly = TRUE)) stop("Please install RANN (install.packages('RANN')).")

  nn <- RANN::nn2(
    data  = as.matrix(spots_tbl[, c("x","y")]),
    query = as.matrix(cells_tbl[, c("centroid_x","centroid_y")]),
    k = 1
  )

  spot_idx <- as.integer(nn$nn.idx[,1])
  dist_px  <- as.numeric(nn$nn.dists[,1])

  # Apply spot radius cutoff (same radius for all spots here; still read from table)
  r_px <- spots_tbl$radius_px[1]
  keep <- dist_px <= r_px

  tibble::tibble(
    cell_id = cells_tbl$cell_id,
    spot_id = ifelse(keep, spots_tbl$spot_id[spot_idx], NA_character_),
    dist_px = dist_px
  )
}

# Summaries per spot: n_cells + fractions per cell type
codex_summarize_spots <- function(cells_tbl, spots_tbl, cell_to_spot, all_types = NULL) {
  stopifnot(all(c("cell_id","assignment") %in% names(cells_tbl)))
  stopifnot(all(c("spot_id") %in% names(spots_tbl)))
  stopifnot(all(c("cell_id","spot_id") %in% names(cell_to_spot)))

  D <- dplyr::left_join(cells_tbl, cell_to_spot, by = "cell_id") %>%
    dplyr::filter(!is.na(.data$spot_id))

  if (is.null(all_types)) all_types <- sort(unique(cells_tbl$assignment))

  counts <- D %>%
    dplyr::count(.data$spot_id, .data$assignment, name = "n_cells_type") %>%
    dplyr::group_by(.data$spot_id) %>%
    dplyr::mutate(
      n_cells_spot = sum(.data$n_cells_type),
      frac_type = .data$n_cells_type / .data$n_cells_spot
    ) %>%
    dplyr::ungroup()

  wide_frac <- counts %>%
    dplyr::select(.data$spot_id, .data$assignment, .data$frac_type) %>%
    tidyr::pivot_wider(
      names_from  = .data$assignment,
      values_from = .data$frac_type,
      values_fill = 0,
      names_prefix = "frac_"
    )

  wide_n <- counts %>%
    dplyr::select(.data$spot_id, .data$assignment, .data$n_cells_type) %>%
    tidyr::pivot_wider(
      names_from  = .data$assignment,
      values_from = .data$n_cells_type,
      values_fill = 0,
      names_prefix = "n_"
    )

  spot_totals <- counts %>%
    dplyr::distinct(.data$spot_id, .data$n_cells_spot)

  spots_tbl %>%
    dplyr::left_join(spot_totals, by = "spot_id") %>%
    dplyr::mutate(n_cells_spot = dplyr::coalesce(.data$n_cells_spot, 0L)) %>%
    dplyr::left_join(wide_n, by = "spot_id") %>%
    dplyr::left_join(wide_frac, by = "spot_id")
}

# Export a QuPath-friendly CSV: centers + diameter in pixels + optional measurements
# IMPORTANT: if prep_obj$flip_used==TRUE, we "unflip" y back into QuPath coordinates using flip_y_to_ref()
codex_export_visium_grid_for_qupath <- function(prep_obj, spots_tbl, out_csv,
                                                include_measurements = TRUE) {
  .codex_use_planar()
  stopifnot(is.list(prep_obj), !is.null(prep_obj$prep$cells_sf))
  stopifnot(all(c("spot_id","x","y","diameter_px","core_id") %in% names(spots_tbl)))

  pts_sf <- sf::st_as_sf(spots_tbl, coords = c("x","y"), crs = sf::st_crs(prep_obj$prep$cells_sf))
  coords <- sf::st_coordinates(pts_sf)

  out <- spots_tbl %>%
    dplyr::mutate(x_sf = coords[,1], y_sf = coords[,2])

  # Convert to QuPath coordinate frame if needed
  if (isTRUE(prep_obj$flip_used)) {
    pts_unflip <- flip_y_to_ref(pts_sf, prep_obj$prep$cells_sf)
    cc <- sf::st_coordinates(pts_unflip)
    out <- out %>%
      dplyr::mutate(x_qp = cc[,1], y_qp = cc[,2], y_was_flipped = TRUE)
  } else {
    out <- out %>%
      dplyr::mutate(x_qp = .data$x_sf, y_qp = .data$y_sf, y_was_flipped = FALSE)
  }

  # keep only useful columns for QuPath
  keep <- c("core_id","spot_id","x_qp","y_qp","diameter_px","diameter_um",
            "n_cells_spot","rotation_deg","offset_x_px","offset_y_px","y_was_flipped")
  keep <- intersect(keep, names(out))

  if (!include_measurements) {
    out2 <- out[, keep, drop = FALSE]
  } else {
    # also keep n_* and frac_* columns if they exist
    extra <- grep("^(n_|frac_)", names(out), value = TRUE)
    out2 <- out[, c(keep, extra), drop = FALSE]
  }

  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(out2, out_csv, row.names = FALSE, quote = FALSE)
  invisible(out2)
}


extract_top_mps_from_visium_sim <- function(rds_file) {
  obj   <- readRDS(rds_file)
  spots <- obj$spots

  # 1) Extract fraction columns
  frac_cols <- grep("^frac_", colnames(spots), value = TRUE)
  if (length(frac_cols) == 0) {
    warning("No frac_ columns found in: ", basename(rds_file))
    return(NULL)
  }

  frac_mat <- as.matrix(spots[, frac_cols, drop = FALSE])
  rownames(frac_mat) <- spots$spot_id
  colnames(frac_mat) <- sub("^frac_", "", frac_cols)

  # 2) Drop Undefined + NA->0 + renormalize row-wise
  if ("Undefined" %in% colnames(frac_mat)) {
    frac_mat <- frac_mat[, colnames(frac_mat) != "Undefined", drop = FALSE]
  }
  frac_mat[is.na(frac_mat)] <- 0

  row_sums <- rowSums(frac_mat)
  ok <- row_sums > 0  # spots that have any non-Undefined mass

  if (!any(ok)) {
    warning("All spots are empty after dropping Undefined in: ", basename(rds_file))
    return(NULL)
  }

  frac_mat <- frac_mat[ok, , drop = FALSE]
  row_sums <- row_sums[ok]
  frac_mat <- sweep(frac_mat, 1, row_sums, "/")  # renormalize to sum=1

  # 3) Top-1 and top-2 per row (vectorized)
  n <- nrow(frac_mat)
  if (ncol(frac_mat) == 1) {
    idx1 <- rep(1L, n)
    val1 <- as.numeric(frac_mat[, 1])
    idx2 <- rep(1L, n)
    val2 <- rep(0, n)
    nm1  <- colnames(frac_mat)[1]
    nm2  <- NA_character_
  } else {
    idx1 <- max.col(frac_mat, ties.method = "first")
    val1 <- frac_mat[cbind(seq_len(n), idx1)]
    nm1  <- colnames(frac_mat)[idx1]

    frac2 <- frac_mat
    frac2[cbind(seq_len(n), idx1)] <- -Inf

    idx2 <- max.col(frac2, ties.method = "first")
    val2 <- frac2[cbind(seq_len(n), idx2)]
    val2[is.infinite(val2)] <- 0
    nm2  <- colnames(frac_mat)[idx2]
    nm2[val2 == 0] <- NA_character_
  }

  out <- data.frame(
    core_id      = obj$core_id %||% NA_character_,
    spot_id      = rownames(frac_mat),
    First_MP     = nm1,
    First_Score  = as.numeric(val1),
    Second_MP    = nm2,
    Second_Score = as.numeric(val2),
    Top_Two      = as.numeric(val1 + val2),
    stringsAsFactors = FALSE
  )
  out
}


# Plotting --------------------------------------------------------------------

codex_plot_core <- function(cells_sf_core,
                            boundaries          = NULL,
                            cell_palette        = NULL,
                            legend_boundaries   = FALSE,
                            boundary_colours    = NULL,
                            boundary_linewidths = NULL,
                            boundary_linetypes  = NULL,
                            cell_stroke_size    = 0.02,
                            title               = NULL,
                            legend_title        = "Cell type") {
  .codex_use_planar()
  stopifnot(inherits(cells_sf_core, "sf"),
            "assignment" %in% names(cells_sf_core))
  if (is.null(cell_palette)) {
    levs <- sort(unique(cells_sf_core$assignment))
    cell_palette <- setNames(viridis::viridis(length(levs)), levs)
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = cells_sf_core, ggplot2::aes(fill = .data$assignment),
            colour = "black", linewidth = cell_stroke_size) +
    ggplot2::scale_fill_manual(values = cell_palette, name = legend_title)

  if (!is.null(boundaries)) {
    if (inherits(boundaries, "sf")) boundaries <- list(boundary = boundaries)
    else if (!is.list(boundaries)) stop("`boundaries` must be sf or list of sf.")
    n_b <- length(boundaries)
    if (is.null(boundary_colours))
      boundary_colours <- RColorBrewer::brewer.pal(max(3, n_b), "Set1")[seq_len(n_b)]
    if (is.null(boundary_linewidths))
      boundary_linewidths <- rep(0.6, n_b)
    if (is.null(boundary_linetypes))
      boundary_linetypes <- rep("solid", n_b)
    if (length(boundary_linetypes) == 1L)
      boundary_linetypes <- rep(boundary_linetypes, n_b)
    nm <- names(boundaries) %||% paste0("b", seq_len(n_b))
    for (i in seq_len(n_b)) {
      bsf <- boundaries[[i]]
      if (!inherits(bsf, "sf")) stop(sprintf("Element %s is not an sf object.", nm[i]))
      if (is.na(sf::st_crs(bsf))) bsf <- sf::st_set_crs(bsf, NA)
      p <- p + ggplot2::geom_sf(data = bsf, fill = NA,
                       colour = boundary_colours[[i]],
                       linewidth = boundary_linewidths[[i]],
                       linetype = boundary_linetypes[[i]],
                       inherit.aes = FALSE, show.legend = legend_boundaries)
    }
  }
  p + ggplot2::coord_sf(datum = NA, expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background  = ggplot2::element_rect(fill = "black", colour = NA),
          panel.background = ggplot2::element_rect(fill = "black", colour = NA),
          legend.text      = ggplot2::element_text(colour = "white", size = 14),
          legend.title     = ggplot2::element_text(colour = "white", size = 16),
          legend.position  = "right") +
    { if (!is.null(title)) ggplot2::labs(title = title) else list() }
}


common_theme <- ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.text  = ggplot2::element_text(colour = "black"),
    axis.title = ggplot2::element_text(colour = "black"),
    axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.9),
    axis.line  = ggplot2::element_line(colour = "black", linewidth = 0.9),
    legend.position = "none",
    plot.margin = ggplot2::margin(4, 4, 2, 4, "pt")
  )


codex_make_boxdot <- function(
    dat, yvar, ylab, binwidth,
    box_width = 0.6,
    dot_mm = 2,
    x_var = NULL,
    x_labels = NULL,                 # optional named vector for relabels
    x_order  = NULL,                 # optional character vector for order
    tilt_x_text = FALSE,
    y_limits = NULL,
    panel_w = grid::unit(30, "mm"),
    panel_h = grid::unit(60, "mm"),
    expend_mult_y = 0.005,
    expand_top = NULL,
    wrap_y_title = TRUE,
    ylab_width   = NULL,
    percent_y = FALSE,
    percent_accuracy = 1,
    dot_fill_fixed = "#0c2a50ff",
    dot_outline    = "#E7E9EB",
    dot_alpha      = 1,
    box_palette    = c(pEMT_core = "#B24745FF", pEMT_edge = "#DF8F44FF"),
    test_method      = c("wilcox.test", "t.test"),
    test_method_args = list()   # e.g., list(exact=TRUE) or list(var.equal=FALSE, alternative="greater")
) {

  test_method <- match.arg(test_method)

  # pick x variable if not supplied
  if (is.null(x_var)) {
    x_var <- intersect(c("Set","pEMT_Pattern","coherent","group"), names(dat))[1]
    if (length(x_var) == 0) stop("Provide x_var or include a column named Set, pEMT_Pattern, coherent, or group.")
  }
  if (!x_var %in% names(dat)) stop("x_var='", x_var, "' not found in data.")

  # map x values to nice labels
  dat2 <- dplyr::mutate(
    dat,
    .x_raw = as.character(.data[[x_var]]),
    x_plot = dplyr::case_when(
      stringr::str_detect(.x_raw, stringr::regex("core", ignore_case = TRUE)) ~ "Core",
      stringr::str_detect(.x_raw, stringr::regex("edge|margin", ignore_case = TRUE)) ~ "Edge",
      .x_raw %in% c("TRUE","T","1")  ~ "Coherent",
      .x_raw %in% c("FALSE","F","0") ~ "Non-coherent",
      TRUE ~ .x_raw
    )
  )

  if (!is.null(x_labels))  dat2$x_plot <- dplyr::recode(dat2$x_plot, !!!x_labels)
  if (!is.null(x_order))   dat2$x_plot <- factor(dat2$x_plot, levels = x_order)
  else if (all(c("Non-coherent","Coherent") %in% unique(dat2$x_plot)))
    dat2$x_plot <- factor(dat2$x_plot, levels = c("Non-coherent","Coherent"))
  else if (all(c("Core","Edge") %in% unique(dat2$x_plot)))
    dat2$x_plot <- factor(dat2$x_plot, levels = c("Core","Edge"))
  else dat2$x_plot <- factor(dat2$x_plot)

  # clean NAs in y
  dat2 <- dat2[!is.na(dat2[[yvar]]), ]

  # y label wrap
  ylab_final <- ylab %||% yvar
  if (isTRUE(wrap_y_title) && is.character(ylab_final)) {
    w <- ylab_width
    if (is.null(w)) w <- if (nchar(ylab_final) > 28) 28 else Inf
    if (is.finite(w)) ylab_final <- stringr::str_wrap(ylab_final, width = w)
  }

  # constant dot size (mm) regardless of binwidth
  panel_h_mm <- as.numeric(grid::convertUnit(panel_h, "mm", valueOnly = TRUE))
  yrange <- if (is.null(y_limits)) diff(range(dat2[[yvar]], na.rm = TRUE)) else diff(y_limits)
  if (yrange <= 0) yrange <- 1
  dotsize_eff <- (dot_mm / panel_h_mm) * (yrange / binwidth)

  # stats selection
  lv <- levels(dat2$x_plot)
  n_groups <- length(lv)
  if (n_groups == 2) {
    comps <- list(lv)
    use_ggpubr <- TRUE
  } else if (n_groups > 2) {
    use_ggpubr <- FALSE
    kw_pval <- stats::kruskal.test(dat2[[yvar]] ~ dat2$x_plot)$p.value
    kw_label <- if (kw_pval < 1e-3) "p < 0.001" else paste0("p = ", formatC(kw_pval, format="fg", digits=2))
  } else use_ggpubr <- FALSE

  # ---- plot ----
  p <- ggplot2::ggplot(dat2, ggplot2::aes(x = x_plot, y = .data[[yvar]])) +
    # boxes colored by group
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = x_plot),
      linewidth = 0.9, width = box_width, fatten = 1,
      outlier.shape = NA
    ) +
    # dots with fixed fill & soft outline
    ggplot2::geom_dotplot(
      binaxis = "y", method = "dotdensity",
      stackdir = "center", stackgroups = TRUE, binpositions = "all",
      binwidth = binwidth, dotsize = dotsize_eff,
      fill = dot_fill_fixed, colour = dot_outline, alpha = dot_alpha, stroke = 0.9
    ) +
    ggplot2::scale_fill_manual(values = box_palette, guide = "none")

  # y scale & labels
  top_expansion <- if (is.null(expand_top)) if (n_groups > 1) 0.12 else 0.06 else expand_top
  p <- p +
    ggplot2::labs(x = NULL, y = ylab_final) +
    ggplot2::scale_y_continuous(
      limits = y_limits,
      expand = ggplot2::expansion(mult = c(expend_mult_y, top_expansion)),
      labels = if (percent_y) scales::percent_format(accuracy = percent_accuracy) else ggplot2::waiver()
    )

  # stats label
  if (use_ggpubr && n_groups == 2) {
    p <- p + ggpubr::stat_compare_means(
      ggplot2::aes(group = x_plot),
      method      = test_method,
      method.args = test_method_args,
      comparisons = comps,
      label.y = if (is.null(y_limits)) max(dat2[[yvar]], na.rm = TRUE) else y_limits[2],
      size = 4
    )
  } else if (n_groups > 2) {
    y_pos <- if (is.null(y_limits)) max(dat2[[yvar]], na.rm = TRUE) else y_limits[2]
    x_pos <- (n_groups + 1) / 2
    p <- p + ggplot2::annotate("text", x = x_pos, y = y_pos, label = kw_label,
                               size = 4, hjust = 0.5, vjust = -0.5)
  }

  p +
    (if (exists("common_theme")) common_theme else ggplot2::theme_classic()) +
    ggplot2::coord_cartesian(clip = "off") +
    ggh4x::force_panelsizes(rows = panel_h, cols = panel_w, respect = TRUE) +
    (if (tilt_x_text) ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) else ggplot2::theme())
}


plot_proximity_heatmap <- function(res,
                                   order = NULL,
                                   value_col = "z_score",
                                   fdr_thresh = 0.05,
                                   show_signif_marks = TRUE,
                                   agg_value_fun = function(x) mean(x, na.rm = TRUE),
                                   agg_padj_fun  = function(x) min(x, na.rm = TRUE),
                                   drop_self = TRUE,
                                   title = NULL) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
  })

  # Accept either the saved list (with $summary) or a tibble/data.frame
  df <- res
  if (is.list(res) && !is.data.frame(res) && "summary" %in% names(res)) {
    df <- res$summary
  }
  if (!is.data.frame(df)) stop("Input must be a data.frame/tibble or a list with $summary.")

  required <- c("from_cell_type_figure", "to_cell_type_figure")
  miss_req <- setdiff(required, colnames(df))
  if (length(miss_req) > 0) stop("Missing required columns: ", paste(miss_req, collapse = ", "))

  if (!(value_col %in% colnames(df))) {
    stop("value_col='", value_col, "' not found. Available cols: ",
         paste(colnames(df), collapse = ", "))
  }

  df <- df %>%
    dplyr::mutate(
      from_cell_type_figure = as.character(from_cell_type_figure),
      to_cell_type_figure   = as.character(to_cell_type_figure)
    )

  # Ensure padj exists (compute if missing)
  if (!("padj" %in% colnames(df))) {
    if ("p_value" %in% colnames(df)) {
      df <- df %>% dplyr::mutate(padj = stats::p.adjust(p_value, method = "BH"))
    } else {
      # If neither exists, create NA padj (significance marks will be disabled)
      df <- df %>% dplyr::mutate(padj = NA_real_)
    }
  }

  # Optional: drop self interactions
  if (drop_self) {
    df <- df %>% dplyr::filter(from_cell_type_figure != to_cell_type_figure)
  }

  # Aggregate duplicates per (from,to)
  df_pair <- df %>%
    dplyr::group_by(from_cell_type_figure, to_cell_type_figure) %>%
    dplyr::summarise(
      value = agg_value_fun(.data[[value_col]]),
      padj  = agg_padj_fun(padj),
      .groups = "drop"
    )

  # Determine axis order
  types <- sort(unique(c(df_pair$from_cell_type_figure, df_pair$to_cell_type_figure)))

  if (!is.null(order)) {
    if (any(duplicated(order))) {
      warning("Your `order` vector contains duplicates. Using unique(order).")
    }
    ord <- unique(as.character(order))
    ord <- ord[ord %in% types]
    ord <- c(ord, setdiff(types, ord))  # append anything missing
  } else {
    ord <- types
  }

  # Build full grid to include missing pairs as NA
  grid <- tidyr::expand_grid(
    from_cell_type_figure = ord,
    to_cell_type_figure   = ord
  ) %>%
    dplyr::left_join(df_pair, by = c("from_cell_type_figure", "to_cell_type_figure")) %>%
    dplyr::mutate(
      from_cell_type_figure = factor(from_cell_type_figure, levels = ord),
      to_cell_type_figure   = factor(to_cell_type_figure,   levels = ord)
    )

  # Base plot
  p <- ggplot2::ggplot(grid, ggplot2::aes(x = to_cell_type_figure, y = from_cell_type_figure, fill = value)) +
    ggplot2::geom_tile(color = NA) +
    ggplot2::scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      na.value = "grey90",
      name = value_col
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid  = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL, fill = value_col, title = title)

  # Add significance marks if we can
  if (show_signif_marks) {
    if (!all(is.na(grid$padj))) {
      p <- p + ggplot2::geom_point(
        data = grid %>% dplyr::filter(!is.na(padj) & padj < fdr_thresh),
        ggplot2::aes(x = to_cell_type_figure, y = from_cell_type_figure),
        inherit.aes = FALSE,
        shape = 8,
        size = 2
      )
    } else {
      warning("padj is NA for all rows (no p_value/padj in input). Significance marks disabled.")
    }
  }

  p
}
