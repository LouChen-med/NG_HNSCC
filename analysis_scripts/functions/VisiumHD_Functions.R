# Load Spatial Image ------------------------------------------------------

# Function to load the spatial image with scale factors and coordinates
load_spatial_image <- function(path_to_samp_dir, image_res = "lowres", filter = TRUE, header = FALSE) {
  # Initialize file paths
  spatial_dir <- file.path(path_to_samp_dir, "spatial")
  if (dir.exists(spatial_dir)) {
    # Look in the "spatial" subdirectory
    spatial_files <- list.files(spatial_dir, full.names = TRUE)
  } else {
    # Look in the root of the sample directory
    message("     === WARNING: Spatial directory not found. Searching in the sample directory ===")
    spatial_files <- list.files(path_to_samp_dir, full.names = TRUE)
  }

  # Match files using regex
  tissue_positions_file <- grep("*tissue_positions*", spatial_files, value = TRUE)
  scalefactors_file <- grep("scalefactors_json\\.json$", spatial_files, value = TRUE)
  image_file <- grep("image\\.png$", spatial_files, value = TRUE)

  # Check for missing files
  missing_files <- c()
  if (length(tissue_positions_file) == 0) missing_files <- c(missing_files, "tissue_positions_file")
  if (length(scalefactors_file) == 0) missing_files <- c(missing_files, "scalefactors_file")
  if (length(image_file) == 0) missing_files <- c(missing_files, "image_file")

  # Handle missing files
  if (length(missing_files) > 0) {
    stop("     === ERROR: Missing required spatial data files: ", paste(missing_files, collapse = ", "), " ===")
  }

  # Prioritize the desired resolution for the image file
  if (image_res == "lowres") {
    if (length(grep("lowres", image_file)) > 0) {
      image_file <- grep("lowres", image_file, value = TRUE)
    } else if (length(grep("hires", image_file)) > 0) {
      image_file <- grep("hires", image_file, value = TRUE)
      message("     === WARNING: Low-resolution image not found. Using high-resolution image instead ===")
    }
  } else if (image_res == "hires") {
    if (length(grep("hires", image_file)) > 0) {
      image_file <- grep("hires", image_file, value = TRUE)
    } else if (length(grep("lowres", image_file)) > 0) {
      image_file <- grep("lowres", image_file, value = TRUE)
      message("     === WARNING: High-resolution image not found. Using low-resolution image instead ===")
    }
  }

  # Load data
  image <- png::readPNG(image_file[1])

  scale_factors <- jsonlite::fromJSON(scalefactors_file[1])
  if (any(grepl("bin_size_um", names(scale_factors)))) {
    assay <- "VisiumHD"
  } else {
    assay <- "Visium"
  }
  scale_factors <- scale_factors[c("spot_diameter_fullres", "tissue_hires_scalef", "fiducial_diameter_fullres", "tissue_lowres_scalef")]
  names(scale_factors) <- c("spot", "hires", "fiducial", "lowres")

  if (tools::file_ext(tissue_positions_file[1]) == "parquet") {
    # Check if the `arrow` package is available, this package must be installed to read parquet files
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Please install arrow to read parquet files")
    }
    coordinates <- as.data.frame(arrow::read_parquet(tissue_positions_file[1]))
    colnames(coordinates) <- c("SpotID", "in_tissue", "array_row", "array_col", "pxl_in_rows", "pxl_in_cols")
  } else {
    coordinates <- read.csv(tissue_positions_file[1], header = header)
    if (!header && all(sapply(coordinates, class) == "character")) {
      coordinates <- read.csv(tissue_positions_file[1], header = TRUE)
    }
    colnames(coordinates) <- c("SpotID", "in_tissue", "array_row", "array_col", "pxl_in_rows", "pxl_in_cols")
  }

  # Optionally filter to in-tissue spots
  if (filter) {
    coordinates <- coordinates[which(coordinates$in_tissue == 1), , drop = FALSE]
  }

  # Calculate spot radius
  if (assay == "VisiumHD") {
    radius <- scale_factors$spot * scale_factors[[image_res]]
  } else {
    radius <- scale_factors$fiducial * scale_factors[[image_res]]
  }
  spot_radius <- radius / max(dim(image))

  return(list(Image = image, Scale = scale_factors, Coordinates = coordinates, Spot_Radius = spot_radius))
}



# Spatial Scale Bar --------------------------------------------------------

# Add scale bar to the spatial plot
add_spatial_scale_bar <- function(
    p,
    image_obj,
    microns_per_pixel,
    bar_um = 1000,
    label = "1 mm",
    x_offset = 0.06,
    y_offset = 0.06
) {
  coords <- .get_coords(image_obj, scale = "lowres")

  x_vals <- coords[, "pxl_in_cols"]
  y_vals <- coords[, "pxl_in_rows"]

  x_range <- range(x_vals, na.rm = TRUE)
  y_range <- range(y_vals, na.rm = TRUE)

  # Length in the same low-resolution coordinate units used by the plot
  bar_length <- bar_um * image_obj$Scale[["lowres"]] / microns_per_pixel

  x_start <- x_range[1] + x_offset * diff(x_range)
  y_start <- y_range[1] + y_offset * diff(y_range)

  p +
    annotate(
      "segment",
      x = x_start,
      xend = x_start + bar_length,
      y = y_start,
      yend = y_start,
      linewidth = 1.2,
      colour = "black"
    ) +
    annotate(
      "text",
      x = x_start + bar_length / 2,
      y = y_start + 0.025 * diff(y_range),
      label = label,
      size = 3.5,
      vjust = 0,
      colour = "black"
    )
}



# Spatial Custom Geom -------------------------------------------------------

# Spatial custom geom
geom_spatial <- function(mapping = NULL,
                         data = NULL,
                         image_obj = image_obj,
                         stat = "identity",
                         position = "identity",
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = TRUE,
                         image.alpha = image.alpha,
                         crop = crop,
                         ...) {

  GeomCustom <- ggproto("GeomCustom",
                        Geom,
                        required_aes = c("x", "y"),
                        extra_params = c("na.rm", "image_obj", "image.alpha", "crop"),
                        default_aes = aes(
                          shape = 21,
                          colour = "black",
                          point.size.factor = 1.0,
                          fill = NA,
                          alpha = NA,
                          stroke = NA
                        ),

                        setup_data = function(self, data, params) {
                          data <- ggproto_parent(Geom, self)$setup_data(data, params)
                          # We need to flip the image as the Y coordinates are reversed
                          data$y = max(data$y) - data$y + min(data$y)
                          data
                        },

                        draw_key = draw_key_point,
                        draw_group = function(data, panel_scales, coord, image_obj, image.alpha, crop) {
                          if (!crop) {
                            y.transform <- c(0, nrow(image_obj$Image)) - panel_scales$y.range
                            data$y <- data$y + sum(y.transform)
                            panel_scales$x$continuous_range <- c(0, ncol(image_obj$Image))
                            panel_scales$y$continuous_range <- c(0, nrow(image_obj$Image))
                            panel_scales$y.range <- c(0, nrow(image_obj$Image))
                            panel_scales$x.range <- c(0, ncol(image_obj$Image))
                          }

                          z <- coord$transform(
                            data.frame(x = c(0, ncol(image_obj$Image)), y = c(0, nrow(image_obj$Image))),
                            panel_scales)
                          z$y <- -rev(z$y) + 1
                          wdth <- z$x[2] - z$x[1]
                          hgth <- z$y[2] - z$y[1]
                          vp <- grid::viewport(x = grid::unit(z$x[1], units = "npc"),
                                               y = grid::unit(z$y[1], units = "npc"),
                                               width = grid::unit(wdth, units = "npc"),
                                               height = grid::unit(hgth, units = "npc"),
                                               just = c("left", "bottom"))

                          i_grob <- grid::rasterGrob(image_obj$Image, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
                          img_grob <- grid::editGrob(i_grob, vp = vp)
                          spot.size <- image_obj$Spot_Radius
                          coords <- coord$transform(data, panel_scales)
                          pts <- grid::pointsGrob(x = coords$x,
                                                  y = coords$y,
                                                  pch = data$shape,
                                                  size = grid::unit(spot.size, "npc") * data$point.size.factor,
                                                  gp = gpar(col = alpha(colour = coords$colour, alpha = coords$alpha),
                                                            fill = alpha(colour = coords$fill, alpha = coords$alpha),
                                                            lwd = coords$stroke)
                          )

                          vp <- grid::viewport()
                          gt <- grid::gTree(vp = vp)
                          if (image.alpha > 0) {
                            if (image.alpha != 1) {
                              img_grob$raster = as.raster(matrix(
                                data = alpha(colour = img_grob$raster, alpha = image.alpha),
                                nrow = nrow(img_grob$raster),
                                ncol = ncol(img_grob$raster),
                                byrow = TRUE)
                              )
                            }
                            gt <- grid::addGrob(gTree = gt, child = img_grob)
                          }
                          gt <- grid::addGrob(gTree = gt, child = pts)
                          # Replacement for ggname
                          gt$name <- grobName(grob = gt, prefix = 'geom_spatial')
                          return(gt)
                          # ggplot2:::ggname("geom_spatial", gt)
                        }
  )

  layer(geom = GeomCustom,
        mapping = mapping,
        data = data,
        stat = stat,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, image_obj = image_obj, image.alpha = image.alpha, crop = crop, ...)
  )
}



# Coordinate Extraction Helper ---------------------------------------------

# Helper function to extract the coordinate from the loaded image
.get_coords <- function(image_obj, scale = "lowres") {

  if (!is.null(scale)) {
    coords <- image_obj$Coordinates[, c("pxl_in_rows", "pxl_in_cols")]
    rownames(coords) <- image_obj$Coordinates$SpotID
    scale <- match.arg(arg = scale, choices = c("spot", "hires", "fiducial", "lowres"))
    scale_use <- image_obj$Scale[[scale]]
    coords <- coords * scale_use
  } else {
    coords <- image_obj$Coordinates[, c("pxl_in_rows", "pxl_in_cols")]
  }
  return(coords)
}



# Plot Spatial Features ------------------------------------------------------

# Plot spatial features (metaprogram-IDs, tissue zones, etc')
# Numeric (such as, expression values) can also be ploted as a contineous variable bu must be first added to the metadata
plot_spatial_features <- function(metadata,
                                  image_obj,
                                  color_by,
                                  barcodes = "SpotID",
                                  cols = NULL,
                                  magma = FALSE,
                                  magma_rev = FALSE,
                                  limits = NULL,
                                  midpoint = NULL,
                                  image.alpha = 0.6,
                                  crop = TRUE,
                                  scale = "lowres",
                                  pt.size = 1.6,
                                  pt.alpha = 1,
                                  stroke = 0,
                                  na.value = 'transparent',
                                  rm_spots = metadata[[barcodes]][is.na(metadata[[color_by]])]) {

  # Spatial image object most be filtered to contain only existing spots
  if (any(image_obj$Coordinates$in_tissue == 0)) {
    stop("The image object provided most be filtered! use - load_spatial_image(path_to_image_dir, filter = TRUE)")
  }

  data <- metadata[, c(barcodes, color_by)]
  if (is.character(data[[color_by]])) {
    data[, color_by] <- factor(x = data[[color_by]])
  }

  coordinates <- .get_coords(image_obj, scale = scale) %>%
    tibble::rownames_to_column(var = barcodes)
  data <- merge(coordinates, data, by = barcodes, all = TRUE) %>%
    tibble::column_to_rownames(var = barcodes)

  if (is.numeric(data[[color_by]])) {
    plot <- ggplot(data = data, aes_string(x = colnames(data)[2], y = colnames(data)[1], fill = color_by, alpha = color_by))
  } else {
    plot <- ggplot(data = data, aes_string(x = colnames(data)[2], y = colnames(data)[1], fill = color_by))
  }
  plot <- plot + geom_spatial(data = data,
                              point.size.factor = ifelse(rownames(data) %in% rm_spots, yes = 0, no = pt.size),
                              image_obj = image_obj,
                              image.alpha = image.alpha,
                              crop = crop,
                              stroke = ifelse(rownames(data) %in% rm_spots, yes = 0, no = stroke))
  plot <- plot + geom_spatial(data = data,
                              point.size.factor = ifelse(rownames(data) %in% rm_spots, yes = 0, no = pt.size),
                              image_obj = image_obj,
                              image.alpha = 0,
                              crop = crop,
                              stroke = ifelse(rownames(data) %in% rm_spots, yes = 0, no = stroke),
                              alpha = pt.alpha)

  if (!is.numeric(data[[color_by]])) {
    if (is.null(cols)) {
      cols <- gg_color_hue(length(levels(data[, color_by])))
    } else if (!is.null(cols)) {
      if (length(cols) < length(levels(data[, color_by]))) {
        warning("Not enough colors provided! seting default colors...")
        cols <- gg_color_hue(length(levels(data[, color_by])))
      } else {
        cols <- cols[1:length(levels(data[, color_by]))]
      }
    }
    scale <- scale_fill_manual(name = , values = cols, na.value = na.value, na.translate = FALSE)
    plot <- plot + scale
  }

  if (is.numeric(data[[color_by]])) {
    if (purrr::is_null(limits)) {
      limits <- range(data[[color_by]], na.rm = TRUE)
    }
    if (purrr::is_null(midpoint)) {
      midpoint <- round(mean(limits), digits = 1)
    }
    if (magma) {
      if (magma_rev) {
        # c_cols <- rev(c(colorRampPalette(c("white", rev(viridis::magma(323, begin = 0.15))[1]))(10), rev(viridis::magma(323, begin = 0.18))))
        c_cols <- rev(viridis::magma(n = 6))
      } else {
        # c_cols <- c(colorRampPalette(c("white", rev(viridis::magma(323, begin = 0.15))[1]))(10), rev(viridis::magma(323, begin = 0.18)))
        c_cols <- viridis::magma(n = 6)
      }
    } else if (isTRUE(min(na.omit(data[[color_by]])) < 0)) {
      c_cols <- grDevices::colorRampPalette(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100)
    } else if (isTRUE(min(na.omit(data[[color_by]])) >= 0)) {
      c_cols <- grDevices::colorRampPalette(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")[1:6]))(50)
    }

    if (length(midpoint) == 1 && midpoint == 0) {
      grad_steps <- (max(limits) - min(limits)) / (length(c_cols) - 1)
      col_vals <- seq(from = min(limits), to = max(limits), by = grad_steps)
    } else if (length(midpoint) == 1 && midpoint != 0) {
      zero_col_val <- round(length(c_cols) / 2)
      low_grad_steps <- (midpoint - min(limits)) / (zero_col_val - 1)
      high_grad_steps <- (max(limits) - midpoint) / (zero_col_val - 1)
      col_vals <- c(seq(from = min(limits), to = midpoint, by = low_grad_steps), seq(from = midpoint, to = max(limits), by = high_grad_steps)[-1])
    } else if (length(midpoint) > 1) {
      zero_col_val <- round(length(c_cols) / 2)
      c_cols <- c(c_cols[1:(zero_col_val - 1)], rep(c_cols[zero_col_val], 2), c_cols[(zero_col_val + 1):length(c_cols)])
      low_grad_steps <- (min(midpoint) - min(limits)) / (zero_col_val - 1)
      high_grad_steps <- (max(limits) - max(midpoint)) / (zero_col_val - 1)
      col_vals <- c(seq(from = min(limits), to = min(midpoint), by = low_grad_steps),
                    seq(from = max(midpoint), to = max(limits), by = high_grad_steps))
    }

    scale <- scale_fill_gradientn(colors = c_cols, values = scales::rescale(col_vals), limits = limits, oob = scales::squish, na.value = na.value)
    plot <- plot + eval(scale) + scale_alpha_continuous(range = c(0, 0.1)) + guides(alpha = "none")
  }

  plot <- plot + coord_fixed() + theme(aspect.ratio = 1) + NoAxes() + theme(panel.background = element_blank())
  return(plot)
}



# Metaprogram Scoring (Top-2 Assignment) ------------------------------------

metaprog_score2 <- function(
    expr_mat,
    metaprog_genes,
    score_diff      = 0,
    score_cut       = NULL,
    bin_size        = 100,
    bin_num         = 30,
    center_rows     = TRUE,
    conserved.genes = 0.7,
    scale           = TRUE
) {
  ## 1. Sanity check
  if (!is.list(metaprog_genes) || is.null(names(metaprog_genes))) {
    stop("`metaprog_genes` must be a named list of character vectors.")
  }

  ## 2. Build a compact two-column map: Program - Gene
  metaprog_df <- stack(metaprog_genes)
  names(metaprog_df) <- c("Gene", "Program")
  programs <- unique(metaprog_df$Program)

  ## 3. Prepare three streaming vectors of length = n_cells
  cell_ids     <- colnames(expr_mat)
  n_cells      <- length(cell_ids)
  first_score  <- rep(-Inf, n_cells)
  second_score <- rep(-Inf, n_cells)
  prog_label   <- rep(NA_character_, n_cells)
  names(first_score)  <- names(second_score) <- names(prog_label) <- cell_ids

  ## 4. Loop over each program, scoring and updating top-2 in one pass
  for (prog in programs) {
    genes <- metaprog_df$Gene[metaprog_df$Program == prog]
    sc    <- score(
      expr_mat,
      genes,
      bin_size        = bin_size,
      bin_num         = bin_num,
      center_rows     = center_rows,
      conserved.genes = conserved.genes
    )
    # Ensure same ordering
    sc <- sc[cell_ids]

    # Which cells beat the current first_score?
    better <- sc > first_score
    # Update second_score for those beaten by sc
    second_score[better] <- first_score[better]
    # Update first_score & label
    first_score[better]  <- sc[better]
    prog_label[better]   <- prog

    # Now see if sc slots into 2nd place for the rest
    runner_up <- sc > second_score & !better
    second_score[runner_up] <- sc[runner_up]
  }

  ## 5. Apply score_diff / score_cut logic
  if (score_diff != 0) {
    # Per-cell difference threshold
    diff_thr <- if (scale) score_diff * first_score else rep(score_diff, n_cells)
    ok_diff  <- (first_score - second_score) > diff_thr
    ok_cut   <- is.null(score_cut) | (first_score > score_cut)
    final_label <- ifelse(ok_diff & ok_cut, prog_label, "Unresolved")
  } else {
    if (is.null(score_cut)) {
      final_label <- prog_label
    } else {
      final_label <- ifelse(first_score > score_cut, prog_label, "Unresolved")
    }
  }

  ## 6. Return only the small vectors
  list(
    meta_program = final_label,
    # If you really need the raw first/second scores, you can return them too
    first_score  = first_score,
    second_score = second_score
  )
}
