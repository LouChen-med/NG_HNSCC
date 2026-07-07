library(FNN)
library(kohonen)

# Clustering and Subclustering --------------------------------------------

make_hclust <- function(method = "average",
                        cor_matrix = marker_cor) {
  m_dist <- as.dist(1 - (cor_matrix))
  m_hc  <- hclust(m_dist, method = method)
  m_hc  <- reorder(as.dendrogram(m_hc), colMeans(cor_matrix))
}



Run_subclust <- function(matrix = amat,
                         cluster = "gating",
                         markers = c("CD3", "CD4", "CD8", "CD15"),
                         cell_type_state = "vasc",
                         x_grid = 5,
                         y_grid = 5,
                         color_min = -2,
                         color_max = 2) {

  # Filter cells based on the specified cluster and cell_type_state
  selected_cells <- cells %>%
    dplyr::filter(!!rlang::sym(cluster) %in% c(cell_type_state)) %>%
    dplyr::pull(cell_name)

  print(paste("subsetting", length(selected_cells), "cells of type", cell_type_state, sep = " "))

  tmp_mat <- matrix[markers, selected_cells] %>%
    as.matrix() %>%
    t()

  print(paste("dimension of matrix is", dim(tmp_mat)[1], "x", dim(tmp_mat)[2], sep = " "))

  # Initialize SOM with specified grid dimensions
  set.seed(128)  # Set the seed
  som_model <- kohonen::som(
    tmp_mat,
    grid = kohonen::somgrid(x_grid, y_grid, "rectangular"),
    rlen = 10,
    alpha = c(0.05, 0.01),
    dist.fcts = "euclidean"
  )

  # Extract winning units for each data point
  winning_units <- som_model$unit.classif

  # Create a data frame with cell names and sub-cluster assignments
  tmp_df <- tibble::tibble(cell_name = rownames(tmp_mat),
                           sub_clust = winning_units)

  # Run UMAP for visualization
  set.seed(77)
  tmp_umap <- uwot::umap(tmp_mat)

  # Run_means for sub-clusters
  tmp_Means <- Run_means(matrix = tmp_mat %>% t(), cell_table = tmp_df, cluster = "sub_clust")

  # Run_means_hm for sub-clusters
  tmp_means_hm <- Run_means_hm(matrix = tmp_mat %>% t(), cell_table = tmp_df, cluster = "sub_clust")
  colnames(tmp_means_hm) <- c(1:dim(tmp_means_hm)[2])

  # Run DEPs for the sub-clusters
  tmp_DEPs <- Run_DEPs(matrix = tmp_mat %>% t(), cell_table = tmp_df, cluster = "sub_clust")

  # Run DEPs hm
  tmp_DEPs_hm <- Run_DEPs_hm(matrix = tmp_mat %>% t(), cell_table = tmp_df, cluster = "sub_clust")
  colnames(tmp_DEPs_hm) <- c(1:dim(tmp_DEPs_hm)[2])
  tmp_zDEPs <- t(scale(t(tmp_DEPs_hm), center = TRUE, scale = TRUE))
  # tmp_zDEPs <- pmin(tmp_zDEPs, 3)

  # Create a UMAP plot
  tmp_ggplot <- ggplot2::ggplot(as.data.frame(tmp_umap), ggplot2::aes(x=V1, y=V2, color=as.factor(tmp_df$sub_clust))) +
    ggrastr::geom_point_rast(size = 1, alpha = 0.8) +  # Use rasterized points
    # directlabels::geom_dl(ggplot2::aes(label = as.factor(tmp_df$sub_clust)), method = "smart.grid") +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=6, alpha = 1))) +
    ggplot2::labs(color = "cluster", x = "UMAP_1", y = "UMAP_2") +
    dark_theme_gray()

  tmp_umap_sample <- tmp_umap %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "cell_name") %>%
    dplyr::left_join(cells, by = "cell_name") %>%
    tibble::as_tibble() %>%
    ggplot2::ggplot(., ggplot2::aes(x=V1, y=V2, color=sample)) +
    ggplot2::geom_point(size=0.5, alpha=0.5) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=6, alpha = 1))) +
    ggplot2::labs(color = "sample", x = "UMAP_1", y = "UMAP_2") +
    dark_theme_gray()

  # Create hm
  tmp_hm_DEPs <- ComplexHeatmap::Heatmap(tmp_DEPs_hm %>% as.matrix(),
                                         name = "score",
                                         col = circlize::colorRamp2(breaks = seq(color_min, color_max, length.out = 11), hotmap),
                                         column_km = 0,
                                         row_km = 0,
                                         row_names_side = "left",
                                         row_dend_side = "left",
                                         column_names_rot = 0)

  tmp_hm_means <- ComplexHeatmap::Heatmap(tmp_means_hm %>% as.matrix(),
                                          name = "score",
                                          col = circlize::colorRamp2(breaks = seq(color_min, color_max, length.out = 11), hotmap),
                                          column_km = 0,
                                          row_km = 0,
                                          row_names_side = "left",
                                          row_dend_side = "left",
                                          column_names_rot = 0)

  # Return the results as a list
  return(list(df = tmp_df,
              Means = tmp_Means,
              DEPs = tmp_DEPs,
              DEPs_hm = tmp_DEPs_hm,
              zDEPs = tmp_zDEPs,
              UMAP = tmp_umap,
              UMAP_plot_clust = tmp_ggplot,
              UMAP_plot_sample = tmp_umap_sample,
              hm_DEPs = tmp_hm_DEPs,
              hm_means = tmp_hm_means
  )
  )
}


train_som <- function(df, map_size = 10, map_shape = "rectangular", t_training = 10, learn_rate = c(0.05, 0.01)) {
  set.seed(5410)
  normalized_data_mat <- as.matrix(df)

  grid_size <- kohonen::somgrid(xdim = map_size, ydim = map_size, topo = map_shape)

  # Initialize SOM with desired grid size.
  som_model <- kohonen::som(
    normalized_data_mat,
    grid = grid_size,
    rlen = t_training,
    alpha = learn_rate,
    dist.fcts = "euclidean"
  )

  # Extract codebook vectors and flatten them
  codebook_vectors <- som_model$codes
  flattened_codebook <- do.call(rbind, codebook_vectors)

  # Find the nearest node for each data point
  nearest_node_indices <- data.frame(CellID = rownames(normalized_data_mat),
                                     SOM = FNN::knnx.index(flattened_codebook, normalized_data_mat, k = 1))

  # Extract winning units for each data point
  winning_units <- som_model$unit.classif

  # Calculate mean expression of cells per 100 clusters
  mean_expression_per_cluster <- sapply(unique(winning_units), function(cluster) {
    cluster_indices <- which(winning_units == cluster)
    # If there is only one point in the cluster, its values are the mean expression
    if(length(cluster_indices) == 1){
      cluster_data <- normalized_data_mat[cluster_indices, ]
      cluster_data
    }
    # Otherwise, calculate mean expression per cluster
    else{
      cluster_data <- normalized_data_mat[cluster_indices, ]
      colMeans(cluster_data)
    }
  })
  colnames(mean_expression_per_cluster) <- as.character(unique(winning_units))

  # Return results
  return(list(nearest_node_indices = nearest_node_indices,
              winning_units = winning_units,
              mean_expression_per_cluster = mean_expression_per_cluster))
}


#' Merge SOM clusters based on correlation threshold
#' @param som_results List containing SOM results (nearest_node_indices, winning_units, mean_expression_per_cluster)
#' @param cor_thresh Correlation threshold above which clusters will be merged
#' @return List containing merged clusters info and cell assignments
merge_som_clusters <- function(som_results, cor_thresh) {
  # Calculate correlation matrix from mean expression
  cor_matrix <- cor(som_results$mean_expression_per_cluster,
                    use = "pairwise.complete.obs",
                    method = "pearson")

  # Initialize variables
  current_matrix <- cor_matrix
  cluster_members <- as.list(colnames(cor_matrix))
  names(cluster_members) <- colnames(cor_matrix)

  while(TRUE) {
    # Set diagonal to -1 to ignore self-correlations
    diag(current_matrix) <- -1

    # Find maximum correlation
    max_cor <- max(current_matrix)

    # If no correlations above threshold, break
    if(max_cor <= cor_thresh) {
      break
    }

    # Find indices of clusters to merge
    max_indices <- which(current_matrix == max_cor, arr.ind = TRUE)[1,]
    cluster1 <- rownames(current_matrix)[max_indices[1]]
    cluster2 <- colnames(current_matrix)[max_indices[2]]

    # Create temporary new cluster name
    new_cluster <- paste(cluster1, cluster2, sep="_")

    # Update cluster members
    cluster_members[[new_cluster]] <- c(
      cluster_members[[cluster1]],
      cluster_members[[cluster2]]
    )
    cluster_members[[cluster1]] <- NULL
    cluster_members[[cluster2]] <- NULL

    # Create new correlation matrix excluding merged clusters
    remaining_clusters <- setdiff(rownames(current_matrix), c(cluster1, cluster2))
    new_matrix <- current_matrix[remaining_clusters, remaining_clusters, drop=FALSE]

    # Calculate correlations for new merged cluster
    if(length(remaining_clusters) > 0) {
      # Calculate mean correlation with remaining clusters
      new_cors <- (current_matrix[remaining_clusters, cluster1] +
                     current_matrix[remaining_clusters, cluster2]) / 2

      # Add new cluster to matrix
      new_matrix <- cbind(new_matrix, new_cors)
      new_matrix <- rbind(new_matrix, c(new_cors, 1))
      rownames(new_matrix)[nrow(new_matrix)] <- new_cluster
      colnames(new_matrix)[ncol(new_matrix)] <- new_cluster
    } else {
      # If no remaining clusters, create 1x1 matrix for final cluster
      new_matrix <- matrix(1, nrow=1, ncol=1)
      rownames(new_matrix) <- new_cluster
      colnames(new_matrix) <- new_cluster
    }

    current_matrix <- new_matrix
  }

  # Create mapping from original to new clusters
  cluster_mapping <- list()
  for(i in seq_along(cluster_members)) {
    new_cluster_id <- i
    original_clusters <- sort(as.numeric(strsplit(names(cluster_members)[i], "_")[[1]]))
    for(orig_cluster in original_clusters) {
      cluster_mapping[[as.character(orig_cluster)]] <- new_cluster_id
    }
  }

  # Update cell assignments
  cell_assignments <- data.frame(
    CellID = som_results$nearest_node_indices$CellID,
    original_cluster = som_results$nearest_node_indices$SOM,
    new_cluster = sapply(som_results$nearest_node_indices$SOM,
                         function(x) cluster_mapping[[as.character(x)]])
  )

  # Calculate new mean expression per merged cluster
  new_mean_expression <- matrix(0,
                                nrow = nrow(som_results$mean_expression_per_cluster),
                                ncol = length(unique(cell_assignments$new_cluster)))

  rownames(new_mean_expression) <- rownames(som_results$mean_expression_per_cluster)
  colnames(new_mean_expression) <- paste0("cluster_", 1:ncol(new_mean_expression))

  for(i in 1:ncol(new_mean_expression)) {
    original_clusters <- names(cluster_mapping)[sapply(cluster_mapping, function(x) x == i)]
    new_mean_expression[,i] <- rowMeans(som_results$mean_expression_per_cluster[, original_clusters, drop=FALSE])
  }

  # Prepare final output
  final_clusters <- list()
  for(i in seq_along(cluster_members)) {
    final_clusters[[i]] <- list(
      members = cell_assignments$CellID[cell_assignments$new_cluster == i],
      original_clusters = sort(as.numeric(strsplit(names(cluster_members)[i], "_")[[1]])),
      size = sum(cell_assignments$new_cluster == i)
    )
  }
  names(final_clusters) <- paste0("cluster_", seq_along(final_clusters))

  return(list(
    clusters = final_clusters,
    cell_assignments = cell_assignments,
    mean_expression = new_mean_expression,
    correlation_threshold = cor_thresh
  ))
}

# Example usage:
# merged_results <- merge_som_clusters(som_results, cor_thresh = 0.7)
#
# # Access results
# head(merged_results$cell_assignments)  # Show first few cell assignments
# print(merged_results$clusters)         # Show cluster information
# dim(merged_results$mean_expression)    # Show dimensions of new expression matrix



# Cluster Plotting --------------------------------------------------------

Plot_hierarchical_clust <- function(matrix = rmat,
                                    marker = marker_mac,
                                    cluster = "gating",
                                    cell_types = c("Mac_CD206", "Inf_Mac", "Mic", "Mac"),
                                    slice_size = 2000,
                                    color_min = 0,
                                    color_max = 0.9,
                                    n_clusters = 6,
                                    cell_table = cells) {

  cell_table <- cells

  # Step 1: Subset and subsample the matrix
  tmp <- matrix[marker, cell_table %>% dplyr::filter(!!rlang::sym(cluster) %in% cell_types) %>% dplyr::pull(cell_name)]
  set.seed(77)
  tmp <- tmp[, sample(colnames(tmp), size = slice_size, replace = FALSE)] # Subsample

  # Step 2: Create the column annotation
  an_sample <- ComplexHeatmap::HeatmapAnnotation(an_sample = cell_table %>% dplyr::filter(cell_name %in% colnames(tmp)) %>% dplyr::pull(sample),
                                                 which = "column",
                                                 annotation_label = "sample",
                                                 annotation_legend_param = list(border = "black"),
                                                 border = TRUE,
                                                 col = list(an_sample = an_cols_sample))

  # Step 3: Generate the main heatmap object
  ht <- ComplexHeatmap::Heatmap(tmp,
                                name = "score",
                                col = circlize::colorRamp2(breaks = seq(color_min, color_max, length.out = 11), hotmap),
                                row_names_side = "left",
                                row_dend_side = "left",
                                show_column_names = FALSE,
                                clustering_method_rows = "ward.D",
                                clustering_method_columns = "ward.D",
                                column_title_side = "bottom",
                                column_split = n_clusters,
                                border = TRUE,
                                top_annotation = an_sample)

  # Step 4: Get the column indices for each group defined by column_split
  column_groups <- ComplexHeatmap::column_order(ht)

  # Step 5: Compute the average marker expression for each group
  avg_expression <- sapply(column_groups, function(cols) {
    rowMeans(tmp[, cols, drop = FALSE])
  })

  colnames(avg_expression) <- as.character(1:dim(avg_expression)[2])

  # Step 6: Create the heatmap for the averaged data
  avg_heatmap <- ComplexHeatmap::Heatmap(avg_expression,
                                         name = "avg expression",
                                         col = circlize::colorRamp2(breaks = seq(color_min, color_max, length.out = 11), hotmap),
                                         cluster_columns = FALSE,  # Do not cluster columns since they represent groups
                                         cluster_rows = FALSE,     # Optional: keep rows unclustered
                                         show_column_names = TRUE,
                                         column_names_rot = 0,    # Rotate column names by 90 degrees
                                         show_row_names = TRUE,
                                         width = grid::unit(4, "cm"),
                                         border = TRUE)

  # Step 7: Combine the two heatmaps
  ht + avg_heatmap
}


plot_cluster_expression <- function(data_list, raw_expression_matrix,
                                    expression_threshold = 0.1,
                                    cluster_rows = TRUE, cluster_cols = TRUE,
                                    clustering_method = "ward.D2",
                                    cluster_labels = NULL, marker_order = NULL) {

  # Extract data components
  cell_assignments <- data_list$cell_assignments
  mean_expression_matrix <- data_list$mean_expression

  # Ensure row and column names exist
  if (is.null(rownames(raw_expression_matrix)) || is.null(colnames(raw_expression_matrix))) {
    stop("Raw expression matrix must have row and column names.")
  }
  if (is.null(rownames(mean_expression_matrix)) || is.null(colnames(mean_expression_matrix))) {
    stop("Mean expression matrix must have row and column names.")
  }

  # Compute fraction of expressing cells per cluster
  is_expressed <- raw_expression_matrix > expression_threshold

  # Merge cell assignments with expression data
  cell_cluster_mapping <- cell_assignments %>%
    dplyr::select(CellID, new_cluster) %>%
    dplyr::mutate(CellID = as.character(CellID))

  expr_df <- as.data.frame(t(is_expressed))
  expr_df$CellID <- rownames(expr_df)
  expr_df <- merge(expr_df, cell_cluster_mapping, by = "CellID")

  fraction_cells <- expr_df %>%
    dplyr::group_by(new_cluster) %>%
    dplyr::summarise(dplyr::across(where(is.logical), mean)) %>%
    tibble::column_to_rownames(var = "new_cluster")

  # Check if cluster names are numeric and if so add cluster_ prefix
  cluster_types <- sapply(rownames(fraction_cells), function(x) {
    if (grepl("^[0-9]+$", x)) {
      return("numeric")
    } else {
      return("character")
    }
  })
  rownames(fraction_cells) <- ifelse(cluster_types == "numeric",
                                     paste0("cluster_", rownames(fraction_cells)),
                                     rownames(fraction_cells))
  colnames(fraction_cells) <- rownames(mean_expression_matrix)

  # Normalize marker expression
  min_vals <- apply(mean_expression_matrix, 1, min)
  max_vals <- apply(mean_expression_matrix, 1, max)
  norm_expression <- sweep(mean_expression_matrix, 1, min_vals, "-")
  norm_expression <- sweep(norm_expression, 1, max_vals - min_vals, "/")

  # Perform hierarchical clustering if enabled
  marker_order_default <- rownames(norm_expression)
  cluster_order <- colnames(norm_expression)

  if (cluster_rows) {
    hc_rows <- hclust(dist(norm_expression), method = clustering_method)
    marker_order_default <- rownames(norm_expression)[hc_rows$order]
  }
  if (cluster_cols) {
    hc_cols <- hclust(dist(t(norm_expression)), method = clustering_method)
    cluster_order <- colnames(norm_expression)[hc_cols$order]
  }

  # If marker_order is provided, override the computed order
  if (!is.null(marker_order)) {
    if (!all(marker_order %in% marker_order_default)) {
      stop("Some markers in 'marker_order' are not found in the dataset.")
    }
    marker_order_default <- marker_order
  }

  # Reorder matrices based on the determined order
  norm_expression <- norm_expression[marker_order_default, cluster_order]
  fraction_cells <- fraction_cells[cluster_order, marker_order_default]

  # Convert matrices to long format
  df_expr <- reshape2::melt(as.data.frame(norm_expression), variable.name = "Cluster", value.name = "MeanExpression")
  df_expr$Marker <- rownames(norm_expression)

  df_fraction <- fraction_cells %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Cluster") %>%
    reshape2::melt(id.vars = "Cluster", variable.name = "Marker", value.name = "Fraction")

  df_plot <- merge(df_expr, df_fraction, by = c("Cluster", "Marker"))

  # Apply marker order to x-axis
  if (is.null(marker_order)) {
    df_plot$Marker <- factor(df_plot$Marker, levels = rev(marker_order_default))  # Reverse order
  } else {
    df_plot$Marker <- factor(df_plot$Marker, levels = marker_order)  # Keep given order
  }

  if (!is.null(cluster_labels)) {
    # Create a mapping of generic names to labels
    cluster_mapping <- data.frame(
      generic = cluster_order,  # "cluster_1", "cluster_2"
      label = cluster_labels[cluster_order],  # "Epithelial", "Fibroblast"
      stringsAsFactors = FALSE
    )

    # Ensure uniqueness without visible numbering
    label_counts <- table(cluster_mapping$label)
    duplicated_labels <- names(label_counts[label_counts > 1])

    # Append an invisible Unicode character for duplicates
    for (dup in duplicated_labels) {
      idx <- which(cluster_mapping$label == dup)
      invisible_spaces <- c("", "\u200B", "\u200C", "\u200D", "\uFEFF")  # Zero-width characters
      cluster_mapping$label[idx] <- paste0(dup, invisible_spaces[seq_along(idx)])
    }

    # Merge df_plot with the cluster mapping to ensure proper labeling
    df_plot <- df_plot %>%
      dplyr::left_join(cluster_mapping, by = c("Cluster" = "generic")) %>%
      dplyr::mutate(Cluster = factor(label, levels = unique(label)))  # Keep correct order

  } else {
    df_plot$Cluster <- factor(df_plot$Cluster, levels = cluster_order)
  }

  # Plot
  ggplot2::ggplot(df_plot, ggplot2::aes(x = Marker, y = Cluster, size = Fraction * 100, fill = MeanExpression)) +
    ggplot2::geom_point(shape = 21, stroke = 0.5, color = "gray") +
    ggplot2::scale_size_continuous(range = c(1, 10), limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
    ggplot2::scale_fill_gradientn(colours = c("#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f")) +
    ggplot2::theme_minimal(base_size = 14) +  # Increase base font size
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.98, vjust = 0.5, size = 14),  # Increase axis text size
      axis.text.y = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 14),
      panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 1),
      axis.ticks = ggplot2::element_line(color = "black", size = 0.8)) +
    ggplot2::labs(x = "", y = "", size = "Fraction of cells (%)", fill = "Mean Expression") +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(order = 1, position = "right", frame.colour = "black", ticks.colour = "black"),
      size = ggplot2::guide_legend(order = 2, override.aes = list(fill = "black"))
    )
}



plot_marker_by_celltype <- function(
    rmat,
    cell_assignments,
    markers             = NULL,
    cell_types          = NULL,
    expression_threshold = 0,         # number, named vector, df, or "auto"
    auto_ref            = NULL,       # named vector marker→cell_type, required if expression_threshold=="auto"
    norm                = FALSE,      # if TRUE, color‐scale the mean to [0,1]
    cluster_rows        = TRUE,
    cluster_cols        = TRUE,
    clustering_method   = "ward.D2",
    marker_order        = NULL,
    celltype_order      = NULL,
    mean_score_max = FALSE
) {
  # 1) Sanity checks
  if (!is.matrix(rmat))
    stop("`rmat` must be a matrix (markers × CellIDs).")
  if (!all(c("CellID","cell_type") %in% colnames(cell_assignments)))
    stop("`cell_assignments` needs columns CellID + cell_type.")
  if (!all(cell_assignments$CellID %in% colnames(rmat)))
    stop("Some CellID in `rmat` are missing from `cell_assignments`.")

  # 2) Select markers
  if (is.null(markers)) {
    markers <- rownames(rmat)
  } else {
    miss <- setdiff(markers, rownames(rmat))
    if (length(miss)) stop("Markers not found: ", paste(miss, collapse=", "))
  }
  if (!is.null(marker_order)) {
    if (!all(marker_order %in% markers))
      stop("`marker_order` has unknown markers.")
    markers <- marker_order
  }
  if (identical(expression_threshold, "auto")) {
    # only keep markers that we have an auto_ref for
    markers <- intersect(markers, names(auto_ref))
    if (length(markers)==0) stop("No markers left after auto‐pruning.")
  }

  # 3) Select cell types
  all_ct <- unique(cell_assignments$cell_type)

  if (is.null(cell_types)) {
    cell_types <- all_ct
  } else {
    miss <- setdiff(cell_types, all_ct)
    if (length(miss)) stop("cell_types not found in assignments: ", paste(miss, collapse = ", "))
  }

  if (!is.null(celltype_order)) {
    # Warn about any requested groups that aren't present
    missing_ct <- setdiff(celltype_order, all_ct)
    if (length(missing_ct)) {
      warning("These cell types from `celltype_order` are not present and will be ignored: ",
              paste(missing_ct, collapse = ", "))
    }
    present_ct <- intersect(celltype_order, all_ct)
    if (length(present_ct) == 0)
      stop("None of the `celltype_order` groups are present in the data.")

    # Use only the requested (and present) groups; this also controls Y-axis order
    cell_types     <- intersect(cell_types, present_ct)  # respect any `cell_types` filter too
    celltype_lvls  <- present_ct
  } else {
    # No explicit order supplied; keep whatever `cell_types` resolved to
    celltype_lvls <- sort(unique(cell_types))
  }

  # 4) Long‐format all scores
  df_long <- as.data.frame(
    t(rmat[markers, , drop=FALSE]), stringsAsFactors=FALSE
  ) %>%
    tibble::rownames_to_column("CellID") %>%
    dplyr::inner_join(cell_assignments, by="CellID") %>%
    dplyr::filter(cell_type %in% cell_types) %>%
    tidyr::pivot_longer(
      cols     = all_of(markers),
      names_to = "Marker",
      values_to= "Score"
    )

  # 5) Build threshold table thr_df(marker, threshold)
  if (identical(expression_threshold, "auto")) {
    # require a named vector auto_ref of length==length(markers)
    if (is.null(auto_ref) || !is.character(auto_ref) ||
        length(auto_ref)!=length(markers) ||
        is.null(names(auto_ref)) ||
        !all(markers %in% names(auto_ref)))
      stop("When `expression_threshold='auto'`, supply named vector `auto_ref` of length(markers) mapping each marker → reference cell_type.")

    thr_df <- tibble::tibble(Marker = markers) %>%
      dplyr::mutate(
        threshold = purrr::map_dbl(Marker, function(m) {
          ref_ct  <- auto_ref[m]
          # scores of that marker within that reference cell type:
          sc_ref  <- df_long$Score[
            df_long$Marker==m & df_long$cell_type==ref_ct
          ]
          # 10th percentile:
          quantile(sc_ref, 0.1, na.rm=TRUE)
        })
      )

  } else if (is.data.frame(expression_threshold)) {
    thr_df <- expression_threshold %>% dplyr::select(Marker, threshold)

  } else if (length(expression_threshold)==1) {
    thr_df <- tibble::tibble(
      Marker    = markers,
      threshold = as.numeric(expression_threshold)
    )

  } else if (is.numeric(expression_threshold) &&
             length(expression_threshold)==length(markers)) {
    vec <- expression_threshold
    if (!is.null(names(vec))) vec <- vec[markers]
    thr_df <- tibble::tibble(
      Marker    = markers,
      threshold = as.numeric(vec)
    )

  } else {
    stop("`expression_threshold` must be 'auto', scalar, length==n_markers, or df.")
  }

  df_long <- df_long %>% dplyr::left_join(thr_df, by="Marker")

  # 6) Summarise raw mean & fraction(Score>threshold) — dot‐sizes use this
  summary_df <- df_long %>%
    dplyr::group_by(cell_type, Marker) %>%
    dplyr::summarise(
      MeanScore = mean(Score, na.rm=TRUE),
      Fraction  = mean(Score > threshold, na.rm=TRUE)*100,
      .groups   = "drop"
    )

  # 7) Optionally normalize MeanScore (for color)
  if (norm) {
    summary_df <- summary_df %>%
      dplyr::group_by(Marker) %>%
      dplyr::mutate(
        MeanScore = (MeanScore - min(MeanScore,na.rm=TRUE)) /
          (max(MeanScore,na.rm=TRUE) - min(MeanScore,na.rm=TRUE))
      ) %>%
      dplyr::ungroup()
    fill_label <- "Normalized\nmean score"
  } else {
    fill_label <- "Mean score"
  }

  # 8) Build mean_mat for clustering
  mean_mat <- summary_df %>%
    tidyr::pivot_wider(
      id_cols     = Marker,
      names_from  = cell_type,
      values_from = MeanScore,
      values_fn   = mean,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames("Marker") %>%
    as.matrix()

  # 9) Compute row/col ordering
  mo <- rownames(mean_mat)
  co <- colnames(mean_mat)
  if (cluster_rows) {
    hr <- hclust(dist(mean_mat), method=clustering_method)
    mo <- mo[hr$order]
  }
  if (cluster_cols) {
    hc <- hclust(dist(t(mean_mat)), method=clustering_method)
    co <- co[hc$order]
  }
  if (!is.null(marker_order))   mo <- marker_order
  if (!is.null(celltype_order)) co <- celltype_order

  # 10) Prepare for plotting
  plot_df <- summary_df %>%
    dplyr::filter(Marker %in% mo, cell_type %in% co) %>%       # NEW: filter to requested levels
    dplyr::mutate(
      Marker    = factor(Marker,    levels=mo),
      cell_type = factor(cell_type, levels=co)
    )

  # 11) Draw
  fill_limits <- NULL
  if (!identical(mean_score_max, FALSE)) {
    if (!is.numeric(mean_score_max) || length(mean_score_max) != 1 || is.na(mean_score_max)) {
      stop("`mean_score_max` must be FALSE or a single numeric value.")
    }
    fill_limits <- c(0, mean_score_max)
  }

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = Marker, y = cell_type, size = Fraction, fill = MeanScore)
  ) +
    ggplot2::geom_point(shape = 21, color = "gray50", stroke = 0.5) +
    ggplot2::scale_size_continuous(
      range  = c(1,10),
      limits = c(0,100),
      breaks = seq(0,100,25),
      name   = "Fraction of cells (%)"
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c("#f7f7f7","#fddbc7","#f4a582",
                  "#d6604d","#b2182b","#67001f"),
      name    = fill_label,
      limits  = fill_limits,
      oob     = scales::squish
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 90, hjust = 0.98, vjust = 0.5, size = 14),
      axis.text.y  = ggplot2::element_text(face = "bold", size = 14),
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 14),
      panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 1),
      axis.ticks = ggplot2::element_line(color = "black", size = 0.8)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(order = 1, position = "right", frame.colour = "black", ticks.colour = "black"),
      size = ggplot2::guide_legend(order = 2, override.aes = list(fill = "black"))
    )
}



plot_marker_sets_by_celltype <- function(
    rmat,
    cell_assignments,
    marker_sets,                         # named list. Each entry can be:
    #  - character vector (pos markers)
    #  - list(pos=..., neg=...) for negative gating
    cell_types               = NULL,
    expression_threshold      = 0,        # pos threshold: scalar / named vec / df / "auto"
    expression_threshold_neg  = NULL,     # neg threshold: defaults to expression_threshold if NULL
    auto_ref                 = NULL,      # named vec Signature -> reference cell_type for auto thresholds
    auto_pos_quantile        = 0.1,
    auto_neg_quantile        = 0.9,

    # --- how to compute dot size ("Fraction") ---
    fraction_rule            = c("pos_only", "pos_and_neg", "combined"),
    # pos_only      : % cells with pos_score > thr_pos
    # pos_and_neg   : % cells with pos_score > thr_pos AND neg_score < thr_neg  (recommended for negative gating)
    # combined      : % cells with combined_score > thr_pos   (combined_score = pos - w*neg)

    # --- combined score options ---
    neg_weight               = 1,        # combined_score = pos_score - neg_weight*neg_score

    ### ============ fill scaling ============
    fill_mode                = c("mean", "log2fc", "z", "softmax"),
    clip_negative            = TRUE,
    softmax_temp             = 0.05,
    norm                     = FALSE,    # only used when fill_mode=="mean": normalize per Signature to [0,1]
    mean_score_max           = FALSE,    # FALSE or numeric cap (squish)

    ### ============ ordering ============
    cluster_rows             = TRUE,
    cluster_cols             = TRUE,
    clustering_method        = "ward.D2",
    set_order                = NULL,
    celltype_order           = NULL,

    ### ============ robustness ============
    case_insensitive_markers = TRUE,
    min_markers_present      = 1,
    eps                      = 1e-6
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(scales)
  })

  fraction_rule <- match.arg(fraction_rule)
  fill_mode     <- match.arg(fill_mode)

  # ---------- checks ----------
  if (!is.matrix(rmat) && !inherits(rmat, "Matrix"))
    stop("`rmat` must be a base matrix or Matrix::Matrix (markers × CellIDs).")
  if (!all(c("CellID","cell_type") %in% colnames(cell_assignments)))
    stop("`cell_assignments` needs columns: CellID + cell_type.")
  if (!is.list(marker_sets) || is.null(names(marker_sets)) || any(!nzchar(names(marker_sets))))
    stop("`marker_sets` must be a named list: Signature -> markers or list(pos,neg).")

  # Map CellIDs to cell_type aligned to rmat columns
  ct_map <- cell_assignments$cell_type
  names(ct_map) <- cell_assignments$CellID
  ct <- unname(ct_map[colnames(rmat)])

  # Select cell types
  all_ct <- sort(unique(cell_assignments$cell_type))
  if (is.null(cell_types)) cell_types <- all_ct
  miss <- setdiff(cell_types, all_ct)
  if (length(miss)) stop("cell_types not found: ", paste(miss, collapse=", "))

  if (!is.null(celltype_order)) {
    present_ct <- intersect(celltype_order, all_ct)
    if (!length(present_ct)) stop("None of the `celltype_order` groups are present.")
    cell_types <- intersect(cell_types, present_ct)
    celltype_lvls <- present_ct
  } else {
    celltype_lvls <- sort(unique(cell_types))
  }

  keep <- !is.na(ct) & (ct %in% cell_types)
  if (!any(keep)) stop("No cells left after filtering to requested cell_types.")
  idx_by_ct <- split(which(keep), ct[keep])  # cell_type -> indices
  keep_idx  <- which(keep)

  # ---------- marker matching ----------
  available <- rownames(rmat)
  match_markers <- function(markers) {
    markers <- unique(markers)
    if (!case_insensitive_markers) return(intersect(markers, available))
    lut <- setNames(available, tolower(available))
    matched <- unname(lut[tolower(markers)])
    matched <- matched[!is.na(matched)]
    unique(matched)
  }

  # Parse marker_sets into pos/neg lists
  parse_set <- function(x) {
    if (is.character(x)) return(list(pos = x, neg = character(0)))
    if (is.list(x) && !is.null(x$pos)) {
      neg <- if (!is.null(x$neg)) x$neg else character(0)
      return(list(pos = x$pos, neg = neg))
    }
    stop("Each marker_sets[[Signature]] must be a character vector (pos) or list(pos=..., neg=...).")
  }

  parsed_sets <- lapply(marker_sets, parse_set)
  pos_markers <- lapply(parsed_sets, \(s) match_markers(s$pos))
  neg_markers <- lapply(parsed_sets, \(s) match_markers(s$neg))

  # Drop signatures with too few pos markers
  npos <- vapply(pos_markers, length, integer(1))
  bad  <- names(npos)[npos < min_markers_present]
  if (length(bad)) {
    warning("Dropping signatures with < min_markers_present POS markers: ", paste(bad, collapse=", "))
    pos_markers <- pos_markers[npos >= min_markers_present]
    neg_markers <- neg_markers[npos >= min_markers_present]
  }
  sig_names <- names(pos_markers)
  if (!length(sig_names)) stop("No marker sets left after matching to rmat rows.")

  if (!is.null(set_order)) {
    sig_names <- intersect(set_order, sig_names)
    if (!length(sig_names)) stop("None of the `set_order` signatures are present.")
    pos_markers <- pos_markers[sig_names]
    neg_markers <- neg_markers[sig_names]
  }

  # ---------- threshold resolution ----------
  # Helpers: make a per-signature numeric threshold vector for pos/neg
  resolve_thresholds <- function(thr, kind = c("pos","neg")) {
    kind <- match.arg(kind)

    if (identical(thr, "auto")) {
      if (is.null(auto_ref) || is.null(names(auto_ref)))
        stop("When threshold is 'auto', provide named vector `auto_ref`: Signature -> reference cell_type.")
      if (!all(sig_names %in% names(auto_ref)))
        stop("`auto_ref` missing signatures: ", paste(setdiff(sig_names, names(auto_ref)), collapse=", "))
      return(list(mode="auto", ref=auto_ref[sig_names]))
    }

    if (is.data.frame(thr)) {
      # support df with Signature + threshold, and optionally threshold_neg
      if (!("Signature" %in% names(thr))) stop("threshold df must include column 'Signature'.")
      if (kind == "pos") {
        if (!("threshold" %in% names(thr))) stop("pos threshold df must include column 'threshold'.")
        v <- thr$threshold; names(v) <- thr$Signature
      } else {
        coln <- if ("threshold_neg" %in% names(thr)) "threshold_neg" else "threshold"
        v <- thr[[coln]]; names(v) <- thr$Signature
      }
      if (!all(sig_names %in% names(v)))
        stop("threshold df missing signatures: ", paste(setdiff(sig_names, names(v)), collapse=", "))
      return(list(mode="fixed", value=as.numeric(v[sig_names])))
    }

    if (length(thr) == 1) {
      v <- rep(as.numeric(thr), length(sig_names)); names(v) <- sig_names
      return(list(mode="fixed", value=v))
    }

    if (is.numeric(thr) && !is.null(names(thr))) {
      if (!all(sig_names %in% names(thr)))
        stop("Named threshold missing signatures: ", paste(setdiff(sig_names, names(thr)), collapse=", "))
      return(list(mode="fixed", value=as.numeric(thr[sig_names])))
    }

    if (is.numeric(thr) && length(thr) == length(sig_names)) {
      v <- as.numeric(thr); names(v) <- sig_names
      return(list(mode="fixed", value=v))
    }

    stop("Threshold must be 'auto', scalar, named vector, length==n_signatures, or df(Signature, threshold[,threshold_neg]).")
  }

  thr_pos <- resolve_thresholds(expression_threshold, "pos")

  if (is.null(expression_threshold_neg)) expression_threshold_neg <- expression_threshold
  thr_neg <- resolve_thresholds(expression_threshold_neg, "neg")

  # Matrix-safe colMeans
  colMeans_any <- function(m) if (inherits(m, "Matrix")) Matrix::colMeans(m) else colMeans(m)

  # ---------- compute per-signature summaries ----------
  summary_list <- vector("list", length(sig_names))
  names(summary_list) <- sig_names

  for (s in sig_names) {
    mk_pos <- pos_markers[[s]]
    mk_neg <- neg_markers[[s]]

    pos_sc <- colMeans_any(rmat[mk_pos, , drop = FALSE])

    has_neg <- length(mk_neg) > 0
    neg_sc <- if (has_neg) colMeans_any(rmat[mk_neg, , drop = FALSE]) else rep(0, length(pos_sc))

    combined_sc <- pos_sc - neg_weight * neg_sc

    # thresholds (pos / neg) for this signature
    get_thr <- function(thr_obj, sc_vec, q) {
      if (thr_obj$mode == "fixed") return(thr_obj$value[[s]])
      ref_ct <- unname(thr_obj$ref[[s]])
      ref_idx <- idx_by_ct[[ref_ct]]
      if (is.null(ref_idx) || !length(ref_idx))
        stop("Auto threshold failed: reference cell_type '", ref_ct, "' has no cells after filtering.")
      as.numeric(quantile(sc_vec[ref_idx], q, na.rm = TRUE))
    }

    thrP <- get_thr(thr_pos, pos_sc, auto_pos_quantile)
    thrN <- if (has_neg) get_thr(thr_neg, neg_sc, auto_neg_quantile) else NA_real_

    tmp <- lapply(names(idx_by_ct), function(ct_name) {
      idx_in   <- idx_by_ct[[ct_name]]
      idx_rest <- setdiff(keep_idx, idx_in)

      mean_in_pos  <- mean(pos_sc[idx_in],      na.rm = TRUE)
      mean_in_neg  <- if (has_neg) mean(neg_sc[idx_in], na.rm = TRUE) else 0
      mean_in_comb <- mean(combined_sc[idx_in], na.rm = TRUE)

      mean_rest_comb <- mean(combined_sc[idx_rest], na.rm = TRUE)

      frac <- switch(
        fraction_rule,
        pos_only = mean(pos_sc[idx_in] > thrP, na.rm = TRUE) * 100,
        combined = mean(combined_sc[idx_in] > thrP, na.rm = TRUE) * 100,
        pos_and_neg = {
          if (!has_neg) {
            mean(pos_sc[idx_in] > thrP, na.rm = TRUE) * 100
          } else {
            mean((pos_sc[idx_in] > thrP) & (neg_sc[idx_in] < thrN), na.rm = TRUE) * 100
          }
        }
      )

      data.frame(
        cell_type     = ct_name,
        Signature     = s,
        MeanPos       = mean_in_pos,
        MeanNeg       = mean_in_neg,
        MeanScore     = mean_in_comb,          # <-- default "mean" uses combined score
        MeanRest      = mean_rest_comb,        # needed for log2fc
        Fraction      = frac,
        thr_pos       = thrP,
        thr_neg       = thrN,
        stringsAsFactors = FALSE
      )
    }) %>% dplyr::bind_rows()

    summary_list[[s]] <- tmp
  }

  summary_df <- dplyr::bind_rows(summary_list)

  # ---------- fill value from MeanScore with modes ----------
  summary_df <- summary_df %>% dplyr::mutate(FillValue = MeanScore)

  if (fill_mode == "mean") {
    if (norm) {
      summary_df <- summary_df %>%
        dplyr::group_by(Signature) %>%
        dplyr::mutate(
          FillValue = {
            mn <- min(MeanScore, na.rm = TRUE)
            mx <- max(MeanScore, na.rm = TRUE)
            if (is.finite(mx - mn) && (mx - mn) > 0) (MeanScore - mn) / (mx - mn) else 0
          }
        ) %>% dplyr::ungroup()
      fill_label <- "Normalized\ncombined score"
    } else {
      fill_label <- "Combined score\n(pos - w·neg)"
    }
  }

  if (fill_mode == "log2fc") {
    summary_df <- summary_df %>%
      dplyr::mutate(FillValue = log2((MeanScore + eps) / (MeanRest + eps)))
    if (clip_negative) summary_df$FillValue <- pmax(summary_df$FillValue, 0)
    fill_label <- if (clip_negative) "log2 FC vs rest\n(clipped)" else "log2 FC vs rest"
  }

  if (fill_mode == "z") {
    summary_df <- summary_df %>%
      dplyr::group_by(Signature) %>%
      dplyr::mutate(
        FillValue = {
          mu <- mean(MeanScore, na.rm = TRUE)
          sdv <- stats::sd(MeanScore, na.rm = TRUE)
          z <- (MeanScore - mu) / sdv
          z[!is.finite(z)] <- 0
          z
        }
      ) %>% dplyr::ungroup()
    if (clip_negative) summary_df$FillValue <- pmax(summary_df$FillValue, 0)
    fill_label <- if (clip_negative) "Z-score\n(clipped)" else "Z-score"
  }

  if (fill_mode == "softmax") {
    T <- softmax_temp
    if (!is.numeric(T) || length(T) != 1 || is.na(T) || T <= 0)
      stop("`softmax_temp` must be a single positive number.")
    summary_df <- summary_df %>%
      dplyr::group_by(Signature) %>%
      dplyr::mutate(
        FillValue = {
          x <- MeanScore / T
          x <- x - max(x, na.rm = TRUE)
          ex <- exp(x)
          ex / sum(ex, na.rm = TRUE)
        }
      ) %>% dplyr::ungroup()
    fill_label <- paste0("Softmax prob\n(T=", softmax_temp, ")")
  }

  # ---------- ordering based on MeanScore (not FillValue) ----------
  mean_mat <- summary_df %>%
    dplyr::select(Signature, cell_type, MeanScore) %>%
    tidyr::pivot_wider(id_cols = Signature, names_from = cell_type, values_from = MeanScore, values_fill = 0) %>%
    tibble::column_to_rownames("Signature") %>%
    as.matrix()

  mo <- rownames(mean_mat)
  co <- colnames(mean_mat)
  if (cluster_rows) {
    hr <- hclust(dist(mean_mat), method = clustering_method)
    mo <- mo[hr$order]
  }
  if (cluster_cols) {
    hc <- hclust(dist(t(mean_mat)), method = clustering_method)
    co <- co[hc$order]
  }
  if (!is.null(set_order))      mo <- set_order[set_order %in% mo]
  if (!is.null(celltype_order)) co <- celltype_order[celltype_order %in% co]

  plot_df <- summary_df %>%
    dplyr::filter(Signature %in% mo, cell_type %in% co) %>%
    dplyr::mutate(
      Signature = factor(Signature, levels = mo),
      cell_type = factor(cell_type, levels = co)
    )

  # ---------- fill squish limits ----------
  fill_limits <- NULL
  if (!identical(mean_score_max, FALSE)) {
    if (!is.numeric(mean_score_max) || length(mean_score_max) != 1 || is.na(mean_score_max))
      stop("`mean_score_max` must be FALSE or a single numeric value.")
    fill_limits <- c(min(plot_df$FillValue, na.rm = TRUE), mean_score_max)
    if (clip_negative && fill_mode %in% c("log2fc","z")) fill_limits[1] <- 0
    if (fill_mode == "softmax") fill_limits <- c(0, 1)
    if (fill_mode == "mean" && norm) fill_limits <- c(0, 1)
  }

  ggplot2::ggplot(plot_df, ggplot2::aes(x = Signature, y = cell_type, size = Fraction, fill = FillValue)) +
    ggplot2::geom_point(shape = 21, color = "gray50", stroke = 0.5) +
    ggplot2::scale_size_continuous(
      range  = c(1, 10),
      limits = c(0, 100),
      breaks = seq(0, 100, 25),
      name   = "Fraction of cells (%)"
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c("#f7f7f7","#fddbc7","#f4a582","#d6604d","#b2182b","#67001f"),
      name   = fill_label,
      limits = fill_limits,
      oob    = scales::squish
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 90, hjust = 0.98, vjust = 0.5, size = 14),
      axis.text.y  = ggplot2::element_text(face = "bold", size = 14),
      legend.text  = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 14),
      panel.grid   = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 1),
      axis.ticks   = ggplot2::element_line(color = "black", size = 0.8)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(order = 1, position = "right", frame.colour = "black", ticks.colour = "black"),
      size = ggplot2::guide_legend(order = 2, override.aes = list(fill = "black"))
    )
}



# Utility Functions -------------------------------------------------------

`%ni%` <- Negate(`%in%`)


outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}


Count_cell_type <- function(cell_type = cell_type) {
  cells %>% dplyr::count({{cell_type}}) %>% print(n=Inf)
}


Unlist = function(L, nested.names = FALSE) {
  if (nested.names) {
    Names = unlist(sapply(L, names), use.names = F)
  } else {
    Names = rep(names(L), lengths(L))
  }
  stats::setNames(unlist(L), Names)
}


# Color Palettes ----------------------------------------------------------

fig.colors <- c(
  "#e64b35", "#4dbbd6", "#00a086", "#3d5488", "#f29b7f", "#8491b4", "#91d1c1", "#dc0000", "#7f6149", "#F0C976",
  "#0044A2", "#777F49", "#b56f58", "#00833e", "#fbb045", "#a5c2fc", "#3a8ea3", "#7f4019", "#00d1bc", "#ab65bc"
)

dark.fig.colors <- c(
  "#b33a29", "#3a8ea3", "#006d5b", "#263455", "#bf7a64", "#5e6781", "#6d9e91", "#a90000", "#4c3a2b", "#bd9e5c",
  "#002e6e", "#474c2b", "#824f3f", "#004f25", "#c78c36", "#839ac9", "#27616f", "#4c260e", "#009e8e", "#7c4989"
)

hotmap <- c(
  "#053061",  # Very dark blue
  "#2166ac",  # Medium-dark blue
  "#4393c3",  # Medium blue
  "#92c5de",  # Light blue
  "#d1e5f0",  # Very light blue
  "#f7f7f7",  # Light grey (almost white)
  "#fddbc7",  # Very light peach
  "#f4a582",  # Light reddish-orange
  "#d6604d",  # Medium-dark reddish-orange
  "#b2182b",  # Dark red
  "#67001f"   # Very dark red
)

# hotmap <- c(
#   "#2166ac",  # Medium-dark blue
#   "#3399FF",  # Vibrant sky blue
#   "#66CCFF",  # Light blue
#   "#99E5FF",  # Pale bright blue
#   "#CCF5FF",  # Soft cyan
#   "#FFFFFF",  # Pure white
#   "#FFCCB3",  # Bright peach
#   "#FF9966",  # Vibrant orange
#   "#FF5733",  # Bright red-orange
#   "#FF0000",  # Bold bright red
#   "#990000"   # Deep bright red
# )

# Emulate ggplot2's default color palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


umap_pal <- c("#88CCEE", "pink2", "sandybrown", "#AA4499", "#44AA99", "#332288",
              "#DDCC77", "#882255", "#E4CAB5", "#888888", "#CC6677", "#999933",
              "#117733", "#0c2a50ff", "darkgoldenrod", "violetred", "#6699CC",
              "skyblue4", "cyan3", "yellowgreen", "burlywood3", "tan3")



# Plotting Utils ----------------------------------------------------------

flip_y <- function(sf_obj) {
  # Get bounding box
  bb <- sf::st_bbox(sf_obj)
  ymax <- as.numeric(bb["ymax"])
  ymin <- as.numeric(bb["ymin"])

  # Function to flip coordinates in a polygon
  flip_coords <- function(coords) {
    coords[, 2] <- ymax + ymin - coords[, 2]  # Flip y-coordinates
    return(coords)
  }

  # Apply to each geometry
  geom_flipped <- sf::st_geometry(sf_obj)

  for (i in seq_along(geom_flipped)) {
    # Handle MULTIPOLYGON
    if (inherits(geom_flipped[[i]], "MULTIPOLYGON")) {
      for (j in seq_along(geom_flipped[[i]])) {
        for (k in seq_along(geom_flipped[[i]][[j]])) {
          geom_flipped[[i]][[j]][[k]] <- flip_coords(geom_flipped[[i]][[j]][[k]])
        }
      }
    }
  }

  # Rebuild sf object
  result <- sf_obj
  sf::st_geometry(result) <- geom_flipped

  return(result)
}

library(sf)
library(purrr)

flip_y_to_ref <- function(sf_obj, ref_sf) {
  bb   <- sf::st_bbox(ref_sf)          # use the cells bbox (already flipped)
  ymin <- as.numeric(bb["ymin"])
  ymax <- as.numeric(bb["ymax"])

  flip_coords <- function(m) {
    m[,2] <- ymin + ymax - m[,2]
    m
  }

  g <- sf::st_geometry(sf_obj)

  g_flip <- lapply(g, function(geom) {
    cls <- class(geom)[2]  # "POLYGON", "MULTIPOLYGON", "LINESTRING", etc.

    switch(cls,
           "POINT" = {
             m <- matrix(sf::st_coordinates(geom), ncol = 2)
             sf::st_point(flip_coords(m)[1,], dim = "XY")
           },
           "MULTIPOINT" = {
             m <- sf::st_coordinates(geom)[,1:2, drop = FALSE]
             sf::st_multipoint(flip_coords(m))
           },
           "LINESTRING" = {
             m <- sf::st_coordinates(geom)[,1:2, drop = FALSE]
             sf::st_linestring(flip_coords(m))
           },
           "MULTILINESTRING" = {
             parts <- geom
             parts <- lapply(parts, function(ls) sf::st_linestring(flip_coords(ls)))
             sf::st_multilinestring(parts)
           },
           "POLYGON" = {
             rings <- geom
             rings <- lapply(rings, function(r) flip_coords(r))
             sf::st_polygon(rings)
           },
           "MULTIPOLYGON" = {
             polys <- geom
             polys <- lapply(polys, function(pg) {
               lapply(pg, function(r) flip_coords(r))
             })
             sf::st_multipolygon(polys)
           },
           { stop("Geometry type not supported: ", cls) }
    )
  })

  sf::st_geometry(sf_obj) <- sf::st_sfc(g_flip, crs = sf::st_crs(sf_obj))
  sf_obj
}



plot_core <- function(sf_obj) {
  ggplot2::ggplot(sf_obj) +
    ggplot2::geom_sf(ggplot2::aes(fill = assignment), colour = "black", linewidth = 0.02) +
    ggplot2::coord_sf(expand = FALSE, datum = NA) +
    ggplot2::scale_fill_manual(values = my_pal, name = "Cell type", na.translate = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "black", color = NA), panel.background = ggplot2::element_rect(fill = "black", color = NA),
                   legend.text = ggplot2::element_text(color = "white", size = 14), legend.title = ggplot2::element_text(color = "white", size = 16))
}


plot_elbow <- function(mat_to_clust, k = 10) {
  nn <- dbscan::kNNdist(sf::st_coordinates(mat_to_clust), k = k)
  nn_vec <- sort(nn)
  plot(nn_vec, type = "l", log = "y",
       xlab = "Point index (sorted)", ylab = sprintf("%d‑NN distance (µm)", k))
  abline(h = median(nn_vec), col = "red", lty = 2)
}



