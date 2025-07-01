#' Detect Bilateral Asymmetry in 2D Landmark Data (Region-wise)
#'
#' This function analyzes 2D landmark coordinate data for bilateral asymmetry on a per-region basis.
#' It performs Generalized Procrustes Analysis (GPA), calculates Procrustes distances,
#' performs Procrustes ANOVA, and generates PCA plots, Mahalanobis distances, TPS deformation plots,
#' asymmetry vectors, and heatmaps. All outputs are saved to a timestamped folder.
#'
#' @param file_path Character. Path to the input CSV file containing 2D landmark data.
#'   The file must include the columns: \code{Specimen}, \code{Landmark}, \code{Midline},
#'   \code{Region}, \code{Side}, \code{X}, and \code{Y}.
#' @param output_base Character. Optional base name for the output directory.
#'   A timestamp will be appended to create a unique folder. Default is \code{"asymmetry_output"}.
#' @param tps_mag Numeric. Magnification factor for the thin-plate spline (TPS) deformation grid
#'   in the TPS plot. Default is \code{2}.
#'
#' @return (Invisibly) the full path to the output directory. Also writes multiple output files per region:
#' \itemize{
#'   \item \code{centroid_sizes.csv} - Centroid size per configuration.
#'   \item \code{mean_procrustes_distances.csv} - Mean shape difference between left and right sides.
#'   \item \code{mahalanobis_distances.csv} - Mahalanobis distances per specimen.
#'   \item \code{procrustes_anova_summary.txt} - Procrustes ANOVA result.
#'   \item \code{pca_scores.csv} - Principal component scores for each configuration.
#'   \item \code{pca_variance_explained.csv} - Percentage and cumulative variance explained by PCs.
#'   \item \code{pca_plot.png} - 2D PCA scatterplot.
#'   \item \code{asymmetry_vectors.png} - Landmark displacement vectors (Left vs Right).
#'   \item \code{tps_plot.png} - Thin-plate spline deformation plot.
#'   \item \code{asymmetry_heatmap.png} - Heatmap showing average landmark-wise asymmetry.
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom utils read.csv write.csv
#' @importFrom stats prcomp cov
#' @importFrom MASS ginv
#' @importFrom geomorph gpagen plotRefToTarget gridPar bilat.symmetry two.d.array
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_segment geom_tile
#' @importFrom ggplot2 labs theme_minimal ggsave scale_fill_gradient theme element_blank
#' @importFrom grid arrow unit
#' @importFrom grDevices png dev.off
#' @importFrom graphics title
#'
#' @export


detect_asymmetry_2D <- function(file_path, output_base = "asymmetry_output", tps_mag = 2) {
  `%>%` <- magrittr::`%>%`

  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M")
  output_dir <- paste0(output_base, "_", timestamp)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Output will be saved in: ", output_dir)

  data <- utils::read.csv(file_path)
  required_cols <- c("Specimen", "Landmark", "Midline", "Region", "Side", "X", "Y")
  if (!all(required_cols %in% names(data))) {
    stop("Missing required columns: ", paste(setdiff(required_cols, names(data)), collapse = ", "))
  }

  data <- dplyr::mutate(
    data,
    Specimen = as.factor(Specimen),
    Landmark = as.factor(Landmark),
    Region = as.factor(Region),
    Side = factor(Side, levels = c("Left", "Right")),
    Midline = as.numeric(Midline)
  )

  regions <- unique(data$Region)

  for (reg in regions) {
    message("Processing region: ", reg)
    region_data <- dplyr::filter(data, Region == reg)

    wide_data <- region_data %>%
      tidyr::pivot_wider(names_from = Side, values_from = c(X, Y)) %>%
      tidyr::drop_na() %>%
      dplyr::arrange(Specimen, Landmark)

    specimens <- unique(wide_data$Specimen)
    n_landmarks <- length(unique(wide_data$Landmark))
    n_specimens <- length(specimens)

    if (n_specimens < 1 || n_landmarks < 1) {
      warning("Skipping region ", reg, ": insufficient data.")
      next
    }

    array_data <- array(NA, dim = c(n_landmarks, 2, n_specimens * 2))
    config_names <- c()

    for (i in seq_along(specimens)) {
      spec <- specimens[i]
      spec_data <- dplyr::filter(wide_data, Specimen == spec)
      array_data[,,(2 * i - 1)] <- as.matrix(spec_data[, c("X_Left", "Y_Left")])
      array_data[,,(2 * i)]     <- as.matrix(spec_data[, c("X_Right", "Y_Right")])
      config_names <- c(config_names, paste0(spec, "_L"), paste0(spec, "_R"))
    }

    landmark_labels <- wide_data %>%
      dplyr::distinct(Landmark) %>%
      dplyr::pull(Landmark) %>%
      as.character()
    dimnames(array_data) <- list(landmark_labels, c("X", "Y"), config_names)

    is_midline <- region_data %>%
      dplyr::distinct(Landmark, Midline) %>%
      dplyr::arrange(match(Landmark, landmark_labels)) %>%
      dplyr::pull(Midline)

    for (i in seq(2, dim(array_data)[3], by = 2)) {
      array_data[which(is_midline == 0), 1, i] <- -array_data[which(is_midline == 0), 1, i]
    }

    gpa <- geomorph::gpagen(array_data, print.progress = FALSE)

    region_dir <- file.path(output_dir, paste0("Region_", reg))
    dir.create(region_dir, showWarnings = FALSE)

    centroid_df <- data.frame(Specimen_Config = names(gpa$Csize), Centroid_Size = gpa$Csize)
    utils::write.csv(centroid_df, file.path(region_dir, "centroid_sizes.csv"), row.names = FALSE)

    asymmetry_list <- list()
    for (i in seq(1, dim(gpa$coords)[3], by = 2)) {
      left <- gpa$coords[,,i]
      right <- gpa$coords[,,i+1]
      pd <- sqrt(rowSums((left - right)^2))
      asymmetry_list[[specimens[ceiling(i / 2)]]] <- pd
    }

    mean_asym_df <- data.frame(
      Specimen = specimens,
      Mean_Procrustes_Distance = sapply(asymmetry_list, mean)
    )
    utils::write.csv(mean_asym_df, file.path(region_dir, "mean_procrustes_distances.csv"), row.names = FALSE)

    ind_factors <- rep(specimens, each = 2)
    side_factors <- rep(c("Left", "Right"), times = n_specimens)
    symm <- geomorph::bilat.symmetry(gpa$coords,
                                     ind = factor(ind_factors),
                                     side = factor(side_factors),
                                     land.pairs = NULL,
                                     iter = 999)

    sink(file.path(region_dir, "procrustes_anova_summary.txt"))
    summary(symm)
    sink()

    pca_coords <- geomorph::two.d.array(gpa$coords)
    pca <- stats::prcomp(pca_coords)

    pca_all_df <- as.data.frame(pca$x)
    pca_all_df$Config <- rownames(pca_all_df)
    utils::write.csv(pca_all_df, file.path(region_dir, "pca_scores.csv"), row.names = FALSE)

    variance <- pca$sdev^2
    percent_var <- 100 * variance / sum(variance)
    cumulative_var <- cumsum(percent_var)
    variance_df <- data.frame(
      PC = paste0("PC", 1:length(percent_var)),
      Percent_Variance = percent_var,
      Cumulative_Variance = cumulative_var
    )
    utils::write.csv(variance_df, file.path(region_dir, "pca_variance_explained.csv"), row.names = FALSE)

    maha_df <- data.frame(Specimen = character(), Mahalanobis_Distance = numeric())
    for (i in seq(1, dim(gpa$coords)[3], by = 2)) {
      left <- as.vector(gpa$coords[,,i])
      right <- as.vector(gpa$coords[,,i + 1])
      pooled_cov <- stats::cov(rbind(left, right))
      inv_cov <- tryCatch(solve(pooled_cov), error = function(e) MASS::ginv(pooled_cov))
      d2 <- sqrt(t(left - right) %*% inv_cov %*% (left - right))
      maha_df <- rbind(maha_df, data.frame(Specimen = specimens[ceiling(i / 2)], Mahalanobis_Distance = d2))
    }
    utils::write.csv(maha_df, file.path(region_dir, "mahalanobis_distances.csv"), row.names = FALSE)

    gg_pca <- ggplot2::ggplot(pca_all_df, ggplot2::aes(x = PC1, y = PC2, label = Config)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_text(nudge_y = 0.02, size = 3) +
      ggplot2::labs(title = paste("PCA: Region", reg)) +
      ggplot2::theme_minimal()
    ggplot2::ggsave(file.path(region_dir, "pca_plot.png"), gg_pca, width = 6, height = 5)

    base_spec <- 1
    left <- gpa$coords[,,(2 * base_spec - 1)]
    right <- gpa$coords[,,(2 * base_spec)]
    vectors_df <- data.frame(
      X1 = left[,1], Y1 = left[,2],
      X2 = right[,1], Y2 = right[,2],
      Landmark = rownames(left)
    )

    gg_vec <- ggplot2::ggplot(vectors_df) +
      ggplot2::geom_segment(ggplot2::aes(x = X1, y = Y1, xend = X2, yend = Y2),
                            arrow = grid::arrow(length = grid::unit(0.1, "cm")), color = "blue") +
      ggplot2::geom_point(ggplot2::aes(x = X1, y = Y1), color = "black", size = 2) +
      ggplot2::geom_point(ggplot2::aes(x = X2, y = Y2), color = "red", size = 2) +
      ggplot2::labs(title = paste0("Asymmetry Vectors: ", specimens[base_spec])) +
      ggplot2::theme_minimal()
    ggplot2::ggsave(file.path(region_dir, "asymmetry_vectors.png"), gg_vec, width = 6, height = 5)

    grDevices::png(file.path(region_dir, "tps_plot.png"), width = 600, height = 600)
    geomorph::plotRefToTarget(left, right, method = "TPS", mag = tps_mag,
                              gridPars = geomorph::gridPar(grid.col = "gray"))
    graphics::title(paste("TPS Deformation Plot - Region", reg))
    grDevices::dev.off()

    mean_pd <- rowMeans(do.call(cbind, asymmetry_list))
    heatmap_df <- data.frame(Landmark = names(mean_pd), Distance = mean_pd)
    heatmap_df$Landmark <- factor(heatmap_df$Landmark, levels = heatmap_df$Landmark)

    gg_heat <- ggplot2::ggplot(heatmap_df, ggplot2::aes(x = Landmark, y = 1, fill = Distance)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient(low = "lightyellow", high = "red") +
      ggplot2::labs(title = paste("Asymmetry Heatmap: Region", reg), y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())
    ggplot2::ggsave(file.path(region_dir, "asymmetry_heatmap.png"), gg_heat, width = 7, height = 3)
  }

  message("Region-wise asymmetry analysis completed and saved to: ", output_dir)
  return(invisible(output_dir))
}
