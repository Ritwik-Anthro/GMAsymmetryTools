#' Detect Bilateral Asymmetry in 3D Landmark Data
#'
#' This function performs region-wise bilateral asymmetry analysis on 3D landmark data.
#' It includes Procrustes superimposition, PCA, Mahalanobis distances, Procrustes ANOVA,
#' and visualization of asymmetry vectors and residuals. Outputs are saved in a timestamped
#' folder with region-specific subfolders.
#'
#' @param file_path Path to the input CSV file containing 3D landmark data. Required columns:
#'   \code{Specimen}, \code{Landmark}, \code{Region}, \code{Midline}, \code{Side}, \code{X}, \code{Y}, \code{Z}.
#' @param output_base Optional base name for the output directory. Default is \code{"asymmetry_output"}.
#'
#' @return A list with:
#' \describe{
#'   \item{OutputFolder}{Path to the output directory where results are saved}
#'   \item{RegionResults}{A named list for each region containing:
#'     \code{MeanShape}, \code{GPA}, \code{Procrustes_ANOVA}, \code{PerSpecimenAsymmetry},
#'     \code{Mahalanobis}, and \code{PCA} objects}
#' }
#'
#' @details
#' This function mirrors only the non-midline landmarks along the X-axis for the right side.
#' It automatically separates landmarks by region and performs shape analysis and asymmetry
#' quantification within each region. Interactive 3D plots are saved as HTML widgets.
#'
#' @importFrom utils read.csv write.csv
#' @importFrom stats prcomp cov
#' @importFrom MASS ginv
#' @importFrom tidyr pivot_wider drop_na
#' @importFrom ggplot2 ggplot aes geom_point geom_text labs theme_minimal ggsave
#' @importFrom rgl open3d plot3d points3d segments3d rglwidget
#' @importFrom htmlwidgets saveWidget
#' @importFrom geomorph gpagen mshape bilat.symmetry two.d.array
#'
#' @export

detect_asymmetry_3D <- function(file_path, output_base = "asymmetry_output") {
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M")
  output_dir <- paste0(output_base, "_", timestamp)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Output will be saved in: ", output_dir)

  data <- utils::read.csv(file_path)
  required_cols <- c("Specimen", "Landmark", "Region", "Midline", "Side", "X", "Y", "Z")
  if (!all(required_cols %in% colnames(data))) {
    stop("CSV must include columns: Specimen, Landmark, Region, Midline, Side, X, Y, Z")
  }

  data <- data.frame(data)
  data$Specimen <- as.factor(data$Specimen)
  data$Landmark <- as.factor(data$Landmark)
  data$Region <- as.factor(data$Region)
  data$Side <- factor(data$Side, levels = c("Left", "Right"))
  data$Midline <- as.numeric(data$Midline)

  regions <- unique(data$Region)
  region_results <- list()

  for (reg in regions) {
    message("Processing region: ", reg)
    region_data <- data[data$Region == reg, ]

    wide_data <- tidyr::pivot_wider(region_data, names_from = Side, values_from = c(X, Y, Z))
    wide_data <- tidyr::drop_na(wide_data)
    wide_data <- wide_data[order(wide_data$Specimen, wide_data$Landmark), ]

    specimens <- unique(wide_data$Specimen)
    n_landmarks <- length(unique(wide_data$Landmark))
    n_specimens <- length(specimens)

    if (n_specimens < 1 || n_landmarks < 1) {
      warning("Skipping region ", reg, ": insufficient data.")
      next
    }

    array_data <- array(NA, dim = c(n_landmarks, 3, n_specimens * 2))
    config_names <- c()

    for (i in seq_along(specimens)) {
      spec <- specimens[i]
      spec_data <- wide_data[wide_data$Specimen == spec, ]
      array_data[, , (2 * i - 1)] <- as.matrix(spec_data[, c("X_Left", "Y_Left", "Z_Left")])
      array_data[, , (2 * i)]     <- as.matrix(spec_data[, c("X_Right", "Y_Right", "Z_Right")])
      config_names <- c(config_names, paste0(spec, "_L"), paste0(spec, "_R"))
    }

    landmark_labels <- sort(unique(wide_data$Landmark))
    dimnames(array_data) <- list(landmark_labels, c("X", "Y", "Z"), config_names)

    is_midline <- region_data[!duplicated(region_data$Landmark), c("Landmark", "Midline")]
    is_midline <- is_midline[match(landmark_labels, is_midline$Landmark), "Midline"]

    for (i in seq(2, dim(array_data)[3], by = 2)) {
      array_data[which(is_midline == 0), 1, i] <- -array_data[which(is_midline == 0), 1, i]
    }

    gpa <- geomorph::gpagen(array_data, print.progress = FALSE)

    region_dir <- file.path(output_dir, paste0("Region_", reg))
    dir.create(region_dir, showWarnings = FALSE)

    mean_shape <- geomorph::mshape(gpa$coords)
    mean_shape_df <- as.data.frame(mean_shape)
    mean_shape_df$Landmark <- rownames(mean_shape_df)
    mean_shape_df <- mean_shape_df[, c("Landmark", "X", "Y", "Z")]
    utils::write.csv(mean_shape_df, file.path(region_dir, "mean_shape.csv"), row.names = FALSE)

    centroid_df <- data.frame(Specimen_Config = names(gpa$Csize), Centroid_Size = gpa$Csize)
    utils::write.csv(centroid_df, file.path(region_dir, "centroid_sizes.csv"), row.names = FALSE)

    # Procrustes distance
    asymmetry_list <- list()
    for (i in seq(1, dim(gpa$coords)[3], by = 2)) {
      left <- gpa$coords[, , i]
      right <- gpa$coords[, , i + 1]
      pd <- sqrt(rowSums((left - right)^2))
      asymmetry_list[[specimens[ceiling(i / 2)]]] <- pd
    }

    mean_asym_df <- data.frame(
      Specimen = specimens,
      Mean_Procrustes_Distance = sapply(asymmetry_list, mean)
    )
    utils::write.csv(mean_asym_df, file.path(region_dir, "mean_procrustes_distances.csv"), row.names = FALSE)

    # Mahalanobis distance
    maha_df <- data.frame(Specimen = character(), Mahalanobis_Distance = numeric())
    for (i in seq(1, dim(gpa$coords)[3], by = 2)) {
      left <- as.vector(gpa$coords[, , i])
      right <- as.vector(gpa$coords[, , i + 1])
      pooled_cov <- stats::cov(rbind(left, right))
      inv_cov <- tryCatch(solve(pooled_cov), error = function(e) MASS::ginv(pooled_cov))
      d2 <- sqrt(t(left - right) %*% inv_cov %*% (left - right))
      maha_df <- rbind(maha_df, data.frame(Specimen = specimens[ceiling(i / 2)], Mahalanobis_Distance = d2))
    }
    utils::write.csv(maha_df, file.path(region_dir, "mahalanobis_distances.csv"), row.names = FALSE)

    # Procrustes ANOVA
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

    # PCA
    pca_coords <- geomorph::two.d.array(gpa$coords)
    pca <- stats::prcomp(pca_coords)
    pca_all_df <- as.data.frame(pca$x)
    pca_all_df$Config <- rownames(pca_all_df)
    utils::write.csv(pca_all_df, file.path(region_dir, "pca_scores.csv"), row.names = FALSE)

    variance <- pca$sdev^2
    percent_var <- 100 * variance / sum(variance)
    cumulative_var <- cumsum(percent_var)
    variance_df <- data.frame(
      PC = paste0("PC", seq_along(percent_var)),
      Percent_Variance = percent_var,
      Cumulative_Variance = cumulative_var
    )
    utils::write.csv(variance_df, file.path(region_dir, "pca_variance_explained.csv"), row.names = FALSE)

    # PCA plot
    gg_pca <- ggplot2::ggplot(pca_all_df, ggplot2::aes(x = PC1, y = PC2, label = Config)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_text(nudge_y = 0.02, size = 3) +
      ggplot2::labs(title = paste("PCA: Region", reg)) +
      ggplot2::theme_minimal()
    ggplot2::ggsave(file.path(region_dir, "pca_plot.png"), gg_pca, width = 6, height = 5)

    # 3D arrows
    base_spec <- 1
    left <- gpa$coords[, , (2 * base_spec - 1)]
    right <- gpa$coords[, , (2 * base_spec)]

    rgl::open3d()
    rgl::plot3d(left, type = "s", col = "blue", radius = 0.01)
    rgl::points3d(right, col = "red", radius = 0.01)
    for (j in 1:n_landmarks) {
      rgl::segments3d(rbind(left[j, ], right[j, ]), col = "green")
    }
    htmlwidgets::saveWidget(rgl::rglwidget(), file = file.path(region_dir, "landmark_displacement_arrows.html"), selfcontained = TRUE)

    rgl::open3d()
    rgl::plot3d(mean_shape, type = "s", col = "gray", radius = 0.01)
    for (j in 1:n_landmarks) {
      rgl::segments3d(rbind(mean_shape[j, ], gpa$coords[j, , 1]), col = "purple")
    }
    htmlwidgets::saveWidget(rgl::rglwidget(), file = file.path(region_dir, "procrustes_residual_vectors.html"), selfcontained = TRUE)

    region_results[[paste0("Region_", reg)]] <- list(
      MeanShape = mean_shape_df,
      GPA = gpa,
      Procrustes_ANOVA = symm,
      PerSpecimenAsymmetry = mean_asym_df,
      Mahalanobis = maha_df,
      PCA = pca
    )
  }

  message("Region-wise 3D asymmetry analysis complete: ", output_dir)

  return(list(
    OutputFolder = output_dir,
    RegionResults = region_results
  ))
}
