#' Title of the function
#'
#' @param file_path Path to input file
#' @param output_base Folder to store output
#' @return A list containing results
#' @export

detect_asymmetry_3D <- function(file_path, output_base = "asymmetry_output") {

  # Step 0: Create output directory
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M")
  output_dir <- paste0(output_base, "_", timestamp)
  dir.create(output_dir, recursive = TRUE)
  message("Output will be saved in: ", output_dir)

  # Step 1: Read and validate input
  data <- utils::read.csv(file_path)
  required_cols <- c("Specimen", "Landmark", "Side", "X", "Y", "Z")
  if (!all(required_cols %in% colnames(data))) {
    stop("CSV must contain: Specimen, Landmark, Side, X, Y, Z")
  }
  if (!all(data$Side %in% c("Left", "Right"))) {
    stop("'Side' column must only contain 'Left' or 'Right'")
  }

  # Step 2: Preprocess
  data <- dplyr::mutate(
    data,
    Landmark = as.factor(Landmark),
    Side = factor(Side, levels = c("Left", "Right")),
    Specimen = as.factor(Specimen)
  )

  n_landmarks <- length(unique(data$Landmark))
  n_specimens <- length(unique(data$Specimen))
  specimen_ids <- unique(data$Specimen)

  # Step 3: Wide format
  wide_data <- data |>
    tidyr::pivot_wider(names_from = Side, values_from = c(X, Y, Z)) |>
    tidyr::drop_na() |>
    dplyr::arrange(Specimen, Landmark)

  # Step 4: Create 3D array
  array_data <- array(NA, dim = c(n_landmarks, 3, n_specimens * 2))
  config_names <- c()
  for (i in seq_along(specimen_ids)) {
    spec <- specimen_ids[i]
    spec_data <- dplyr::filter(wide_data, Specimen == spec)
    array_data[,,(2 * i - 1)] <- as.matrix(spec_data[, c("X_Left", "Y_Left", "Z_Left")])
    array_data[,,(2 * i)]     <- as.matrix(spec_data[, c("X_Right", "Y_Right", "Z_Right")])
    config_names <- c(config_names, paste0(spec, "_L"), paste0(spec, "_R"))
  }
  dimnames(array_data) <- list(levels(data$Landmark), c("X", "Y", "Z"), config_names)

  # Step 5: Mirror right side (X only)
  for (i in seq(2, dim(array_data)[3], by = 2)) {
    array_data[,1,i] <- -array_data[,1,i]
  }

  # Step 6: GPA
  gpa <- geomorph::gpagen(array_data, print.progress = FALSE)

  # Save mean shape
  mean_shape <- geomorph::mshape(gpa$coords)
  mean_shape_df <- as.data.frame(mean_shape)
  mean_shape_df$Landmark <- rownames(mean_shape_df)
  mean_shape_df <- mean_shape_df[, c("Landmark", "X", "Y", "Z")]
  utils::write.csv(mean_shape_df, file.path(output_dir, "mean_shape.csv"), row.names = FALSE)

  # Step 7: Per-specimen Procrustes distance
  asymmetry_list <- list()
  for (i in seq(1, dim(gpa$coords)[3], by = 2)) {
    left <- gpa$coords[,,i]
    right <- gpa$coords[,,i+1]
    pd <- sqrt(rowSums((left - right)^2))
    asymmetry_list[[specimen_ids[ceiling(i/2)]]] <- pd
  }

  mean_asymmetry_df <- data.frame(
    Specimen = specimen_ids,
    Mean_Procrustes_Distance = sapply(asymmetry_list, mean)
  )
  utils::write.csv(mean_asymmetry_df, file.path(output_dir, "mean_procrustes_distances.csv"), row.names = FALSE)

  # Step 8: Procrustes ANOVA
  ind_factors <- rep(specimen_ids, each = 2)
  side_factors <- rep(c("Left", "Right"), times = n_specimens)
  symm <- geomorph::bilat.symmetry(gpa$coords,
                                   ind = factor(ind_factors),
                                   side = factor(side_factors),
                                   land.pairs = NULL,
                                   iter = 999)
  if (!is.null(symm)) {
    sink(file.path(output_dir, "procrustes_anova_summary.txt"))
    print(summary(symm))
    sink()
  } else {
    warning("Procrustes ANOVA result is empty or malformed. No summary file saved.")
  }

  # Step 9: PCA plot
  pca_coords <- geomorph::two.d.array(gpa$coords)
  pca <- stats::prcomp(pca_coords)
  pca_df <- data.frame(pca$x[, 1:2], Config = dimnames(gpa$coords)[[3]])

  gg_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, label = Config)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(nudge_y = 0.02, size = 3) +
    ggplot2::labs(title = "PCA (2D projection) of 3D Aligned Configurations") +
    ggplot2::theme_minimal()

  ggplot2::ggsave(file.path(output_dir, "pca_plot_3D_projection.png"), gg_pca, width = 6, height = 5)

  # Step 10: Displacement Arrows (on 3D plot)
  base_spec <- 1
  left <- gpa$coords[,,(2 * base_spec - 1)]
  right <- gpa$coords[,,(2 * base_spec)]

  rgl::open3d()
  rgl::plot3d(left, type = "s", col = "blue", radius = 0.01)
  rgl::points3d(right, col = "red", radius = 0.01)
  for (j in 1:n_landmarks) {
    rgl::segments3d(rbind(left[j,], right[j,]), col = "green")
  }
  htmlwidgets::saveWidget(rgl::rglwidget(), file = file.path(output_dir, "landmark_displacement_arrows.html"), selfcontained = TRUE)

  # Step 11: Residual Vectors to Mean Shape
  rgl::open3d()
  rgl::plot3d(mean_shape, type = "s", col = "gray", radius = 0.01)
  for (j in 1:n_landmarks) {
    rgl::segments3d(rbind(mean_shape[j,], gpa$coords[j,,1]), col = "purple")
  }
  htmlwidgets::saveWidget(rgl::rglwidget(), file = file.path(output_dir, "procrustes_residual_vectors.html"), selfcontained = TRUE)

  # Final message
  message("3D asymmetry analysis complete. Files saved in: ", output_dir)

  return(list(
    OutputFolder = output_dir,
    GPA = gpa,
    MeanShape = mean_shape_df,
    PerSpecimenAsymmetry = mean_asymmetry_df,
    Procrustes_ANOVA = symm
  ))
}
