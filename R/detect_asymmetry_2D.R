#' Title of the function
#'
#' @param file_path Path to input file
#' @param output_base Folder to store output
#' @return A list containing results
#' @export

detect_asymmetry_2D <- function(file_path, output_base = "asymmetry_output") {

  # Step 0: Create output directory
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M")
  output_dir <- paste0(output_base, "_", timestamp)
  dir.create(output_dir, recursive = TRUE)
  message("Output will be saved in: ", output_dir)

  # Step 1: Read and validate input
  data <- utils::read.csv(file_path)
  required_cols <- c("Specimen", "Landmark", "Side", "X", "Y")
  if (!all(required_cols %in% colnames(data))) {
    stop("Input CSV must contain columns: Specimen, Landmark, Side, X, Y")
  }
  if (!all(data$Side %in% c("Left", "Right"))) {
    stop("'Side' column must contain only 'Left' or 'Right'")
  }

  # Step 2: Preprocess
  data <- dplyr::mutate(data,
                        Landmark = as.factor(Landmark),
                        Side = factor(Side, levels = c("Left", "Right")),
                        Specimen = as.factor(Specimen))

  n_landmarks <- length(unique(data$Landmark))
  n_specimens <- length(unique(data$Specimen))
  specimen_ids <- unique(data$Specimen)

  # Step 3: Reshape
  wide_data <- data |>
    tidyr::pivot_wider(names_from = Side, values_from = c(X, Y)) |>
    tidyr::drop_na() |>
    dplyr::arrange(Specimen, Landmark)

  # Step 4: Build array
  array_data <- array(NA, dim = c(n_landmarks, 2, n_specimens * 2))
  config_names <- c()
  for (i in seq_along(specimen_ids)) {
    spec <- specimen_ids[i]
    spec_data <- dplyr::filter(wide_data, Specimen == spec)
    array_data[,,(2 * i - 1)] <- as.matrix(spec_data[, c("X_Left", "Y_Left")])
    array_data[,,(2 * i)]     <- as.matrix(spec_data[, c("X_Right", "Y_Right")])
    config_names <- c(config_names, paste0(spec, "_L"), paste0(spec, "_R"))
  }
  dimnames(array_data) <- list(levels(data$Landmark), c("X", "Y"), config_names)

  # Step 5: Mirror right side
  for (i in seq(2, dim(array_data)[3], by = 2)) {
    array_data[,1,i] <- -array_data[,1,i]  # Mirror X
  }

  # Step 6: GPA
  gpa <- geomorph::gpagen(array_data, print.progress = FALSE)

  # Step 7: Procrustes distance
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
    warning("No Procrustes ANOVA result generated.")
  }

  # Step 9: PCA Plot
  pca_coords <- geomorph::two.d.array(gpa$coords)
  pca <- stats::prcomp(pca_coords)
  pca_df <- data.frame(pca$x[, 1:2], Config = dimnames(gpa$coords)[[3]])

  gg_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, label = Config)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(nudge_y = 0.02, size = 3) +
    ggplot2::labs(title = "PCA of Aligned Configurations") +
    ggplot2::theme_minimal()

  ggplot2::ggsave(file.path(output_dir, "pca_plot.png"), gg_pca, width = 6, height = 5)

  # Step 10: Asymmetry Vector Plot
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
    ggplot2::labs(title = paste0("Asymmetry Vectors: Specimen ", specimen_ids[base_spec])) +
    ggplot2::theme_minimal()

  ggplot2::ggsave(file.path(output_dir, "asymmetry_vectors.png"), gg_vec, width = 6, height = 5)

  # Step 11: Heatmap
  mean_pd <- rowMeans(do.call(cbind, asymmetry_list))
  heatmap_df <- data.frame(Landmark = names(mean_pd), Distance = mean_pd)
  heatmap_df$Landmark <- factor(heatmap_df$Landmark, levels = heatmap_df$Landmark)

  gg_heat <- ggplot2::ggplot(heatmap_df, ggplot2::aes(x = Landmark, y = 1, fill = Distance)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient(low = "lightyellow", high = "red") +
    ggplot2::labs(title = "Heatmap of Asymmetry Magnitude", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())

  ggplot2::ggsave(file.path(output_dir, "asymmetry_heatmap.png"), gg_heat, width = 7, height = 3)

  # Step 12: TPS plot
  grDevices::png(file.path(output_dir, "tps_plot.png"), width = 600, height = 600)
  geomorph::plotRefToTarget(left, right, method = "TPS", mag = 2,
                            gridPars = geomorph::gridPar(grid.col = "gray"))
  graphics::title("Thin-Plate Spline (TPS) Deformation Plot")
  grDevices::dev.off()

  message("2D asymmetry analysis complete. Results saved in: ", output_dir)

  return(list(
    OutputFolder = output_dir,
    GPA = gpa,
    PerSpecimenAsymmetry = mean_asymmetry_df,
    Procrustes_ANOVA = symm
  ))
}
