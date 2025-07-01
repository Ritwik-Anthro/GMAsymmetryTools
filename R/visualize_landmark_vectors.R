#' Visualize 3D Landmark Displacement Vectors
#'
#' Creates an interactive HTML plot showing displacement vectors between two sets of 3D landmarks.
#'
#' @param source_lm_file CSV file containing source landmarks (columns: X, Y, Z).
#' @param target_lm_file CSV file containing target landmarks (columns: X, Y, Z).
#' @param output_file Path to save the interactive HTML file. Default is "landmark_vectors.html".
#'
#' @return None. Saves an interactive visualization as HTML.
#' @export

visualize_landmark_vectors <- function(source_lm_file, target_lm_file, output_file = "landmark_vectors.html") {
  # Check for required packages
  if (!requireNamespace("rgl", quietly = TRUE)) stop("Package 'rgl' is required but not installed.")
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) stop("Package 'htmlwidgets' is required but not installed.")

  # Read input CSVs
  source <- utils::read.csv(source_lm_file)
  target <- utils::read.csv(target_lm_file)

  # Sanity checks
  if (nrow(source) == 0 || nrow(target) == 0) stop("One or both landmark files are empty.")
  required_cols <- c("X", "Y", "Z")
  if (!all(required_cols %in% names(source)) || !all(required_cols %in% names(target))) {
    stop("Landmark files must contain columns: X, Y, Z")
  }
  if (nrow(source) != nrow(target)) stop("Source and target files must have the same number of rows.")

  # Set up 3D window
  rgl::open3d()
  rgl::par3d(windowRect = c(100, 100, 700, 700))

  # Plot landmarks
  rgl::points3d(source$X, source$Y, source$Z, col = "black", size = 8)
  rgl::points3d(target$X, target$Y, target$Z, col = "red", size = 8)

  # Draw vectors between corresponding points
  for (i in seq_len(nrow(source))) {
    rgl::segments3d(
      x = c(source[i, "X"], target[i, "X"]),
      y = c(source[i, "Y"], target[i, "Y"]),
      z = c(source[i, "Z"], target[i, "Z"]),
      col = "blue", lwd = 2
    )
  }

  # Export interactive widget
  widget <- rgl::rglwidget()
  htmlwidgets::saveWidget(widget, file = output_file, selfcontained = TRUE)
  message("Vector visualization saved to: ", output_file)
}
