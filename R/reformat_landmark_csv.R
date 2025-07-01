#' Reformat Landmark CSV for Asymmetry Analysis
#'
#' Reformats a raw landmark dataset by checking for required columns, sorting landmarks naturally, and saving the cleaned file.
#'
#' @param input_csv Path to the input CSV file.
#' @param output_csv Optional. Name of the output CSV file. Default is "reformatted_output.csv".
#'
#' @return None. Writes a reformatted CSV to disk.
#' @export

reformat_landmark_csv <- function(input_csv, output_csv = "reformatted_output.csv") {
  # Check if gtools is installed
  if (!requireNamespace("gtools", quietly = TRUE)) {
    stop("Package 'gtools' is required but not installed. Please install it using install.packages('gtools').")
  }

  # Read input CSV
  data <- utils::read.csv(input_csv, stringsAsFactors = FALSE)

  # Define required column sets without Slider and CurveID
  required_cols_2D <- c("Specimen", "Landmark", "Midline", "Region", "Side", "X", "Y")
  required_cols_3D <- c(required_cols_2D, "Z")

  # Determine if file is 2D or 3D
  is_3D <- all(required_cols_3D %in% names(data))
  is_2D <- all(required_cols_2D %in% names(data))

  if (!is_2D && !is_3D) {
    stop("CSV must contain either 2D (X, Y) or 3D (X, Y, Z) landmark fields along with metadata columns.")
  }

  # Keep only relevant columns
  data <- if (is_3D) {
    data[, required_cols_3D]
  } else {
    data[, required_cols_2D]
  }

  # Natural sort Landmark labels
  data$Landmark <- factor(
    data$Landmark,
    levels = gtools::mixedsort(unique(data$Landmark)),
    ordered = TRUE
  )

  # Sort by Specimen, Landmark, and Side
  data <- data[order(data$Specimen, data$Landmark, data$Side), ]

  # Write to output CSV
  utils::write.csv(data, file = output_csv, row.names = FALSE)
  message("Reformatted file saved as: ", normalizePath(output_csv))
}
