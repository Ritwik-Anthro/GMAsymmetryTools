#' Slide Semilandmarks for 2D Data
#'
#' Performs Generalized Procrustes Analysis (GPA) with sliding semilandmarks on 2D landmark data.
#'
#' @param csv_file Path to the CSV file containing 2D landmark data. Required columns: Specimen, Landmark, Midline, Region, Side, Slider, CurveID, X, Y.
#'
#' @return A data frame with slid 2D coordinates and associated metadata.
#' @export

slide_semilandmarks_2D <- function(csv_file) {
  # Check for required packages
  if (!requireNamespace("geomorph", quietly = TRUE)) stop("Please install 'geomorph'")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install 'tidyr'")

  # Read and check data
  data <- utils::read.csv(csv_file)
  required_cols <- c("Specimen", "Landmark", "Midline", "Region", "Side", "Slider", "CurveID", "X", "Y")
  if (!all(required_cols %in% colnames(data))) {
    stop("CSV must contain columns: Specimen, Landmark, Midline, Region, Side, Slider, CurveID, X, Y")
  }

  # Add configuration ID
  data$Config <- paste(data$Specimen, data$Side, sep = "_")

  # Extract unique landmarks and configs
  landmarks <- unique(data$Landmark)
  configs <- unique(data$Config)
  n_lm <- length(landmarks)
  n_config <- length(configs)

  # Create array (landmarks × dimensions × configs)
  coords_array <- array(NA, dim = c(n_lm, 2, n_config),
                        dimnames = list(landmarks, c("X", "Y"), configs))

  for (i in seq_along(configs)) {
    config_id <- configs[i]
    subset <- dplyr::filter(data, Config == config_id)
    coords_array[ , , i] <- as.matrix(subset[match(landmarks, subset$Landmark), c("X", "Y")])
  }

  # Use reference config for sliding info
  ref_config <- dplyr::filter(data, Config == configs[1])
  slider_info <- ref_config[match(landmarks, ref_config$Landmark), c("Slider", "CurveID")]

  # Build curve matrix
  curve_matrix <- NULL
  unique_curves <- unique(slider_info$CurveID[slider_info$Slider == 1])
  curve_id <- 1

  for (cid in unique_curves) {
    indices <- which(slider_info$Slider == 1 & slider_info$CurveID == cid)
    if (length(indices) >= 3) {
      n_pairs <- length(indices) - 1
      new_matrix <- matrix(NA, nrow = n_pairs, ncol = 3)
      new_matrix[, 1] <- curve_id
      new_matrix[, 2] <- indices[1:n_pairs]
      new_matrix[, 3] <- indices[2:(n_pairs + 1)]
      curve_matrix <- rbind(curve_matrix, new_matrix)
      curve_id <- curve_id + 1
    }
  }

  # Perform GPA with sliding semilandmarks
  gpa <- geomorph::gpagen(coords_array,
                          curves = curve_matrix,
                          ProcD = TRUE,
                          print.progress = FALSE)

  # Reconstruct full metadata from original data
  metadata_fields <- c("Specimen", "Landmark", "Midline", "Region", "Side", "Slider", "CurveID")
  result_df <- NULL
  for (i in seq_along(configs)) {
    coords <- gpa$coords[ , , i]
    original_subset <- dplyr::filter(data, Config == configs[i])
    meta <- original_subset[match(landmarks, original_subset$Landmark), metadata_fields]
    df <- data.frame(
      meta,
      X = coords[,1],
      Y = coords[,2]
    )
    result_df <- rbind(result_df, df)
  }

  message("Sliding semilandmarks (2D) completed with metadata preserved.")
  return(result_df)
}
