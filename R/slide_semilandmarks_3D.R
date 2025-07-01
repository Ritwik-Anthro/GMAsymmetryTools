#' Slide Semilandmarks for 3D Data
#'
#' Performs Generalized Procrustes Analysis (GPA) with sliding semilandmarks on 3D landmark data.
#'
#' @param csv_file Path to the CSV file containing 3D landmark data. Required columns: Specimen, Landmark, Midline, Region, Side, Slider, CurveID, X, Y, Z.
#'
#' @return A data frame with slid 3D coordinates and associated metadata.
#' @export

slide_semilandmarks_3D <- function(csv_file) {
  # Ensure required packages are installed
  if (!requireNamespace("geomorph", quietly = TRUE)) stop("Please install the 'geomorph' package.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install the 'dplyr' package.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install the 'tidyr' package.")

  # Read and check data
  data <- utils::read.csv(csv_file)
  required_cols <- c("Specimen", "Landmark", "Midline", "Region", "Side", "Slider", "CurveID", "X", "Y", "Z")
  if (!all(required_cols %in% colnames(data))) {
    stop("CSV must contain columns: Specimen, Landmark, Midline, Region, Side, Slider, CurveID, X, Y, Z")
  }

  # Create configuration ID
  data$Config <- paste(data$Specimen, data$Side, sep = "_")

  # Get unique landmarks and configurations
  landmarks <- unique(data$Landmark)
  configs <- unique(data$Config)
  n_lm <- length(landmarks)
  n_config <- length(configs)

  # Build 3D coordinate array
  coords_array <- array(NA, dim = c(n_lm, 3, n_config),
                        dimnames = list(landmarks, c("X", "Y", "Z"), configs))

  for (i in seq_along(configs)) {
    config_id <- configs[i]
    subset <- dplyr::filter(data, Config == config_id)
    coords_array[ , , i] <- as.matrix(subset[match(landmarks, subset$Landmark), c("X", "Y", "Z")])
  }

  # Extract sliding semilandmark info from reference configuration
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

  # Run GPA with sliding semilandmarks
  gpa <- geomorph::gpagen(coords_array,
                          curves = curve_matrix,
                          ProcD = TRUE,
                          print.progress = FALSE)

  # Preserve metadata in output
  metadata_fields <- c("Specimen", "Landmark", "Midline", "Region", "Side", "Slider", "CurveID")
  result_df <- NULL

  for (i in seq_along(configs)) {
    coords <- gpa$coords[ , , i]
    original_subset <- dplyr::filter(data, Config == configs[i])
    meta <- original_subset[match(landmarks, original_subset$Landmark), metadata_fields]

    df <- data.frame(
      meta,
      X = coords[, 1],
      Y = coords[, 2],
      Z = coords[, 3]
    )
    result_df <- rbind(result_df, df)
  }

  message("Sliding semilandmarks (3D) completed. All metadata retained.")
  return(result_df)
}
