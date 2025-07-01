#' Warp Mesh Using Thin Plate Spline (TPS)
#'
#' Warps a 3D mesh using TPS from source to target landmarks and saves static and interactive outputs.
#'
#' @param source_lm_file CSV file with source landmarks (columns: X, Y, Z).
#' @param target_lm_file CSV file with target landmarks (columns: X, Y, Z).
#' @param mesh_file Path to the 3D mesh file (.ply or .stl format).
#' @param output_base Optional. Base name for the output directory. Default is "mesh_warp_output".
#' @param views Optional. Vector of view names (e.g., "frontal", "lateral") for image snapshots.
#' @param warp_magnitude Optional. Magnification factor for TPS deformation. Default is 1.
#'
#' @return A list containing the output directory path and warped mesh object.
#' @export

warp_mesh_TPS <- function(source_lm_file, target_lm_file, mesh_file,
                          output_base = "mesh_warp_output",
                          views = c("frontal", "lateral"),
                          warp_magnitude = 1) {
  # Ensure required packages are available
  if (!requireNamespace("Morpho", quietly = TRUE)) stop("Please install the 'Morpho' package.")
  if (!requireNamespace("Rvcg", quietly = TRUE)) stop("Please install the 'Rvcg' package.")
  if (!requireNamespace("rgl", quietly = TRUE)) stop("Please install the 'rgl' package.")
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) stop("Please install the 'htmlwidgets' package.")
  if (!requireNamespace("tools", quietly = TRUE)) stop("Please install the 'tools' package.")

  # Create output directory
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M")
  output_dir <- paste0(output_base, "_", timestamp)
  dir.create(output_dir, recursive = TRUE)
  message("Output will be saved in: ", output_dir)

  # Load landmark files
  source_lm_df <- utils::read.csv(source_lm_file, header = TRUE)
  target_lm_df <- utils::read.csv(target_lm_file, header = TRUE)

  # Check columns
  required_cols <- c("X", "Y", "Z")
  if (!all(required_cols %in% names(source_lm_df)) || !all(required_cols %in% names(target_lm_df))) {
    stop("Landmark files must contain columns: X, Y, Z")
  }

  # Extract coordinate matrices
  source_lm <- as.matrix(source_lm_df[, required_cols])
  target_lm <- as.matrix(target_lm_df[, required_cols])

  # Check dimensions match
  if (!all(dim(source_lm) == dim(target_lm))) {
    stop("Source and target landmarks must have the same number and dimension of landmarks.")
  }

  # Apply warp magnitude
  delta <- (target_lm - source_lm) * warp_magnitude
  adjusted_target <- source_lm + delta

  # Read mesh
  mesh_ext <- tools::file_ext(mesh_file)
  mesh <- switch(
    mesh_ext,
    ply = Rvcg::vcgPlyRead(mesh_file, updateNormals = TRUE),
    stl = Rvcg::vcgImport(mesh_file),
    stop("Unsupported mesh format. Please use a .ply or .stl file.")
  )

  message("Warping mesh using TPS with warp_magnitude = ", warp_magnitude)
  warped_mesh <- Morpho::tps3d(x = mesh, refmat = source_lm, tarmat = adjusted_target)

  # Save warped mesh
  ply_out <- file.path(output_dir, "warped_mesh.ply")
  Rvcg::vcgPlyWrite(warped_mesh, filename = ply_out, binary = TRUE)
  message("Warped mesh saved to: ", ply_out)

  # Render warped mesh and landmarks
  rgl::open3d()
  rgl::shade3d(warped_mesh, color = "lightblue", alpha = 1)
  rgl::points3d(adjusted_target, col = "red", size = 5)

  # Save static images
  for (view in views) {
    view_file <- file.path(output_dir, paste0("mesh_view_", view, ".png"))
    if (view == "frontal") {
      rgl::view3d(userMatrix = rgl::rotationMatrix(0, 0, 1, 0))  # front
    } else if (view == "lateral") {
      rgl::view3d(userMatrix = rgl::rotationMatrix(pi / 2, 0, 1, 0))  # side
    }
    rgl::rgl.snapshot(view_file)
    message("Saved ", view, " view to: ", view_file)
  }

  # Save interactive viewer
  html_file <- file.path(output_dir, "warped_mesh_viewer.html")
  widget <- rgl::rglwidget()
  htmlwidgets::saveWidget(widget, file = html_file, selfcontained = TRUE)
  message("Interactive HTML viewer saved to: ", html_file)

  message("TPS mesh warping complete.")

  return(list(
    OutputDirectory = output_dir,
    WarpedMesh = warped_mesh
  ))
}
