#' Warp Mesh Using TPS
#'
#' This function warps a mesh from a source landmark configuration to a target configuration using thin-plate spline (TPS) transformation.
#'
#' @param source_lm_file Path to the CSV file with source (reference) landmarks. Must contain X, Y, and Z columns.
#' @param target_lm_file Path to the CSV file with target (specimen) landmarks. Must contain X, Y, and Z columns.
#' @param mesh_file Path to the mesh file (.ply or .stl) corresponding to the source.
#' @param output_base Base name for the output folder. A timestamp will be appended.
#' @param views Character vector of views to save as PNGs (default: c("frontal", "lateral")).
#'
#' @return A list containing the output directory and warped mesh object.
#' @export

warp_mesh_TPS <- function(source_lm_file, target_lm_file, mesh_file,
                          output_base = "mesh_warp_output",
                          views = c("frontal", "lateral")) {

  # Create timestamped output directory
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M")
  output_dir <- paste0(output_base, "_", timestamp)
  dir.create(output_dir, recursive = TRUE)
  message("Output will be saved in: ", output_dir)

  # Load landmarks
  source_lm_df <- utils::read.csv(source_lm_file, header = TRUE)
  target_lm_df <- utils::read.csv(target_lm_file, header = TRUE)

  if (all(c("X", "Y", "Z") %in% names(source_lm_df))) {
    source_lm <- as.matrix(source_lm_df[, c("X", "Y", "Z")])
  } else {
    stop("Source landmark file must contain columns: X, Y, Z")
  }

  if (all(c("X", "Y", "Z") %in% names(target_lm_df))) {
    target_lm <- as.matrix(target_lm_df[, c("X", "Y", "Z")])
  } else {
    stop("Target landmark file must contain columns: X, Y, Z")
  }

  # Dimension check
  if (!all(dim(source_lm) == dim(target_lm))) {
    stop("Source and target landmarks must have the same dimensions and number of landmarks.")
  }

  # Load mesh based on extension
  mesh_ext <- tools::file_ext(mesh_file)
  mesh <- switch(mesh_ext,
                 "ply" = Rvcg::vcgPlyRead(mesh_file, updateNormals = TRUE),
                 "stl" = Rvcg::vcgImport(mesh_file),
                 stop("Unsupported mesh format. Please use .ply or .stl"))

  # Warp the mesh using TPS
  message("Warping mesh using TPS...")
  warped_mesh <- Morpho::tps3d(x = mesh, refmat = source_lm, tarmat = target_lm)

  # Save warped mesh
  ply_out <- file.path(output_dir, "warped_mesh.ply")
  Rvcg::vcgPlyWrite(warped_mesh, filename = ply_out, binary = TRUE)
  message("Warped mesh saved to: ", ply_out)

  # Plot and annotate mesh
  rgl::open3d()
  rgl::shade3d(warped_mesh, color = "lightblue", alpha = 1)
  rgl::points3d(target_lm, col = "red", size = 5)

  for (view in views) {
    view_file <- file.path(output_dir, paste0("mesh_view_", view, ".png"))
    if (view == "frontal") {
      rgl::view3d(userMatrix = rgl::rotationMatrix(0, 0, 1, 0))
    } else if (view == "lateral") {
      rgl::view3d(userMatrix = rgl::rotationMatrix(pi/2, 0, 1, 0))
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
