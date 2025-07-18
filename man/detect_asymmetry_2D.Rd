% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detect_asymmetry_2D.R
\name{detect_asymmetry_2D}
\alias{detect_asymmetry_2D}
\title{Detect Bilateral Asymmetry in 2D Landmark Data (Region-wise)}
\usage{
detect_asymmetry_2D(file_path, output_base = "asymmetry_output", tps_mag = 2)
}
\arguments{
\item{file_path}{Character. Path to the input CSV file containing 2D landmark data.
The file must include the columns: \code{Specimen}, \code{Landmark}, \code{Midline},
\code{Region}, \code{Side}, \code{X}, and \code{Y}.}

\item{output_base}{Character. Optional base name for the output directory.
A timestamp will be appended to create a unique folder. Default is \code{"asymmetry_output"}.}

\item{tps_mag}{Numeric. Magnification factor for the thin-plate spline (TPS) deformation grid
in the TPS plot. Default is \code{2}.}
}
\value{
(Invisibly) the full path to the output directory. Also writes multiple output files per region:
\itemize{
  \item \code{centroid_sizes.csv} - Centroid size per configuration.
  \item \code{mean_procrustes_distances.csv} - Mean shape difference between left and right sides.
  \item \code{mahalanobis_distances.csv} - Mahalanobis distances per specimen.
  \item \code{procrustes_anova_summary.txt} - Procrustes ANOVA result.
  \item \code{pca_scores.csv} - Principal component scores for each configuration.
  \item \code{pca_variance_explained.csv} - Percentage and cumulative variance explained by PCs.
  \item \code{pca_plot.png} - 2D PCA scatterplot.
  \item \code{asymmetry_vectors.png} - Landmark displacement vectors (Left vs Right).
  \item \code{tps_plot.png} - Thin-plate spline deformation plot.
  \item \code{asymmetry_heatmap.png} - Heatmap showing average landmark-wise asymmetry.
}
}
\description{
This function analyzes 2D landmark coordinate data for bilateral asymmetry on a per-region basis.
It performs Generalized Procrustes Analysis (GPA), calculates Procrustes distances,
performs Procrustes ANOVA, and generates PCA plots, Mahalanobis distances, TPS deformation plots,
asymmetry vectors, and heatmaps. All outputs are saved to a timestamped folder.
}
