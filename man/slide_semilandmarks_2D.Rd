% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slide_semilandmarks_2D.R
\name{slide_semilandmarks_2D}
\alias{slide_semilandmarks_2D}
\title{Slide Semilandmarks for 2D Data}
\usage{
slide_semilandmarks_2D(csv_file)
}
\arguments{
\item{csv_file}{Path to the CSV file containing 2D landmark data. Required columns: Specimen, Landmark, Midline, Region, Side, Slider, CurveID, X, Y.}
}
\value{
A data frame with slid 2D coordinates and associated metadata.
}
\description{
Performs Generalized Procrustes Analysis (GPA) with sliding semilandmarks on 2D landmark data.
}
