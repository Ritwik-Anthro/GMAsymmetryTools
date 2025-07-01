if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "Specimen", "Landmark", "Side", "X", "Y", "Z",
    "PC1", "PC2", "Config",
    "X1", "Y1", "X2", "Y2",
    "X_Left", "Y_Left", "X_Right", "Y_Right",
    "Distance", "Midline",
    "Mean_Procrustes_Distance",    # Added for asymmetry analysis output
    "Centroid_Size",               # Added for centroid size export
    "Mahalanobis_Distance",        # Added for Mahalanobis distance output
    "Region",                      # Added for region-specific data
    "tps_mag",                     # Added for TPS deformation magnitude
    "PCA",                         # Added for PCA-related data
    "Landmark_Displacement",       # Added for displacement analysis
    "Residual_Vectors",            # Added for residual vector analysis
    "Procrustes_Residual_Vectors", # Added for Procrustes residual analysis
    "Asymmetry_Vectors"            # Added for visualizing asymmetry vectors
  ))
}
