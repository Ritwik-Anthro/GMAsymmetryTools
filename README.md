
# GMAsymmetryTools

**GMAsymmetryTools** is an R toolkit for geometric morphometric analysis
of bilateral asymmetry in 2D and 3D landmark data. It includes pipelines
for sliding semilandmarks, reformatting, Procrustes ANOVA, PCA,
Mahalanobis distance analysis, TPS mesh warping, and interactive
visualization.

------------------------------------------------------------------------

## Pipeline Overview

1.  **Prepare your landmark dataset** with the following columns:

| Specimen | Landmark | Midline | Region | Side | Slider | CurveID | X | Y | Z\* |
|----|----|----|----|----|----|----|----|----|----|
| *(Z optional: omit for 2D)* |  |  |  |  |  |  |  |  |  |

- The dataset should be a single file containing both left side and
  right side data.
- Specimen contains the unique individual ID (e.g., Spec1)
- Landmark contains the landmark label (e.g., LM1)
- Ensure that landmark label is same for both left and right side for a
  particular landmark.
- Midline can be 1 or 0 (where 1 = midline landmark).
- Region = anatomical region (thumb, palm, etc.). Use numeric values for
  each region (e.g., 1,2,…)
- Ensure the `Side` column has only `"Left"` and `"Right"` entries
  (case-sensitive).
- Slider can be 0 or 1. Slider = 1 indicates the landmark will be slid.
- CurveID = ID for semilandmark curve grouping (e.g., 1,2,…)

2.  **Slide semilandmarks** using:
    - `slide_semilandmarks_2D()` for 2D data
    - `slide_semilandmarks_3D()` for 3D data
3.  **Reformat the dataset** using:
    - `reformat_landmark_csv()` to sort and standardize landmarks.
4.  **Analyze asymmetry** using:
    - `detect_asymmetry_2D()` – for 2D data
    - `detect_asymmetry_3D()` – for 3D data
5.  **Visualize shape difference via TPS warping**:
    - `warp_mesh_TPS()` – warps a 3D mesh based on landmark deformation.
6.  *(Optional)* Visualize landmark displacement:
    - `visualize_landmark_vectors()`

------------------------------------------------------------------------

## Functions

### Sliding Semilandmarks

- `slide_semilandmarks_2D(csv_file)`  
  Applies Generalized Procrustes Analysis with sliding to 2D
  semilandmarks. Preserves full metadata.

- `slide_semilandmarks_3D(csv_file)`  
  Same as above, but for 3D data.

### Reformat

- `reformat_landmark_csv(input_csv, output_csv)`  
  Sorts landmarks by specimen and side using natural order; preserves
  metadata.

### Asymmetry Analysis

- `detect_asymmetry_2D(file_path, output_base = "asymmetry_output", tps_mag = VALUE)`  
  Region-wise 2D asymmetry analysis. Exports:
  - `centroid_sizes.csv`
  - `mean_procrustes_distances.csv`
  - `mahalanobis_distances.csv`
  - `pca_scores.csv`
  - `pca_variance_explained.csv`
  - `asymmetry_heatmap.png`
  - `pca_plot.png`
  - `landmark_displacement_arrows.html`
    - Each landmark is shown in blue (Left) and red (Right).
    - A green arrow connects the Left and Right versions of each
      landmark.
  - `procrustes_residual_vectors.html`
    - Gray spheres represent the mean shape (mshape()).
    - Purple arrows extend from the mean to the corresponding landmark
      in the individual’s configuration.
- `detect_asymmetry_3D(file_path, output_base = "asymmetry_output")`  
  Region-wise 3D asymmetry analysis. Exports:
  - `centroid_sizes.csv`
  - `mean_procrustes_distances.csv`
  - `mahalanobis_distances.csv`
  - `pca_scores.csv`
  - `pca_variance_explained.csv`
  - `mean_shape.csv`
  - `pca_plot.png`
  - `tps_plot.png`
  - `landmark_displacement_arrows.html`
    - Each landmark is shown in blue (Left) and red (Right).
    - A green arrow connects the Left and Right versions of each
      landmark.
  - `procrustes_residual_vectors.html`
    - Gray spheres represent the mean shape (mshape()).
    - Purple arrows extend from the mean to the corresponding landmark
      in the individual’s configuration.

### Mesh Warping

- `warp_mesh_TPS(source_lm_file, target_lm_file, mesh_file, output_base, views, warp_magnitude = VALUE)`  
  TPS-warp a 3D mesh using source and target landmark files. Exports:
  - `warped_mesh.ply`
  - Frontal/lateral snapshots
  - Interactive HTML viewer

WARNING! The **source and target** landmarks must:

- Contain only `X`, `Y`, `Z`

- Be in the same coordinate space.

  - The resulting mean_shape is a shape produced by the average of both
    Left and mirrored Right landmarks across all cubes. - You can use
    either Left or mirrored Right landmarks as the target for warping,
    but you must be consistent in interpreting the coordinate space.
  - If you use original (non-mirrored) Right landmarks, then TPS warping
    will result in distorted or biologically incorrect warps, because
    you’re mixing coordinate systems.

- Match the geometry of the mesh

  - TPS (Thin Plate Spline) warping is a global deformation method. That
    means:

    - The entire mesh is warped, not just the area near your landmarks.
    - It minimizes bending energy across the whole mesh to match all
      corresponding landmarks between the source and target.
    - Even regions far from the landmarks will deform slightly, based on
      how the landmarks move.
    - If you only care about a specific region, crop or segment the mesh
      before TPS.
    - This is by design — TPS assumes landmarks represent the global
      geometry of the shape.

### 3D Landmark Vector Viewer

- `visualize_landmark_vectors(source_lm_file, target_lm_file, output_file)`  
  Visualize 3D displacement vectors between two aligned landmark sets.
  Saves interactive HTML.
  - Black points: source landmarks
  - Red points: target landmarks
  - Blue arrows: direction and magnitude of transformation

------------------------------------------------------------------------

## Folder Outputs

Each asymmetry or TPS function saves a timestamped folder with: - `.csv`
files: shapes, distances, PCA - `.png` plots: PCA, TPS, heatmaps -
`.txt` summaries: Procrustes ANOVA - `.html` interactive widgets (3D
only)

------------------------------------------------------------------------

## Installation

Coming soon on GitHub.

In the meantime:

``` r
# Development installation
devtools::install_github("Ritwik-Anthro/GMAsymmetryTools")
```

------------------------------------------------------------------------

## Citation

A citation template will be made available soon.
