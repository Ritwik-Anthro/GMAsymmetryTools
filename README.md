# GMAsymmetryTools

**GMAsymmetryTools** is an R package for detecting and visualizing bilateral asymmetry in 2D and 3D landmark data using geometric morphometrics. It also supports thin-plate spline (TPS) mesh warping using 3D landmark coordinates.

---

## Functions

- `detect_asymmetry_2D()`  
  Analyze 2D landmark data for bilateral asymmetry with Procrustes ANOVA, PCA, and vector plots.

- `detect_asymmetry_3D()`  
  Analyze 3D landmark data for bilateral asymmetry. Outputs include Procrustes ANOVA summary, mean shape, PCA plot, interactive 3D visualizations of landmark displacements and residuals.

- `warp_mesh_TPS()`  
  Warp 3D mesh models using TPS deformation from a source to a target landmark configuration. Useful for visualizing shape differences between individuals or averaged configurations.

---

## Input Data Format

### For 2D Asymmetry (`detect_asymmetry_2D`)
| Specimen | Landmark | Side  | X     | Y     |
|----------|----------|-------|-------|-------|
| Spec1    | LM1      | Left  | 12.3  | 45.6  |
| Spec1    | LM1      | Right | 14.1  | 43.7  |

### For 3D Asymmetry (`detect_asymmetry_3D`)
| Specimen | Landmark | Side  | X     | Y     | Z     |
|----------|----------|-------|-------|-------|-------|
| Spec1    | LM1      | Left  | 12.3  | 45.6  | 10.1  |
| Spec1    | LM1      | Right | 14.1  | 43.7  | 10.5  |

Ensure the `Side` column has only `"Left"` and `"Right"` entries (case-sensitive).

---

## TPS Mesh Warping Notes

To perform TPS warping with `warp_mesh_TPS()`:

1. First run `detect_asymmetry_3D()` to generate the **mean shape**:
   - Output: `mean_shape.csv` inside the `asymmetry_output_YYYY-MM-DD_HH-MM/` folder.
   - Move this file manually to your working directory.
   - Format it to retain only the **X**, **Y**, and **Z** columns (remove Landmark column).

2. Create a separate target landmark CSV manually:
   - This should have the same number of landmarks and contain only `X`, `Y`, and `Z` columns.

3. Provide a 3D mesh file (`.ply` or `.stl`) that corresponds to the **source** configuration.

4. Run `warp_mesh_TPS(source_lm_file, target_lm_file, mesh_file)`.

5. The warped mesh is saved as `warped_mesh.ply` in the output directory.

 Best viewed using a mesh viewer (e.g., **MeshLab**, **Blender**, or **3D Viewer** in R using `rgl`).

### DO NOT mix coordinate spaces!
The resulting mean_shape is a shape produced by the average of both Left and mirrored Right landmarks across all cubes. You can use either Left or mirrored Right landmarks as the target for warping, but you must be consistent in interpreting the coordinate space. If you use original (non-mirrored) Right landmarks, then TPS warping will result in distorted or biologically incorrect warps, because you're mixing coordinate systems.

### TPS (Thin Plate Spline) warping is a global deformation method. That means:

1. The entire mesh is warped, not just the area near your landmarks.

2. It minimizes bending energy across the whole mesh to match all corresponding landmarks between the source and target.

3. Even regions far from the landmarks will deform slightly, based on how the landmarks move.

4. If you only care about a specific region, crop or segment the mesh before TPS.

This is by design — TPS assumes landmarks represent the global geometry of the shape.

---

## Folder Outputs

Each function creates a timestamped output folder containing:

- `.csv` results (mean shape, distances)
- `.txt` Procrustes ANOVA summary
- `.png` plots (PCA, vectors, heatmaps)
- `.html` interactive 3D widgets (for 3D functions)

---

## Installation

Coming soon via GitHub. For now, you can install manually:

```r
# Development version
devtools::install_github("Ritwik-Anthro/GMAsymmetryTools")
