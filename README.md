# CD3-continuous-distribution-colon
R workflow for continuous spatial analysis of cell distribution along the colonic crypt axis using normalized coordinates, kernel-density estimation, and functional ANOVA.

This repository contains an R workflow and data for analyzing continuous spatial distributions of cells along the colonic crypt axis. The approach preserves each cell’s spatial coordinate rather than assigning cells only to discrete anatomical layers.

The workflow was developed to compare the distribution of CD3+ T cells between young and old nonhuman primates.

Analysis overview

The workflow:

1. Imports cell-level spatial-coordinate data.
2. Aligns and scales crypt-axis positions to a common 0–1 coordinate system.
3. Calculates a kernel-density estimate for each animal across 100 evenly spaced positions.
4. Visualizes mean spatial-density curves by experimental group.
5. Uses functional ANOVA to test whether the complete spatial distributions differ between groups.
6. Identifies locations along the crypt axis where group differences are most pronounced.

Input data

The input file should contain one row per detected cell and include, at minimum:

- "animal_id": unique animal identifier
- "age_group": experimental group, such as "Young" or "Old"
- "y_scaled": normalized cell position along the crypt axis, ranging from 0 to 1

Additional variables, such as sex, tissue compartment, ROI, or marker status, may also be included.

Coordinate orientation:
"0 = crypt base"
"1 = luminal surface"

Main script

"colon_cd3_spatial_functional_analysis_public.R"

The script contains the complete workflow for data preparation, density estimation, statistical testing, and figure generation.

Requirements

The analysis was performed in R using packages including:

- "tidyverse"
- "MASS"
- "fdANOVA"

Install the required packages with:

install.packages(c("tidyverse", "MASS", "fdANOVA"))

Running the analysis

1. Download or clone this repository.
2. Place the input data in the appropriate data folder.
3. Update the input and output paths at the beginning of the R script.
4. Run the script from beginning to end.

Outputs

The workflow generates:

- Per-animal kernel-density estimates
- Group-level spatial-density curves
- Functional ANOVA results
- Location-specific group comparisons
- Sensitivity analysis
- Publication-ready figures

Statistical considerations

Each animal is treated as the independent experimental unit. Kernel-density estimates are calculated separately for each animal before group-level comparisons. Functional ANOVA is then used to test whether spatial-distribution curves differ between experimental groups across the complete normalized axis.

Data availability

"[State whether example data, processed data, or no data are included. If the original data are unpublished, consider providing a small synthetic example dataset.]"

Citation

If you use or adapt this workflow, please cite:
Citation information will be added following publication.

Contact

For questions about this workflow, contact:

Christina Stevens, Wake Forest University School of Medicine, christina.stevens@wfusm.edu
