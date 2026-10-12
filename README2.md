# Continuous spatial analysis of colonic CD3+ cells

This repository contains an R script and example data for comparing the relative distribution of CD3+ cells along the colonic crypt axis between young and old animals. The script analyzes total mucosa, epithelium, and lamina propria separately.

The normalized axis runs from **0 (crypt base)** to **1 (luminal surface)**. Each animal contributes one probability-density curve to each analysis in which it meets the cell-count requirement. The animal, rather than the individual cell, is the independent unit.

**This analysis tests where CD3+ cells are located, conditional on CD3 positivity. It does not test absolute CD3+ abundance, cells per tissue area, or the proportion of all cells that are CD3+.**

## Repository files

| File | Purpose |
| --- | --- |
| `colon_cd3_spatial_functional_analysis_github.R` | Data preparation, spatial analysis, plots, and results |
| `colon_cd3_ycoords_github.csv` | Cell-level example data, including CD3+ and CD3− nuclei |
| `demographics_colon_github.csv` | Animal IDs, age groups, and additional demographic/biomarker data |
| [DATA_DICTIONARY.md](DATA_DICTIONARY.md) | Column definitions, accepted labels, and demographic units |

The example dataset contains 19 animals: 9 young and 10 old. Additional demographic and biomarker variables are provided to accompany the manuscript table and support further exploration. **The current script uses only `ID` and `age_cat` from the demographics file; it does not analyze or adjust for the other variables.**

## Run the example analysis

1. Download the R script and both CSV files into the same folder. Download the actual files rather than copying the GitHub page text.
2. Open R or RStudio and install the required packages once:

   ```r
   install.packages(c("tidyverse", "gridExtra"))
   ```

3. Set the working directory to the folder containing the files. In RStudio, use **Session → Set Working Directory → Choose Directory**, or run:

   ```r
   setwd("/path/to/repository-folder")
   ```

4. Run the script:

   ```r
   source("colon_cd3_spatial_functional_analysis_github.R")
   ```

Results are written to a `spatial_analysis_outputs` folder in the working directory. Re-running the script in the same output folder replaces files with the same names; choose a different `OUTPUT_DIR` to retain a separate run.

The default run includes 9,999 permutations for each primary analysis plus sensitivity analyses. Runtime depends on the data and computer. For a preliminary check, set `B_MAIN <- 999` and `RUN_SENSITIVITY <- FALSE` in **USER SETTINGS**, then restore the intended final settings before producing results.

## Use your own data

Prepare two comma-delimited CSV files with a header row and no introductory text or key rows mixed into the data. Save missing optional demographic measurements as blank cells rather than treating them as zero.

### 1. Cell-level file

Include one row per segmented nucleus, with these required columns:

| Column | Required contents |
| --- | --- |
| `ID` | Subject identifier matching the demographics file |
| `cell_label` | `positive` for CD3+ or `negative` for CD3− |
| `lp_epi` | `epi` for epithelium or `lp` for lamina propria |
| `y_unscaled` | Finite numeric coordinate along the crypt axis, before 0–1 scaling |

**Include all segmented nuclei in the analyzed tissue, including CD3− nuclei.** The script uses all nuclei to establish each animal's coordinate range before selecting CD3+ nuclei. Providing only CD3+ cells changes the normalization reference.

Coordinates must describe a consistently oriented crypt axis. By default, larger raw coordinates must be closer to the luminal surface. If every specimen has the opposite orientation, set `REVERSE_AXIS <- TRUE`. Correct mixed orientations before running the script.

The script pools all rows with the same `ID` and scales them together. It does not align separate images or crypts, apply rotation, or normalize each ROI separately. If you combine multiple images or ROIs for one subject, first ensure their coordinates share a meaningful common crypt axis. Separate specimens or repeated measurements are not modeled as separate levels of sampling.

### 2. Subject-level demographics file

Include one row per subject with these required columns:

| Column | Required contents |
| --- | --- |
| `ID` | Subject identifier matching the cell-level file |
| `age_cat` | Young/Old group label using one of the accepted encodings below |

| Group | Accepted labels |
| --- | --- |
| Young | `0`, `young`, `y` |
| Old | `1`, `old`, `o` |

Capitalization and surrounding spaces are ignored. The script reads group labels **only from the designated age-group column**. It does not search IDs or other columns for zeros, ones, or age labels.

Group assignments must be supplied by the user. The script does not calculate groups from continuous age or impose an age cutoff. Extra demographic columns are allowed and ignored by this analysis. Subject IDs are read as text; keep their formatting consistent across files, especially if IDs contain leading zeros.

### 3. Set filenames and column mappings

Either use the example filenames and headers or edit the settings near the top of the script. For example:

```r
CELL_FILE <- file.path(ANALYSIS_DIR, "my_cells.csv")
COVARIATE_FILE <- file.path(ANALYSIS_DIR, "my_subjects.csv")

CELL_ID_COLUMN <- "subject_id"
CELL_LABEL_COLUMN <- "CD3_status"
COMPARTMENT_COLUMN <- "compartment"
RAW_Y_COLUMN <- "raw_y"
COVARIATE_ID_COLUMN <- "subject_id"
AGE_GROUP_COLUMN <- "age_group"
```

Changing a column mapping changes which column is read, not the accepted values within it. The current implementation compares two groups labeled Young and Old. Other comparisons require deliberate changes to the group recoding and output labels; more than two groups require changes to the statistical test.

Excel input is also supported. Install `readxl`, change the filenames to `.xlsx` or `.xls`, and set `CELL_SHEET` and `COVARIATE_SHEET` to the appropriate sheet names. For the supplied demographics workbook, use the **Data** sheet. CSV files contain a single table and do not use the sheet settings.

## What the script does

1. Checks required columns, IDs, labels, and coordinates; matches subjects between files.
2. Scales coordinates within each animal using all nuclei: `y_scaled = (y_unscaled - minimum) / (maximum - minimum)`. Scaling precedes CD3+ selection and compartment-specific analyses.
3. Selects CD3+ cells and checks the minimum cell count separately for each analysis.
4. Estimates a boundary-corrected kernel density for each included animal using reflection at 0 and 1. Each curve is normalized to integrate to 1.
5. Selects one age-blind bandwidth per analysis: the median of animal-specific Sheather–Jones bandwidths, with `nrd0` as a fallback when needed.
6. Applies a centered-log-ratio (CLR) transformation and computes a pointwise Welch-type statistic (`t²`). Plotted density curves remain on the original density scale.
7. Permutes age labels at the animal level to obtain a global maximum-statistic p-value and pointwise and simultaneous thresholds.
8. Runs the enabled sensitivity analyses and saves results, quality-control summaries, and R session information.

## Settings and interpretation

| Setting | Default | Meaning |
| --- | --- | --- |
| `GRID_N` | `301` | Number of positions evaluated between 0 and 1 |
| `MIN_CD3_CELLS` | `10` | Minimum CD3+ count per animal within each analysis |
| `B_MAIN` | `9999` | Primary permutation count |
| `ALPHA` | `0.05` | Significance level |
| `RANDOM_SEED` | `20260710` | Seed used to make permutation runs reproducible |
| `USE_CLR_PRIMARY` | `TRUE` | Use CLR-transformed density curves for primary tests |
| `RUN_SENSITIVITY` | `TRUE` | Compare bandwidth multipliers of 0.75 and 1.25 and an analysis without CLR |
| `B_SENSITIVITY` | `1999` | Sensitivity-analysis permutation count |
| `RUN_FDANOVA_SENSITIVITY` | `FALSE` | Optional global `fdANOVA` Fmaxb sensitivity analysis |

The cell-count cutoff should be selected before examining group differences. At least two included animals per group are required for each analysis. A successful run with a small sample is not evidence of adequate statistical power.

In the statistic panel, the **dotted curve** is an unadjusted pointwise threshold; the **dashed horizontal line** is the simultaneous maximum-statistic threshold across the crypt axis. Pointwise crossings are exploratory. The simultaneous threshold addresses multiple positions within one analysis, not multiplicity across all three tissue analyses. The global summary provides Holm-adjusted p-values across all three analyses and, separately, across the two compartment analyses.

Permutation inference requires an appropriate subject-level sampling design and justified exchangeability of group labels under the null. The current script does not adjust for covariates, model paired observations, or perform stratified permutations.

The optional `fdANOVA` test requires that package and returns a global result only. It is disabled by default and is not required for the primary analysis.

## Outputs

| Output | Contents |
| --- | --- |
| `qc_by_animal.csv` | Nucleus counts, CD3+ counts by compartment, and raw coordinate ranges |
| `combined_data_with_y_scaled.csv` | Internal analysis columns and scaled coordinates; additional demographic variables are not retained |
| `global_functional_test_summary.csv` | Group sizes, bandwidths, global p-values, Holm adjustments, and interval counts |
| `sensitivity_analysis_summary.csv` | Sensitivity results, when enabled |
| `primary_analysis_objects.rds` | Primary analysis objects for further work in R |
| `R_sessionInfo.txt` | R and loaded package versions for the run |
| `bandwidth_method_check_*.csv` | Per-animal bandwidth diagnostics for total mucosa, epithelium, and lamina propria |

Each tissue-analysis subfolder also contains cell counts and inclusion decisions, subject bandwidths, animal KDE curves, pointwise results, simultaneous significant intervals, and a combined density/statistic figure in PNG and PDF formats. An empty interval table indicates that no grid interval crossed the simultaneous threshold.

## Troubleshooting

- **File not found:** Check the working directory, input filenames, and `ANALYSIS_DIR`.
- **Missing required column:** Check spelling and capitalization of headers or update the column mappings.
- **Unrecognized label:** Check the accepted values in [DATA_DICTIONARY.md](DATA_DICTIONARY.md). Missing age labels stop the script.
- **Animals without demographics:** Every subject in the cell file needs a demographics entry. Demographics-only subjects generate a warning and are not analyzed.
- **Conflicting age labels:** Provide a single age-group assignment per subject.
- **Invalid coordinate or zero coordinate range:** Check missing/nonnumeric values and confirm that each animal has more than one distinct coordinate.
- **Too few included animals:** Inspect cell-count reports and the prespecified inclusion criterion for the affected compartment.

The documentation reflects the supplied script. Input preparation and matching were checked during repository preparation; a complete execution in R has not yet been verified. The script saves the R environment details for each completed run.
