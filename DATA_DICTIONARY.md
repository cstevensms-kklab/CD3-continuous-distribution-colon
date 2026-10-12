# Data dictionary

## Cell-level data: `colon_cd3_ycoords_github.csv`

One row represents one segmented nucleus. Both CD3+ and CD3− nuclei are included because all nuclei establish each animal's coordinate range before CD3+ cells are selected.

| Variable | Units or accepted values | Description |
| --- | --- | --- |
| `ID` | Text identifier | Animal identifier used to match the demographics table |
| `cell_label` | `positive`, `negative` | CD3-positive or CD3-negative nucleus |
| `lp_epi` | `epi`, `lp` | Epithelium or lamina propria |
| `y_unscaled` | Original image-coordinate units | Unscaled coordinate along the crypt axis; the script converts this to a dimensionless 0–1 axis |

Cell labels and compartment labels are case-insensitive, and surrounding spaces are ignored. The script also accepts `epithelium` or `epithelial` for `epi`, and `lamina propria` or `lamina_propria` for `lp`.

With `REVERSE_AXIS <- FALSE`, the minimum coordinate within each animal maps to 0 (base) and the maximum maps to 1 (luminal surface). `REVERSE_AXIS <- TRUE` reverses this mapping.

## Subject-level data: `demographics_colon_github.csv`

One row represents one animal. The current spatial-analysis script uses only `ID` and `age_cat`. Remaining variables are included to accompany the manuscript table and allow additional analyses.

| Variable | Units or coding | Description |
| --- | --- | --- |
| `ID` | Text identifier | Subject identifier matching the cell-level table |
| `age_cat` | Young = `0`; Old = `1` | Age category in the example dataset. The accompanying manuscript defined Young as <18 years and Old as >/=18 years. |
| `age` | years | Continuous age |
| `girth` | cm | Waist circumference |
| `hdlc` | mg/dL | High-density lipoprotein cholesterol |
| `tg` | mg/dL | Triglycerides |
| `glu` | mg/dL | Fasting glucose measured in peripheral blood |
| `MCP1` | pg/mL | Monocyte chemoattractant protein 1 |
| `CRP` | mg/L | C-reactive protein |
| `Pl_Defensin` | mg/mL | Alpha-defensin measured in plasma |
| `LBP1` | μg/mL | LPS-binding protein |

For `age_cat`, the script additionally accepts `young` or `y` for Young, and `old` or `o` for Old. Capitalization and surrounding spaces are ignored.

Keep missing optional measurements blank; a missing measurement is not a zero. Missing required IDs, age categories, cell labels, compartment labels, or coordinates must be resolved before analysis.
