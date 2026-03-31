# Engineering-Thesis
# TIL Texture Analysis via GLCM for Cancer Patient Survival Prediction

MATLAB codebase for analyzing Tumor Infiltrating Lymphocyte (TIL) images using Gray-Level Co-occurrence Matrix (GLCM) texture analysis, combined with clinical survival data from TCGA.

---

## Project Overview

The goal of this project is to evaluate the prognostic value of a structural texture metric (`struct_ness`) derived from TIL images. For each combination of preprocessing and GLCM parameters, the following statistics are computed:

- **Harrell's C-index** — concordance measure between predicted scores and actual survival times,
- **Hazard Ratio** (HR) and its logarithm,
- **P-value** from the Cox proportional hazards model,
- **Confidence intervals** for HR (in both linear and log scale).

Results for each analyzed cancer type are saved to an `.xlsx` file.

---

## Repository Structure

```
DaneTestowe/
├── Kod/
│   ├── Main.m              # Entry point — path configuration and analysis launcher
│   ├── GLCM.m              # Iterates over all parameter combinations, aggregates results
│   ├── Process.m           # Image preprocessing and GLCM computation
│   ├── Analiza.m           # Survival analysis (Cox PH, C-index, HR, p-value)
│   └── concordanceIndex.m  # Harrell's C-index computation
└── Test data/
    ├── TCGA-CDR-SupplementalTableS1.xlsx   # Patient clinical data (TCGA)
    └── Zdjęcia testowe/                    # TIL images in PNG format (TCGA)
```

---

## Requirements

- **MATLAB** (tested on R2021a or newer)
- **Image Processing Toolbox** — `graycomatrix`, `imgaussfilt`, `imopen`, `imclose`, `medfilt2`, `colfilt`
- **Statistics and Machine Learning Toolbox** — `coxphfit`

---

## Getting Started

1. Open `Main.m` in MATLAB.
2. Set the path to the folder containing TIL images:
   ```matlab
   path_to_TIL_folder = 'C:\path\to\images\folder\';
   ```
3. Set the path to the Excel file with clinical data:
   ```matlab
   path_to_excel = 'TCGA-CDR-SupplementalTableS1.xlsx';
   ```
4. Set the output folder for results:
   ```matlab
   save_results_to = 'C:\path\to\results\folder\';
   ```
5. In the `typelist` variable, specify the subfolder names corresponding to cancer types:
   ```matlab
   typelist = {'CancerTypeName'};
   ```
6. Run the script — results will be saved as `Harrells_c_index_data_<type>.xlsx`.

---

## Analysis Parameters

The code automatically tests all combinations of the following parameters:

| Parameter         | Values                                                    | Description                                        |
|-------------------|-----------------------------------------------------------|----------------------------------------------------|
| `NLvl`            | 2, 4, 6, 8                                                | Number of gray-level quantization levels in GLCM   |
| `NBack`           | `'n'`, `'y'`                                              | Background removal (white pixels = 255)            |
| `blurFilterType`  | `'n'`, `'average'`, `'gaussian'`, `'median'`, `'geomean'` | Blur filter type applied before GLCM               |
| `imgcut`          | `'n'`, `'y'`                                              | Split image into 4 quadrants before computing GLCM |
| `open_close`      | `'n'`, `'close'`, `'open'`                                | Morphological opening / closing operation          |
| `a`               | 3, 5, 7                                                   | Window size for the average filter                 |
| `s`               | 1, 2, 3                                                   | Standard deviation for the Gaussian filter         |

---

## Input Data Format

**TIL images** (`*.png`):
- Grayscale or RGB images — only the red channel (`img(:,:,1)`) is used.
- Filenames must start with the 12-character TCGA patient barcode (e.g. `TCGA-V4-A9E5`).
- White pixels (value 255 in all channels) are treated as background.

**Clinical data** (`TCGA-CDR-SupplementalTableS1.xlsx`):
- Downloaded from the TCGA Clinical Data Resource (CDR).
- Required columns: `bcr_patient_barcode`, `age_at_initial_pathologic_diagnosis`, `OS`, `OS_time`, `DSS`, `DSS_time`, `DFI`, `DFI_time`, `PFI`, `PFI_time`, `vital_status`, `race`, `ajcc_pathologic_tumor_stage`.

---

## Output Format

For each cancer type, the script produces `Harrells_c_index_data_<type>.xlsx` with the following columns:

| Column                              | Description                                          |
|-------------------------------------|------------------------------------------------------|
| `Description`                       | Parameter combination label                          |
| `Harrells_c_index`                  | Harrell's C-index                                    |
| `Negatywny_Harrells_c_index`        | C-index for the negated predictor sign               |
| `P_values`                          | P-value from the Cox model                           |
| `Hazard_ratio`                      | Hazard ratio (HR)                                    |
| `Hazard_ratio_in_logarithm`         | Log hazard ratio                                     |
| `Confidence_intervals`              | 95% confidence intervals for HR                      |
| `Confidence_intervals_in_logarithm` | 95% confidence intervals for log(HR)                 |

---

## Data Source

Clinical data and histopathological images are sourced from **The Cancer Genome Atlas (TCGA)**:
- Clinical data: [TCGA Clinical Data Resource (CDR)](https://www.cell.com/cell/fulltext/S0092-8674(18)31033-3)
- TIL images: generated using TIL detection tools applied to TCGA whole-slide histopathology images
