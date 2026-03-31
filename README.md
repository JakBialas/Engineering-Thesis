# Engineering-Thesis
# TIL Histopathological Image Texture Analysis Using GLCM and Harrell's Concordance Index

Source code for an engineering thesis — MATLAB.

This project performs texture analysis of Tumor Infiltrating Lymphocyte (TIL) histopathological images sourced from the TCGA (*The Cancer Genome Atlas*) database. Texture features are extracted using the Gray-Level Co-occurrence Matrix (GLCM) method and combined with patient clinical data in a Cox proportional hazards model. Survival prediction quality is evaluated using Harrell's concordance index (C-index).

---

## Requirements

- **MATLAB** (R2020b or newer recommended)
- Toolboxes:
  - Image Processing Toolbox (`graycomatrix`, `imgaussfilt`, `medfilt2`, etc.)
  - Statistics and Machine Learning Toolbox (`coxphfit`)

---

## Project Structure

```
DaneTestowe/
├── Kod/
│   ├── Main.m              # Main script that runs the analysis
│   ├── GLCM.m              # Iterates over parameter combinations, calls Process
│   ├── Process.m           # Image preprocessing and GLCM feature extraction
│   ├── Analiza.m           # Cox model, hazard ratio, C-index computation
│   └── concordanceIndex.m  # Harrell's concordance index implementation
└── Test data/
    ├── TCGA-CDR-SupplementalTableS1.xlsx   # Patient clinical data (TCGA)
    └── Zdjęcia testowe/                    # TIL images in .png format (TCGA)
```

---

## Usage

1. Open `Main.m` in MATLAB.
2. Set the paths at the top of the script:

```matlab
% Path to the folder containing TIL images (must end with '\')
path_to_TIL_folder = 'C:\...\DaneTestowe\Test data\';

% Path to the Excel file with clinical data
path_to_excel = 'TCGA-CDR-SupplementalTableS1.xlsx';

% Output folder for results (must end with '\')
save_results_to = 'C:\...\DaneTestowe\';
```

3. Run the script. Results will be saved to an Excel file (`Harrells_c_index_data_<type>.xlsx`) in the specified output folder.

> Note: The analysis covers a large number of parameter combinations and may take tens of minutes to complete.

---

## Method Description

### Image Preprocessing (`Process.m`)

The following steps are applied optionally to each TIL image:

- **Morphological filtering** — opening (`open`) or closing (`close`) operations on the background mask
- **Smoothing** — Gaussian, average, median, or geometric mean filter
- **Background removal** — white pixels (255, 255, 255) are treated as background and excluded
- **Image splitting** — optional division into 4 quadrants before GLCM computation

The GLCM is computed across 4 directions (0°, 45°, 90°, 135°), from which texture measures (M1/M2) are derived.

### Survival Analysis (`Analiza.m`)

Extracted texture features are combined with patient age and fed into a Cox proportional hazards model (`coxphfit`). The following statistics are computed:

- **Harrell's C-index** (and its negation for inverse correlation)
- **Hazard ratio** (HR and log-HR)
- **95% confidence intervals** (in linear and logarithmic scale)
- **p-values**

### Parameter Grid (`GLCM.m`)

The analysis iterates over the following parameter combinations:

| Parameter | Values |
|---|---|
| Gray levels (`NLvl`) | 2, 4, 6, 8 |
| Background removal (`NBack`) | n, y |
| Smoothing filter (`blurFilterType`) | none, average, gaussian, median, geomean |
| Image splitting (`imgcut`) | n, y |
| Morphological filter (`open_close`) | none, close, open |
| Gaussian sigma (`s`) | 1, 2, 3 |
| Average filter window size (`a`) | 3, 5, 7 |

---

## Data

TIL images and clinical data are sourced from the [TCGA](https://www.cancer.gov/ccg/research/genome-sequencing/tcga) project. The file `TCGA-CDR-SupplementalTableS1.xlsx` contains patient survival endpoints (OS, DSS, DFI, PFI).

Image filenames follow the TCGA naming convention — the first 12 characters correspond to the patient barcode, enabling automatic matching of images to clinical records.

---

## Author

Jakub Białas
