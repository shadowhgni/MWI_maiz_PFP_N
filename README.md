# 🌽 Partial Factor Productivity of Nitrogen (PFP-N) in Malawi

<div align="center">

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R](https://img.shields.io/badge/R-%3E%3D4.3-276DC3?logo=r)](https://www.r-project.org/)
[![Data: LSMS-ISA](https://img.shields.io/badge/data-LSMS--ISA-0072BC)](https://www.worldbank.org/en/programs/lsms)
[![Vocabulary: terminag](https://img.shields.io/badge/vocabulary-terminag-2E8B57)](https://github.com/controvoc/terminag)
[![Trials: Carob](https://img.shields.io/badge/trials-Carob-8B4513)](https://carob-data.org/)

**How does the partial factor productivity of nitrogen vary across space and time in Malawi?**

*Five LSMS-ISA survey rounds (2010–2025) · Carob on-farm and on-station trials · terminag-standardised variables · Word and HTML report regenerated on every run*

</div>

---

## 📋 Overview

This workflow harmonises Malawi's Integrated Household Surveys (IHS3–IHS6) and the 2013 round of the Integrated Household Panel Survey (IHPS) with agronomic trials compiled by [Carob](https://carob-data.org/). It also has an optional link to [LCAS](https://systems-agronomy.github.io/lcas/) surveys. Everything is put into one plot-level dataset with [terminag](https://github.com/controvoc/terminag) variable names, which is then used to study the **partial factor productivity of nitrogen**:

$$\text{PFP-N} = \frac{\text{grain yield (kg/ha)}}{\text{mineral N applied (kg N/ha)}}$$

Before answering the research question, the workflow assesses how robust the data are: by quality flag, by source (survey, on-farm trial, on-station trial) and by place. It also tests whether trials can serve as a benchmark for the quality of the survey data.

| | |
|---|---|
| 🌾 **Crop / unit** | Sole maize · one row per plot-crop × yield source |
| 🗓️ **Survey rounds** | IHS3 2010/11 · IHPS 2013 · IHS4 2016/17 · IHS5 2019/20 · IHS6 2024/25 |
| 🧪 **Trials** | Carob `agronomy` and `survey` collections, filtered to Malawi |
| 📏 **N window** | 10–300 kg N/ha (PFP-N is unstable below, implausible above) |
| 📐 **Response curves** | Hyperbolic (Eq 02) vs exponential, with marginal response and bootstrap CIs |
| 🤖 **Machine learning** | Random forests (`caret` + `ranger`), random / spatial-block / leave-one-wave-out CV, SHAP for every plot |
| 🎯 **Efficiency** | AUE-N from trial controls and from propensity-score matching of survey plots |
| 📄 **Report** | `pfpn_report.docx` + `pfpn_report.html`, numbers, tables and figures regenerated each run |

---

## 🚀 Quick Start

```r
# 1. Install the packages (once)
install.packages(c("tidyverse", "xml2", "readxl", "haven", "sf", "terra", "geodata",
                   "mgcv", "glmmTMB", "quantreg", "caret", "ranger", "iml",
                   "emmeans", "multcomp", "multcompView", "MatchIt",
                   "patchwork", "scales", "rmarkdown", "knitr"))

# 2. Edit two things (see "What you must edit" below):
#    - run_all.R      : pfpn_root <- ".."
#    - pfpn_compile.R : paths to your five LSMS .zip files

# 3. Run everything from the scripts folder
setwd("path/to/PFP_N_MWI/pfpn_scripts")
source("run_all.R")
```

The report appears at `PFP_N_MWI/pfpn_analysis_2/pfpn_report.docx` (and `.html`).

> ⏱️ The **first run** downloads codebooks, trials and covariates and fits the random forests and SHAP values for every plot. Expect **an hour or more**. Later runs reuse everything already downloaded.

---

## 📁 Folder Architecture

The scripts live in `pfpn_scripts/`. All other folders are **created by the scripts** one level up, in the project root, the first time they run.

```
MWI_maiz_PFP_N/                    ← project root  (pfpn_root = ".." seen from pfpn_scripts/)
│
├── pfpn_scripts/                  ← YOU KEEP THESE · working directory when running
│   ├── run_all.R                  # runs the four steps in order
│   ├── 00.fpn_utils.R             # shared helpers (paths, downloads, letters, SHAP, ...)
│   ├── 01.pfpn_variable_lookup.R  # step 1 · variable spec checked against each codebook
│   ├── 02.pfpn_compile.R          # step 2 · LSMS + LCAS + Carob → unified dataset
│   ├── 03.pfpn_analysis.R         # step 3 · data robustness, GAMs, variance components
│   ├── 04.pfpn_analysis_2.R       # step 4 · curves, random forests, SHAP, trial zones, AUE-N
│   ├── pfpn_report.Rmd            # report template (rendered by step 4)
│
├── pfpn_lookup/                   ← step 1 (rebuildable)
│   ├── raw/                       #   DDI codebooks (ddi_<id>.xml), LCAS forms, terminag zip
│   ├── spec_check.csv             #   every mapped variable: ok / label_mismatch / missing
│   ├── data_dictionary_<wave>.csv #   full variable list of the files used, per wave
│   └── value_labels.csv, *_stub.csv
│
├── pfpn_config/                   ← ✋ CURATED BY YOU · back it up
│   ├── source_map.csv             #   which source variable plays which role
│   ├── harvest_unit_kg.csv        #   harvest units → kg
│   ├── fert_unit_kg.csv           #   fertilizer units → kg
│   ├── area_unit_m2.csv           #   area units → m²
│   ├── intercrop_fraction.csv     #   share-of-plot categories → fraction
│   ├── fertilizer_n.csv           #   fertilizer products → % N
│   ├── crop_map.csv               #   crop labels → terminag crop names
│   └── district_centroids.csv     #   district centroids (IHS6 has no EA coordinates)
│
├── pfpn_compiled/                 ← step 2 (rebuildable)
│   ├── raw/                       #   extracted .dta files, Carob zips, covariate cache
│   ├── pfpn_unified.rds           #   the dataset, every column labelled
│   ├── pfpn_unified.csv           #   same data as CSV
│   ├── pfpn_unified_dictionary.csv#   data dictionary (labels, units, derivation, terminag)
│   ├── data_quality_statement.md  #   data-quality statement with computed indicators
│   └── summary_by_source.csv, missing_conversions.csv, crop_stand_summary.csv
│
├── pfpn_analysis/                 ← step 3 (rebuildable)
│   ├── figures/  tables/          #   every figure has its data in tables/fig_<name>.csv
│   └── summary.md
│
├── pfpn_analysis_2/               ← step 4 (rebuildable)
│   ├── figures/  tables/          #   every figure has its data in tables/fig_<name>.csv
│   ├── objects/                   #   results.rds, rf_final_model.rds
│   ├── pfpn_report.docx           #   📄 THE REPORT
│   └── pfpn_report.html
│
└── covariates/                    ← optional: drop your own GeoTIFFs here (file name = predictor name)
│
├── LICENSE
└── README.md
```

> 🗑️ You can delete `pfpn_lookup/`, `pfpn_compiled/`, `pfpn_analysis/` and `pfpn_analysis_2/` at any time: they are rebuilt.
> ✋ **Do not delete `pfpn_config/` or `covariates/`** without a backup. They hold values you typed in by hand and your own rasters. The scripts would recreate empty stubs, but your entries would be lost.

The LSMS zip files can live **anywhere**, even outside the project. Their paths are set in `pfpn_compile.R`.

---

## ✏️ What You Must Edit

Only two places need changing before the first run.

### 1. `run_all.R`: where the output folders go

Near the top:

```r
pfpn_root <- ".."   # project root, seen from pfpn_scripts/
```

With `".."`, all `pfpn_*` folders are created next to `pfpn_scripts/`, as in the architecture above. With `"."` they would be created inside `pfpn_scripts/`.

### 2. `pfpn_compile.R`: paths to the LSMS zip files (lines ~68–72)

```r
lsms_mwi_2010_zip <- '../../../Farm sizes across Africa/data/raw/web_scrapped/survey_data/LSMS_Malawi_2010/MWI_2010_IHS-III_v01_M_STATA8.zip'
lsms_mwi_2013_zip <- '.../MWI_2010-2013_IHPS_v01_M_Stata.zip'
lsms_mwi_2016_zip <- '.../MWI_2016_IHS-IV_v04_M_STATA14.zip'
lsms_mwi_2019_zip <- '.../MWI_2019_IHS-V_v06_M_Stata.zip'
lsms_mwi_2024_zip <- '.../MWI_2024-2025_IHS-VI_v01_M_STATA14.zip'
```

- **Relative paths are resolved from `pfpn_scripts/`** (the working directory), not from the project root.
- Full paths are the safest choice. Use forward slashes, even on Windows: `'C:/Users/me/data/MWI_2019_IHS-V_v06_M_Stata.zip'`.
- **Keep the zips exactly as downloaded.** The scripts read the `.dta` files directly from them, and handle the IHS3 zip, which holds both a `Full_Sample/` and a `Panel/` copy (the full sample is used).
- A wave whose zip is not found is **skipped with a message**, not silently.

### Optional settings

| Where | Setting | Default | Purpose |
|---|---|---|---|
| `pfpn_compile.R` | `carob_countries` | `"Malawi"` | countries kept from Carob |
| `pfpn_compile.R` | `lcas_datasets` | *(empty)* | add LCAS files here (see below) |
| `pfpn_analysis_2.R` | `N_min`, `N_max` | 10, 300 | N window (kg N/ha) |
| `pfpn_analysis_2.R` | `n_trees`, `n_boot`, `shap_nsim` | 500, 200, 20 | lower for a quick test run (e.g. 100, 20, 5) |
| `pfpn_analysis_2.R` | `match_km`, `match_years` | 30, 2 | trial zones around Carob sites |
| `pfpn_analysis_2.R` | `psm_caliper` | 0.2 | propensity-score caliper (SD of the logit score) |
| `pfpn_analysis_2.R` | `download_covariates`, `download_chirps` | `TRUE` | environmental covariates from the internet |
| `run_all.R` | `pfpn.refresh` | `"if_updated"` | re-download only when the server has a newer file (`"never"`, `"always"`) |

---

## 📊 Data: How to Obtain Every Resource

All resources are **free**. The LSMS microdata need a free account and acceptance of the terms of use. Everything else is downloaded automatically.

### 🔑 Manual download (free account, login required)

**Malawi LSMS-ISA microdata: World Bank Microdata Library**

1. Create a free account at <https://microdata.worldbank.org> and sign in.
2. Open each study page below, go to the **Get Microdata** tab, accept the terms of use, and download the **Stata** version.
3. Save the zip files anywhere, and put their paths in `pfpn_compile.R`.

| Round | Study page | Variable documentation | Used for |
|---|---|---|---|
| IHS3 2010/11 | [catalog/1003](https://microdata.worldbank.org/catalog/1003) | [data dictionary](https://microdata.worldbank.org/catalog/1003/data-dictionary) | `Full_Sample/` only |
| IHPS 2013 | [catalog/2248](https://microdata.worldbank.org/catalog/2248) | [data dictionary](https://microdata.worldbank.org/catalog/2248/data-dictionary) | the 2013 round (`*_13` files) only |
| IHS4 2016/17 | [catalog/2936](https://microdata.worldbank.org/catalog/2936) | [data dictionary](https://microdata.worldbank.org/catalog/2936/data-dictionary) | |
| IHS5 2019/20 | [catalog/3818](https://microdata.worldbank.org/catalog/3818) | [data dictionary](https://microdata.worldbank.org/catalog/3818/data-dictionary) | |
| IHS6 2024/25 | [catalog/8507](https://microdata.worldbank.org/catalog/8507) | [data dictionary](https://microdata.worldbank.org/catalog/8507/data-dictionary) | no EA coordinates → district centroids |

> ℹ️ The IHPS 2010 round (`*_10` files) re-interviews IHS3 households already in catalog 1003, so it is ignored to avoid double counting. IHS2 2004/05 (catalog 2307) is not used, because its harvest is not recorded per plot.

### 🌐 Downloaded automatically (no login)

| Resource | Source | Used by | Cached in |
|---|---|---|---|
| DDI codebooks of each round | `microdata.worldbank.org/metadata/export/<id>/ddi` | step 1 | `pfpn_lookup/raw/` |
| terminag vocabulary | [github.com/controvoc/terminag](https://github.com/controvoc/terminag) | steps 1, 2 | `pfpn_lookup/raw/` |
| Carob trial collections (`agronomy`, `survey`) | [carob-data.org](https://carob-data.org/) · `geodata.ucdavis.edu/carob/` | step 2 | `pfpn_compiled/raw/` |
| LCAS module forms (xlsforms) | [systems-agronomy.github.io/lcas](https://systems-agronomy.github.io/lcas/) | step 1 | `pfpn_lookup/raw/` |
| GADM district boundaries | [gadm.org](https://gadm.org/) via `geodata` | steps 2–4 | `pfpn_compiled/raw/` |
| Elevation, WorldClim climate | [worldclim.org](https://www.worldclim.org/) via `geodata` | step 4 | `pfpn_compiled/raw/covariates/` |
| SoilGrids (pH, SOC, N, clay) | [soilgrids.org](https://soilgrids.org/) via `geodata` | step 4 | `pfpn_compiled/raw/covariates/` |
| CHIRPS monthly rainfall | [chc.ucsb.edu/data/chirps](https://www.chc.ucsb.edu/data/chirps) | step 4 | `pfpn_compiled/raw/covariates/` |

If a download fails (no internet, server down), the step **continues without that resource** and says so; the codebooks are the only exception. For a codebook, the script stops and tells you where to download it by hand: from the study page, save the **DDI/XML** metadata export as `pfpn_lookup/raw/ddi_<id>.xml`.

### 📦 Optional

- **LCAS survey data.** Obtain the datasets as described on the [LCAS website](https://systems-agronomy.github.io/lcas/), rename the columns with its `rename_lcas.R`, and list each file in `lcas_datasets` in `pfpn_compile.R`.
- **Your own rasters.** Any GeoTIFF in `covariates/` becomes a random-forest predictor named after the file.

---

## 🔄 Pipeline

| Step | Script | What it does | Main outputs |
|---|---|---|---|
| 1 | `pfpn_variable_lookup.R` | Lists, per survey round, the variables needed (household, plot, crop, harvest, area, fertilizer type, quantity, unit, enumerator kg, organic inputs, location). It checks each against that round's codebook: does it exist, and does its label match? | `pfpn_lookup/spec_check.csv`, `pfpn_config/source_map.csv` (written only if absent) |
| 2 | `pfpn_compile.R` | Reads the `.dta` files from the zips. Converts harvest, fertilizer and area units, computes N from product N contents, and keeps sole crops. Takes EA coordinates (or district centroids for IHS6), removes coordinates outside Malawi, adds Carob trials, and labels every column. | `pfpn_compiled/pfpn_unified.rds` + dictionary + data-quality statement |
| 3 | `pfpn_analysis.R` | Data robustness: quality flags, reported vs GPS area, fertilizer kg checks, a trial yield frontier, robust GAMs over space and time, variance components and sensitivity analyses | `pfpn_analysis/` |
| 4 | `pfpn_analysis_2.R` | Coverage funnel per wave; violin + box plots with ANOVA/emmeans letters; Eq 02 vs exponential response curves; random forests with three CV schemes; SHAP for every plot; trial zones and AUE-N with propensity-score matching. Renders the report. | `pfpn_analysis_2/`, `pfpn_report.docx/.html` |

Run a single step on its own (once the earlier steps have run at least once):

```r
setwd("path/to/PFP_N_MWI/pfpn_scripts")
options(pfpn.root = "..")          # same as pfpn_root in run_all.R
source("pfpn_analysis_2.R")
```

---

## 🧭 Your Part: Curating `pfpn_config/`

After each run of step 2, check **`pfpn_compiled/missing_conversions.csv`**. It lists every unit, product or category that has no conversion value yet. Fill in the empty cells in the corresponding file in `pfpn_config/`, then run again. New labels found in the data are **appended** to these files; your entries are never overwritten.

| File | Fill in | Example |
|---|---|---|
| `harvest_unit_kg.csv` | kg per harvest unit, by crop and condition | `50 KG BAG` (shelled maize) = 50 |
| `fert_unit_kg.csv` | kg per fertilizer unit | `2 KG BAG` = 2 |
| `area_unit_m2.csv` | m² per area unit | `ACRE` = 4046.86 |
| `intercrop_fraction.csv` | share of plot for each category | `1/2` = 0.5 |
| `fertilizer_n.csv` | % N per product | `UREA` = 46 · `23:21:0+4S` = 23 |
| `crop_map.csv` | terminag crop name for each label | `MAIZE LOCAL` → `maize` |

Labels such as "50 KG BAG", "KILOGRAM", "ACRE" and N:P:K grades are pre-filled automatically. Local units (pails, ox-carts, basins) need the official IHS conversion factors.

Also review **`pfpn_lookup/spec_check.csv`**: `label_mismatch` rows deserve a look (the variable exists but its label is unexpected), and `missing` rows become `TODO` in `source_map.csv` and are skipped.

---

## 📈 Key Outputs

| File | Description |
|---|---|
| `pfpn_compiled/pfpn_unified.rds` | Plot-level dataset; every column has `label`, `unit`, `description` and `derivation` attributes |
| `pfpn_compiled/pfpn_unified_dictionary.csv` | Data dictionary: label, type, unit, terminag term, derivation, observed values, % missing |
| `pfpn_compiled/data_quality_statement.md` | Error sources for N, yield and PFP-N, with indicators computed per source |
| `pfpn_analysis/summary.md` | Robustness findings and space-time results |
| `pfpn_analysis_2/pfpn_report.docx` | **Full report** with summary, methods, results, discussion (Burke et al. and others) and references |
| `pfpn_analysis_2/tables/fig_*.csv` | The exact data behind every figure |
| `pfpn_analysis_2/objects/rf_final_model.rds` | Final random forest (`caret` object) |

---

## 🔬 Methods at a Glance

- **PFP-N and plot area.** For sole crops, yield and N are both expressed per planted hectare, so PFP-N equals harvest kg ÷ N kg. Errors in plot area therefore cancel out of PFP-N, although they still affect yield and N rate separately.
- **Fertilizer kg.** The enumerator's total kg is used first, then quantity × unit weight. Differences above 10 % are flagged.
- **Response curves.** PFP-N = A_min + A_max / (N + N_0) (Eq 02) and PFP-N = A_min + A_max·e^(−N/N_0) are fitted on the log scale. They are compared by AIC and by cross-validation that holds out whole EAs. Each fit also gives the implied marginal yield response (kg grain per extra kg N).
- **Random forests.** Three nested predictor sets: survey; + environment; + spatial signature (coordinates and oblique coordinates). Each is evaluated with random, spatial-block and leave-one-wave-out CV. SHAP values are computed for every plot and cross-checked with `iml::Shapley`.
- **Group comparisons.** One-way ANOVA → `emmeans` → `multcomp::cld` (Bonferroni; "a" = highest mean). Pairwise Wilcoxon letters are kept as tables.
- **Trials as a benchmark.** Trials are researcher-managed and are not ground truth for farmers' fields. They give an attainable-yield frontier, and a like-for-like comparison within 30 km and ±2 years of each trial.
- **AUE-N** = (yield with N − yield without N) / N. In trials, each plot is compared with the N = 0 plots of the same trial; in surveys, fertilized plots are matched to unfertilized plots by propensity score (`MatchIt`), within trial zones and nationwide. Matching balances only observed characteristics.

---

## 🛠️ Troubleshooting

| Message | Cause and fix |
|---|---|
| `LSMS zip not found, wave skipped` | Wrong path in `pfpn_compile.R`. Relative paths start from `pfpn_scripts/`. |
| `No LSMS or LCAS source could be read` | None of the zip paths is right. |
| `No usable DDI codebook for catalog <id>` | Codebook download failed. Save the DDI/XML export by hand as `pfpn_lookup/raw/ddi_<id>.xml`. |
| A wave is missing from the figures | Open `pfpn_analysis_2/figures/1_wave_coverage_funnel.png`: it shows the step where the wave loses its rows. |
| Few or no yields for a wave | Harvest units without kg weights. Fill in `pfpn_config/harvest_unit_kg.csv`. |
| `Districts without a centroid match` | Add the district name to `district_aliases` in `pfpn_utils.R`. |
| `rm(list = ls())` at the top of scripts | Supported: each step runs in its own environment. |

---

## 📚 Key References (for discussuion only)

- Burke, W.J., Snapp, S.S., Jayne, T.S. (2020). An in-depth examination of maize yield response to fertilizer in Central Malawi reveals low profits and too many weeds. *Agricultural Economics* 51(6): 923–940.
- Burke, W.J., Jayne, T.S., Snapp, S.S. (2022). Nitrogen efficiency by soil quality and management regimes on Malawi farms: Can fertilizer use remain profitable? *World Development* 152: 105792.
- Burke, W.J., Snapp, S.S., Peter, B.G., Jayne, T.S. (2022). Sustainable intensification in jeopardy: Transdisciplinary evidence from Malawi. *Science of the Total Environment* 837: 155758.
- Ragasa, C., et al. (2025). Maize yield responsiveness and profitability of fertilizer: New survey evidence from six African countries. *Food Policy* 133: 102815.
- Jayne, T.S., Mason, N.M., Burke, W.J., Ariga, J. (2018). Taking stock of Africa's second-generation agricultural input subsidy programs. *Food Policy* 75: 1–14.

The full reference list is in the report.

---

## 📝 Citation

No DOI has been assigned to this workflow yet; one will be added here if the workflow is deposited. Until then, please cite the repository and the data sources listed below.

The workflow (scripts, report template and README) was developed by D. Hougni (CIMMYT) with assistance from Claude AI (Anthropic, 2026).

---

## 📄 License

The **code** in this repository (R scripts, report template and documentation) is released under the **GNU General Public License v3.0** ([GPL-3.0](https://www.gnu.org/licenses/gpl-3.0)). You may use, modify and redistribute it; derived works must be distributed under the same licence.

### Reused data

**No data are redistributed with this workflow.** The scripts download or read every dataset from its original provider, and each dataset remains subject to its **own licence and terms of use**, which users must accept and respect, including any citation requirements:

| Data | Provider | Terms |
|---|---|---|
| Malawi LSMS-ISA microdata (IHS3–IHS6, IHPS) | World Bank Microdata Library | terms of use accepted at download; no redistribution of the microdata |
| DDI codebooks and data dictionaries | World Bank Microdata Library | see the Microdata Library's terms |
| LCAS module forms and survey data | [LCAS](https://systems-agronomy.github.io/lcas/) | see the LCAS website and each dataset's licence |
| Carob trial collections | [Carob](https://carob-data.org/) | each Carob dataset keeps the licence of its original publication; cite the original datasets |
| terminag vocabulary | [controvoc/terminag](https://github.com/controvoc/terminag) | see the repository's licence |
| GADM boundaries | [GADM](https://gadm.org/) | see the GADM licence |
| Elevation and climate | [WorldClim](https://www.worldclim.org/) | see WorldClim's terms |
| Soil properties | [SoilGrids](https://soilgrids.org/) | see SoilGrids' licence |
| Rainfall | [CHIRPS](https://www.chc.ucsb.edu/data/chirps) | see CHIRPS' terms |

The GPL-3.0 licence of the code does not extend to any of these datasets, nor to outputs derived from them, which inherit the conditions of the data they come from.

---

## 📞 Contact

* 🏢  **D. Hougni (CIMMYT):** d.hougni@cgiar.org
* 📧  **D. Hougni (Personal):** shadowhgni@yahoo.fr

---

<div align="center"><i>Last updated: October 2026</i></div>
