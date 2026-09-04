# Comparative age–period–cohort analysis reveals divergent temporal trajectories of breast and colorectal cancer among women in Singapore

This code contains the means by which all statistical analyses and data visualisation
was performed for the above manuscript, which features an age- period-cohort (APC) 
analysis of breast and colorectal cancer incidence among women aged 25–84 in 
Singapore (1968–2023), stratified by ethnicity (Chinese, Malay, Indian). Code includes 
data cleaning, age-standardised rate (ASR) and age-specific rate calculation, figure 
generation, APC and statistical analyses.

## Repository structure

```
.
├── README.md
├── .gitignore
├── scripts/
│   └── apc_final.R   # main analysis script
├── data/                           # see below
└── results/                        # output PDFs will be generated here
```

## Data

This script expects two input files in `data/`:
- `cancer.csv` — cancer registry incidence counts by sex, age group, race group, and period
- `population.csv` — corresponding population denominators by year, sex, age group, and race group

Raw data is **not included in this repository** due to data-sharing restrictions. Place your own
files in `data/` before running the script.

Here's a paragraph you can drop into the README:

Both `cancer.csv` and `population.csv` should be plain comma-separated files with a header row. 
`cancer.csv` must contain the following columns:
`year`: 5-year diagnosis period as a hyphenated range (e.g. `1968-1972`)
`sex`: (`"F"` for females)
`agedxgrp`: 5-year age bands in the format `"25-29 Years"`
`racegrp`: ethnicity represented by the following codes — `CN` (Chinese), `MY` (Malay), `IN` (Indian), and `XX` for any other/unknown group
`site`: cancer site, noted exactly as `"Breast"` or `"Colon & rectum"`
`n (incidence count)`: incidence count column named (which R reads in as `n..incidence.count.`) — this column can contain either whole numbers or the string `"<5"` for suppressed low counts

`population.csv` must contain the following columns:
`year`: single calendar year (not a range)
`sex`: same `"F"` convention
`agegrp`: same `"25-29 Years"` format as `agedxgrp` above
`racegrp`: same race codes
`population`: the raw population count for that year/age/sex/race combination. 

Age bands and race codes in both files need to match exactly between the two, since the script joins incidence and population data on `period`, `agedxgrp`, `racegrp`, and `sex`.

## Setup

Open the project as an RStudio Project with the repo root as the working
directory (or run `setwd(here::here("scripts"))` using the `here` package).
The script will install any missing required packages on first run:

`dplyr`, `splines`, `tidyr`, `ggplot2`, `Epi`, `stringr`, `colorspace`,
`ggpubr`, `showtext`, `epitools`, `purrr`, `viridis`

It also registers the Calibri font via `showtext` for plot text — update the
font paths in the script if Calibri isn't installed at the same location on
your machine.

## Running

From the `scripts/` folder, run `apc_final.R` top to bottom. It
will:
1. Clean and merge incidence and population data
2. Compute descriptive statistics and age-standardised incidence rates (with a
   sensitivity analysis across all possible imputation values for `<5` counts)
3. Fit APC models (via `Epi::apc.fit`) per race group and cancer site
4. Produce Figures 1–5 and Table 1, saved as PDFs in `results/`

## Output

- `ASR_smoothed.pdf` — Figure 1: age-standardised incidence trends over time
- `smooth_age_specific_rate_with_age.pdf` — Figure 2: age-specific rates
- `smooth_age_specific_rate_overtheyears.pdf` — Figure 3: age-specific rates by period
- `component_effects.pdf` — Figure 4: APC component effects
- `pre_vs_post_cohort_effects.pdf` — Figure 5: cohort effects by menopausal age group
- `diff_suppression_counts_smoothed_asr.pdf` — sensitivity analysis
