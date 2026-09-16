# Public replication package (restricted micro-data excluded) [![DOI](https://zenodo.org/badge/1275071083.svg)](https://doi.org/10.5281/zenodo.20788358)


This is the **public** replication package for the paper on estimating the size of the
unauthorised (irregular) foreign population in Poland from aggregated administrative data.

It mirrors the full replication package, **except that two restricted sources are omitted**:
the **tax register** and the **raw police micro-data**. Neither can be shared publicly. All
code and data in this folder are free of both (see *What is omitted* below). Everything else
— in particular the main ZUS-based population estimates and all diagnostics that do not use
the tax register — is fully reproducible here.

## What is omitted (and why)

**1. Tax register.** It cannot be released, so this package drops it from the inputs and from
every output that used it:

- **Data:** `data-raw/tax-records-2019-2024.csv` and `data/poland-tax-register.csv` are not
  included, and the `pop_tax` column is removed from `data/poland-full-database.csv` and
  `data/poland-for-model.csv`.
- **Figure 2** — the *Tax register* panel is dropped (PESEL and Social Insurance only).
- **Table 2** — the *Tax* column is dropped.
- **Appendix register diagnostics** — the Tax `alpha`/`beta` plots (by sex and by continent)
  are dropped (`figA-tax-*`).
- **Cross-register comparisons** (AIC/BIC across specifications, popsize by register, the
  gap table) are computed over **ZUS and PESEL only**; the Tax rows are absent.

**2. Raw police micro-data.** The raw police extract is individual, proceeding-level
micro-data (one row per person/incident, with a unique proceeding id, and including
criminal-proceeding records), so it is restricted and not shared:

- **Data:** `data-raw/police-data-2016-2024.csv` is not included. Only the **aggregated**
  public table `data/poland-police.csv` (counts by year/country/age/sex, with no id and no
  individual records) is shipped.
- **Code:** `codes/2-prepare-data.R` loads the pre-built `data/poland-police.csv` aggregate
  directly instead of rebuilding it from the raw extract; the aggregation block (subset to
  criminal/process.reg proceedings, deduplicate by id, count by PESEL flag) is omitted but
  preserved in the full, non-public package. All police-based outputs (the `police_id_no` /
  `police_id_yes` counts used downstream) are unaffected, since they only ever used these
  aggregates.

The model in the paper uses the **ZUS (social insurance)** register as the reference
population `N`, so none of these omissions affect the main estimates. The modelling frame
(`pop_insured > 0`) is identical to the full package (1,382 community-year-sex cells).

## Repository structure

```
codes/        R scripts, run in order by 0-run-all.R   (tax-stripped)
data-raw/     raw register extracts, EXCLUDING the tax records and raw police micro-data
data/         cleaned community-level tables, EXCLUDING tax (police shipped as aggregate)
results/      cached model fits/bootstraps (*.rds) + simulation outputs
figs/         main-text figures   (fig1-fig7)
figs-appen/   appendix figures    (figA-*)
tables/       LaTeX tables (*.tex)
paper/        manuscript sources (the published figures/tables use the full data)
```

Figures are written to **`figs/`** (main text) and **`figs-appen/`** (appendix), matching the
`\includegraphics` paths in the manuscript.

## How to run

From this folder, in R:

```r
source("codes/0-run-all.R")
```

This sources, in order: `1-functions.R`, `2-prepare-data.R` (builds `data/` from `data-raw/`),
`3-main-paper-analysis.R`, `4-simulation-study.R`, `5-supplement.R`.

Notes:

- The modelling uses the `uncounted` package:
  `pak::pkg_install("ncn-foreigners/uncounted")`.
- Expensive fits, bootstraps and leave-one-out runs are cached under `results/*.rds` and
  reused; set the `recompute_*` flags to `TRUE` to regenerate them. The cross-register caches
  are intentionally not shipped and are rebuilt (tax-free) on first run.
- `4-simulation-study.R` runs 2000 replications and is **slow**; its raw results are provided
  in `results/simulation-results/`. To skip the re-run, comment out the
  `source("codes/4-simulation-study.R")` line in `0-run-all.R`; `5-supplement.R` builds all
  simulation tables (appendix) from the provided outputs.

## Manuscript

`paper/` contains the manuscript sources and the submitted PDF:
[**paper_counting_shadows_submission.pdf**](paper/paper_counting_shadows_submission.pdf)
(main text with the supplementary materials). The manuscript was produced from the full
data, so the few tax-based items above appear there but cannot be regenerated from this
public package.

## Notes on classifications

+ Croatia joined the Schengen Area in 2023; Bulgaria and Romania joined for
  air/sea borders in March 2024 and fully (land borders) on 1 January 2025
+ Cyprus and Ireland are EU members but not in Schengen
+ Iceland, Liechtenstein, Norway, and Switzerland are non-EU Schengen members
+ The United Kingdom was never in the Schengen Area; Brexit removed it from the EU
+ Modelling population: EU free-movement nationals are excluded from
  `full_database_processed` and everything built from it, because they cannot be
  unauthorised in the sense of the model and have zero apprehensions in every cell.
  The rule (codes/2-prepare-data.R) drops Croatia, Bulgaria, Romania, Ireland and
  Cyprus throughout, and the United Kingdom for 2019-2020 (free movement through the
  end of the Brexit transition period; a visa-free third country from 2021)

Country codes:

+ XXX — Stateless (Art. 1 of the 1954 Convention)
+ UNK — Unknown
+ KOSOVO and RKS — Kosovo

## Session info

```r
> sessionInfo()
R version 4.6.0 (2026-04-24)
Platform: aarch64-apple-darwin23
Running under: macOS Sequoia 15.6.1
```
