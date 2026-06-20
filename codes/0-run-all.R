## Master script: run the full pipeline from the project root.
## Loads packages, then sources every stage in order. Each stage relies on the
## objects created by the previous one (they share this session's environment).

## one-time install of the modelling package (uncomment on first run):
# install.packages("pak")
# pak::pkg_install("ncn-foreigners/uncounted")

dir.create("data",    showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figs",       showWarnings = FALSE)   # main-text figures (fig1-fig7)
dir.create("figs-appen", showWarnings = FALSE)   # appendix figures (figA-*)
dir.create("tables",  showWarnings = FALSE)

library(ggplot2)
library(data.table)
library(countrycode)
library(stringi)
library(uncounted)
library(gt)          # threeparttable LaTeX tables (tab_spanner, footnotes, as_latex)

source("codes/1-functions.R")             # shared helpers (gt_to_tex, register_diag_plots)
source("codes/2-prepare-data.R")          # read raw sources -> full_database_processed (+ data/*.csv)
source("codes/3-simulation-study.R")      # simulation study (Table 1)
source("codes/4-main-paper-analysis.R")   # modelling + main-text figures and tables
source("codes/5-supplement.R")            # appendix figures and tables
