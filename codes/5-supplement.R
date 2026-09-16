## Supplement (appendix) plots and tables. Uses `full_database_processed` from
## 2-prepare-data.R and the helpers (gt_to_tex, register_diag_plots) from 1-functions.R.

## ---- standalone setup -------------------------------------------------------
## When run on its own (not via 0-run-all.R), load packages, the processed data and the
## cached S8 fits / all-spec fits. Guarded so it is a no-op in the full pipeline.
## Helpers are sourced unconditionally so interactive re-runs always pick up the latest
## gt_to_tex()/register_diag_plots().
source("codes/1-functions.R")
if (!exists("full_database_processed")) {
  library(data.table)
  library(ggplot2)
  library(countrycode)
  library(stringi)
  library(uncounted)
  library(gt)
  full_database_processed <- fread("data/poland-for-model.csv")
  fit_all <- readRDS("results/fit-all.rds")     # all 8 specs, Poisson + NB
  fit_po  <- readRDS("results/fit-po-s8.rds")   # S8 Poisson PMLE
  fit_nb  <- readRDS("results/fit-nb-s8.rds")   # S8 NB-MLE
}

## ---- appendix modelling infrastructure --------------------------------------
## Three reference registers as modelling frames (m = border, n = police, N = register
## count) plus the S8 alpha/beta formulas. Used by the covariate-comparison, beta-
## selection, alpha-trajectory and coefficient sections of the appendix. Cheap to build,
## so defined unconditionally -- available whether run standalone or via 0-run-all.R.
registers <- c("ZUS", "PESEL")                                   # Tax omitted (restricted data)
ref_col   <- c(ZUS = "pop_insured", PESEL = "pop_register")
datasets  <- setNames(lapply(registers, function(r) {
  full_database_processed[get(ref_col[[r]]) > 0,
    .(year = as.factor(year), sex, country_code,
      m = border, n = police, N = get(ref_col[[r]]),
      ukr, simplified_proc, continent)]
}), registers)
f_alpha <- ~ year * ukr + sex          # S8 community anchor
f_beta  <- ~ year                      # S8 exposure elasticity
if (!exists("model_zus")) model_zus <- datasets$ZUS   # ZUS frame == main-text model_zus

## Table S9 Country-year cells with a zero reference count ($N = 0$) despite observed apprehensions ($m > 0$), for 18+ non-Schengen countries
cy <- full_database_processed[year %in% 2019:2024 & schengen == "non-Schengen",
                              .(m = sum(border,      na.rm = TRUE),
                                n = sum(police,      na.rm = TRUE),
                                PESEL = sum(pop_register, na.rm = TRUE),
                                ZUS   = sum(pop_insured,  na.rm = TRUE)),
                              by = .(year, country_code, country)]

gaps <- melt(cy, id.vars = c("year", "country_code", "country", "m", "n"),
             measure.vars = c("PESEL", "ZUS"),
             variable.name = "Register", value.name = "N")[N == 0 & m > 0]

## countries apprehended (m>0) but missing from the PESEL population register: only
## Hong Kong SAR (2023); ZUS has more gaps -- see table.
print(gaps[Register == "PESEL", .(year, country_code, country, m, n)])

gaps_tab <- gaps[, .(Countries = uniqueN(country_code), `Sum m` = sum(m), `Sum n` = sum(n)),
                 by = .(Register, Year = year)]
setorder(gaps_tab, Register, Year)

gap_gt <- gaps_tab |>
  gt(groupname_col = "Register") |>
  tab_caption(paste0("Country-year cells with a zero reference count (N = 0) despite observed ",
                     "apprehensions (m > 0), for 18+ non-Schengen countries.")) |>
  fmt_number(columns = c(Countries, `Sum m`, `Sum n`), decimals = 0) |>
  cols_label(Countries = "Number of Countries", `Sum m` = "Total apprehensions (Sum of $m$)",
             `Sum n` = "Sum $n$ (Police)") |>
  cols_align(align = "left",  columns = Year) |>
  cols_align(align = "right", columns = c(Countries, `Sum m`, `Sum n`)) |>
  tab_source_note(source_note = paste0("Note: m = border apprehensions; n = police records; ",
                                       "N = reference register count. Missing countries -- PESEL ",
                                       "register: Only Hong Kong (2023) is absent.")) |>
  tab_options(table.font.size = px(12), source_notes.font.size = px(12))

## label matches the manuscript (appendix, Model assumptions and details: Zero-count prevalence)
gt_to_tex(gap_gt, "tables/tblA0-zero-counts-assumptions.tex", label = "tbl-zero-counts-assumptions")
gap_gt


## Table S10 Zero counts in apprehension and police records by sex and year for 18+ non-Schengen countries with ZUS as a reference population

zc <- full_database_processed[year %in% 2019:2024 & schengen == "non-Schengen" & pop_insured > 0]

## one row per (year, country, sex) with a positive ZUS count, so within (year, sex) the
## number of rows equals the number of countries. The Total rows aggregate over sex to the
## country level first: a country counts as m = 0 (n = 0) only if it has no apprehensions
## (police records) of either sex, and the denominator is the number of countries with a
## positive ZUS count for at least one sex. Summing the sex-specific zero cells and dividing
## by the number of countries would mix cells and countries and give shares above 100%.
zc_sex <- zc[, .(m0 = sum(!(border > 0)), n0 = sum(!(police > 0)),
                 both = sum(border > 0 & police > 0),
                 Countries = uniqueN(country)), by = .(year, sex)]
zc_cty <- zc[, .(border = sum(border), police = sum(police)), by = .(year, country)]
zc_tot <- zc_cty[, .(m0 = sum(!(border > 0)), n0 = sum(!(police > 0)),
                     both = sum(border > 0 & police > 0),
                     Countries = uniqueN(country), sex = "total"), by = year]
stopifnot(zc_sex[, all(m0 <= Countries & n0 <= Countries & both <= Countries)],
          zc_tot[, all(m0 <= Countries & n0 <= Countries & both <= Countries)])

zc_all <- rbind(zc_sex, zc_tot)
zc_all[, ord := fifelse(sex == "m", 1L, fifelse(sex == "f", 2L, 3L))]
setorder(zc_all, year, ord)
zc_all[, Sex := factor(ord, 1:3, c("Males", "Females", "Total"))]
zc_all[, `:=`(m0_pct   = round(100 * m0   / Countries, 1),
              n0_pct   = round(100 * n0   / Countries, 1),
              both_pct = round(100 * both / Countries, 1))]
## Year printed once per block (on the Males row, ord == 1); blank on Females/Total
zc_tab <- zc_all[, .(Year = fifelse(ord == 1L, as.character(year), ""),
                     Sex, m0_pct, m0, n0_pct, n0, both_pct, both, Countries)]

zc_gt <- zc_tab |>
  gt() |>
  tab_caption(paste0("Zero counts in apprehension and police records by sex and year ",
                     "for 18+ non-Schengen countries.")) |>
  cols_label(m0_pct   = "m=0 (%)",    m0   = "m=0 (C)",
             n0_pct   = "n=0 (%)",    n0   = "n=0 (C)",
             both_pct = "Both>0 (%)", both = "Both>0 (C)") |>
  fmt_number(columns = c(m0_pct, n0_pct, both_pct), decimals = 1) |>
  fmt_number(columns = c(m0, n0, both, Countries),  decimals = 0) |>
  cols_align(align = "left",  columns = c(Year, Sex)) |>
  cols_align(align = "right", columns = c(m0_pct, m0, n0_pct, n0, both_pct, both, Countries)) |>
  tab_style(style = cell_text(style = "italic"),
            locations = cells_body(rows = Sex == "Total")) |>
  tab_source_note(source_note = paste0("Note: $m$ denotes border apprehensions; $n$ denotes ",
                                       "police records. Reference population: ZUS. C = number of ",
                                       "countries. The last column denotes the number of countries ",
                                       "with a positive ZUS count for a given sex in a given year. ",
                                       "Total rows aggregate $m$ and $n$ over sex within a country, ",
                                       "so a country counts as $m=0$ ($n=0$) only if it has no ",
                                       "apprehensions (police records) of either sex; the last ",
                                       "column is then the number of countries with a positive ZUS ",
                                       "count for at least one sex. Country includes also categories ",
                                       "such as stateless or unknown.")) |>
  tab_options(table.font.size = px(12), source_notes.font.size = px(12))

gt_to_tex(zc_gt, "tables/tblA1-zero-counts.tex", label = "tbl-zero-counts")
zc_gt


# Section Year and sex -- PESEL register (Tax register omitted -- restricted data)
register_diag_plots("pop_register", "PESEL",
                    "figs-appen/figA-pesel-alpha.pdf", "figs-appen/figA-pesel-beta.pdf",
                     facet = facet_wrap(~year, nrow = 2),
                     extra = scale_color_brewer(type = "qual", palette = "Set1"))

# Section Year and continent (ZUS and PESEL registers; Tax omitted -- restricted data)
register_diag_plots("pop_insured",
                    "ZUS",
                    "figs-appen/figA-zus-alpha-continent.pdf",
                    "figs-appen/figA-zus-beta-continent.pdf",
                     facet = facet_wrap(~continent, nrow = 2),
                     color_var = "year",
                     subset = !continent %in% c("unknown/stateless", "Oceania"),
                     extra = scale_color_brewer(type = "qual", palette = "Paired"))

register_diag_plots("pop_register",
                    "PESEL",
                    "figs-appen/figA-pesel-alpha-continent.pdf",
                    "figs-appen/figA-pesel-beta-continent.pdf",
                     facet = facet_wrap(~continent, nrow = 2),
                     color_var = "year",
                     subset = !continent %in% c("unknown/stateless", "Oceania"),
                     extra = scale_color_brewer(type = "qual", palette = "Paired"))


## ============================================================================
## Full covariate comparison (across the three reference registers)
## Cross-register fits are expensive -> cached to results/; set recompute_appendix
## <- TRUE (or delete the .rds) to refit.
## ============================================================================
recompute_appendix <- FALSE

## fit one (register, method, alpha, beta); NULL on failure
fitm <- function(d, method, a, b) tryCatch(
  estimate_hidden_pop(data = d, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
                      method = method, cov_alpha = a, cov_beta = b, countries = ~ country_code),
  error = function(e) NULL)
## fitted alpha per observation = X_alpha %*% alpha_coefs (original scale)
alpha_vals <- function(fit) as.numeric(fit$X_alpha %*% fit$alpha_coefs)

specs <- list(
  S1 = list(a = ~ 1,                      b = ~ 1,    label = "Intercept only"),   # time-invariant baseline: constant beta too (as in 3-main-paper-analysis.R)
  S2 = list(a = ~ year,                   b = ~ year, label = "Year"),
  S3 = list(a = ~ year + sex,             b = ~ year, label = "Year + sex"),
  S4 = list(a = ~ year * ukr,             b = ~ year, label = "Year x UKR"),
  S5 = list(a = ~ year * simplified_proc, b = ~ year, label = "Year x simplified"),
  S6 = list(a = ~ year + continent,       b = ~ year, label = "Year + continent"),
  S7 = list(a = ~ year + sex + ukr,       b = ~ year, label = "Year + sex + UKR"),
  S8 = list(a = ~ year * ukr + sex,       b = ~ year, label = "Year x UKR + sex"))

## ---- Table: AIC/BIC/QAIC across eight specifications and three registers ---------
## QAIC uses, within each register, the Pearson dispersion c-hat of the most general Poisson
## PMLE specification (S8); see c_hat()/qaic() in 1-functions.R. A cache written before the
## QAIC columns existed is rebuilt.
aic_cache_ok <- file.exists("results/appendix-aic-full.rds") &&
  "QAIC_PO" %in% names(readRDS("results/appendix-aic-full.rds"))
if (!recompute_appendix && aic_cache_ok) {
  aic_full <- readRDS("results/appendix-aic-full.rds")
} else {
  aic_full <- rbindlist(lapply(registers, function(reg) {
    fits <- lapply(specs, function(sp) list(po = fitm(datasets[[reg]], "poisson", sp$a, sp$b),
                                            nb = fitm(datasets[[reg]], "nb",      sp$a, sp$b)))
    chat_reg <- if (!is.null(fits$S8$po)) c_hat(fits$S8$po) else NA_real_
    rbindlist(lapply(names(specs), function(sn) {
      po <- fits[[sn]]$po; nb <- fits[[sn]]$nb
      data.table(Register = reg, Spec = sn, Label = specs[[sn]]$label,
                 AIC_PO   = if (!is.null(po)) round(AIC(po), 1) else NA_real_,
                 BIC_PO   = if (!is.null(po)) round(BIC(po), 1) else NA_real_,
                 QAIC_PO  = if (!is.null(po)) round(qaic(po, chat_reg), 1) else NA_real_,
                 chat_PO  = if (!is.null(po)) round(c_hat(po), 2) else NA_real_,
                 chat_ref = round(chat_reg, 2),
                 AIC_NB   = if (!is.null(nb)) round(AIC(nb), 1) else NA_real_,
                 BIC_NB   = if (!is.null(nb)) round(BIC(nb), 1) else NA_real_)
    }))
  }))
  saveRDS(aic_full, "results/appendix-aic-full.rds")
}

chat_note_aic <- aic_full[Spec == "S8", paste0(Register, ": ", sprintf("%.1f", chat_ref), collapse = "; ")]
aic_gt <- aic_full[, .(Spec, Label, AIC_PO, BIC_PO, QAIC_PO, chat_PO, AIC_NB, BIC_NB, Register)] |>
  gt(groupname_col = "Register") |>
  tab_caption(paste0("AIC, BIC and QAIC for eight covariate specifications (S1--S8) across three ",
                     "registers and two estimation methods.")) |>
  tab_spanner(label = "Poisson PMLE", columns = c(AIC_PO, BIC_PO, QAIC_PO, chat_PO)) |>
  tab_spanner(label = "NB-MLE",       columns = c(AIC_NB, BIC_NB)) |>
  cols_label(AIC_PO = "AIC", BIC_PO = "BIC", QAIC_PO = "QAIC", chat_PO = "$\\hat c$",
             AIC_NB = "AIC", BIC_NB = "BIC") |>
  fmt_number(columns = c(AIC_PO, BIC_PO, QAIC_PO, AIC_NB, BIC_NB), decimals = 1) |>
  fmt_number(columns = chat_PO, decimals = 2) |>
  tab_source_note(source_note = paste0("Note: $\\hat c$ = Pearson dispersion (chi-square over residual df) of the ",
                  "Poisson PMLE fit. QAIC $= -2\\ell/\\hat c_{S8} + 2(K+1)$, with $K$ the number of parameters ",
                  "(including $\\gamma$) and $\\hat c_{S8}$ the dispersion of the most general specification ",
                  "in each register (", chat_note_aic, "). Not defined for NB-MLE, which models the dispersion.")) |>
  tab_options(table.font.size = px(11))
gt_to_tex(aic_gt, "tables/tblA2-full-aic.tex", label = "tbl-full-aic")

## ---- Table: beta covariate selection (S8--S11), shared alpha = year x UKR + sex ----
beta_forms  <- list(S8 = ~ year, S9 = ~ year + sex, S10 = ~ year + ukr, S11 = ~ year + sex + ukr)
beta_labels <- c(S8 = "year", S9 = "year + sex", S10 = "year + UKR", S11 = "year + sex + UKR")
## QAIC (Poisson rows only) uses, within each register, the Pearson dispersion of the adopted
## specification S8. S11 is the most general model of the set, but its lower dispersion partly
## reflects an inadmissible Ukrainian detection rate (see tbl-s11-identification), so S8 is the
## cleaner yardstick. A cache written with a different c-hat convention is rebuilt.
beta_cache_ok <- file.exists("results/appendix-beta-comp.rds") &&
  "chat_source" %in% names(readRDS("results/appendix-beta-comp.rds"))
if (!recompute_appendix && beta_cache_ok) {
  beta_comp <- readRDS("results/appendix-beta-comp.rds")
} else {
  beta_comp <- rbindlist(lapply(registers, function(reg) {
    fits <- lapply(beta_forms, function(bf) lapply(c(poisson = "poisson", nb = "nb"), function(meth)
      fitm(datasets[[reg]], meth, ~ year * ukr + sex, bf)))
    chat_reg <- if (!is.null(fits$S8$poisson)) c_hat(fits$S8$poisson) else NA_real_
    rbindlist(lapply(names(beta_forms), function(bn) rbindlist(lapply(c("poisson", "nb"), function(meth) {
      fit <- fits[[bn]][[meth]]
      if (is.null(fit)) return(NULL)
      data.table(Register = reg, Spec = bn, `Covariates in beta` = beta_labels[bn],
                 Model = c(poisson = "Poisson PMLE", nb = "NB-MLE")[meth],
                 AIC = round(AIC(fit), 1), BIC = round(BIC(fit), 1),
                 QAIC  = if (meth == "poisson") round(qaic(fit, chat_reg), 1) else NA_real_,
                 c_hat = if (meth == "poisson") round(c_hat(fit), 2) else NA_real_,
                 chat_ref = round(chat_reg, 2), chat_source = "S8",
                 k = attr(logLik(fit), "df"))
    }))))
  }))
  setorder(beta_comp, Register, Model, AIC)
  saveRDS(beta_comp, "results/appendix-beta-comp.rds")
}

beta_comp[Model == "POISSON", Model := "Poisson PMLE"]   # normalise labels from older caches
beta_comp[Model == "NB",      Model := "NB-MLE"]
## order S8 -> S11 within each model (not by AIC), registers as ZUS/PESEL
beta_comp[, Spec     := factor(Spec, c("S8", "S9", "S10", "S11"))]
beta_comp[, Register := factor(Register, registers)]
setorder(beta_comp, Register, Model, Spec)
chat_note_beta <- beta_comp[Spec == "S8" & Model == "Poisson PMLE",
                            paste0(Register, ": ", sprintf("%.1f", chat_ref), collapse = "; ")]
beta_gt <- beta_comp[, .(Model, Spec, `Covariates in beta`, AIC, BIC, QAIC, c_hat, k, Register)] |>
  gt(groupname_col = "Register") |>
  tab_caption(paste0("Effect of beta specification on AIC, BIC and QAIC across three registers. S8 is the primary ",
                     "specification. All specifications share alpha = year x UKR + sex; ",
                     "k = number of parameters.")) |>
  cols_label(c_hat = "$\\hat c$") |>
  fmt_number(columns = c(AIC, BIC, QAIC), decimals = 1) |>
  fmt_number(columns = c_hat, decimals = 2) |>
  sub_missing(columns = c(QAIC, c_hat), missing_text = "") |>
  tab_source_note(source_note = paste0("Note: $\\hat c$ = Pearson dispersion of the Poisson PMLE fit; ",
                  "QAIC $= -2\\ell/\\hat c_{S8} + 2(k+1)$ with $\\hat c_{S8}$ the dispersion of the adopted ",
                  "specification S8 in each register (", chat_note_beta, "). ",
                  "Not defined for NB-MLE, which models the dispersion.")) |>
  tab_options(table.font.size = px(11))
gt_to_tex(beta_gt, "tables/tblA3-beta-comparison.tex", label = "tbl-beta-comparison")

## ---- Table: why S11 is not adopted -- identification of the Ukrainian decomposition ----
## S11 adds sex and UKR to beta. Ukraine already carries its own anchor in every year (alpha =
## year x UKR), so a UKR term in beta is identified from the same 12 Ukrainian cells: the fit
## improves but the size/detection split for Ukraine is no longer identified. The table
## contrasts S8 and S11 (Poisson PMLE, ZUS): fit, the beta:ukr coefficient, the range of the
## fitted Ukrainian detection rate (admissible only if <= 1), the cluster-FWB estimates for
## Ukrainians and non-Ukrainians, and the number of failed bootstrap replicates.
## The S11 bootstrap (R = 999, cluster by country, seed 2026) is cached like the S8 one.
fit_s11 <- fitm(model_zus, "poisson", ~ year * ukr + sex, ~ year + sex + ukr)
boot_s11_file <- "results/boot-ukr-po-s11.rds"
if (!recompute_appendix && file.exists(boot_s11_file)) {
  boot_s11 <- readRDS(boot_s11_file)
} else {
  boot_s11 <- bootstrap_popsize(fit_s11, by = ~ year + ukr, R = 999, cluster = ~ country_code,
                                level = 0.95, seed = 2026)
  saveRDS(boot_s11, boot_s11_file)
}
boot_s8 <- readRDS("results/boot-ukr-po.rds")          # S8, same settings (3-main-paper-analysis.R)
s11_row <- function(fit, boot, spec) {
  d   <- as.data.table(fit$data); rho <- fit$rho_values
  ps  <- as.data.table(boot$popsize); ps[, c("year", "ukr") := tstrsplit(group, ", ", fixed = TRUE)]
  ci  <- function(y, u) ps[year == y & ukr == u, sprintf("%s (%s, %s)", formatC(round(estimate / 1e3), format = "d", big.mark = ","),
                                                          formatC(round(pmax(lower, 0) / 1e3), format = "d", big.mark = ","),
                                                          formatC(round(upper / 1e3), format = "d", big.mark = ","))]
  b_ukr <- if ("beta:ukr" %in% names(coef(fit))) sprintf("%.2f (%.2f)", coef(fit)["beta:ukr"], sqrt(diag(fit$vcov))["beta:ukr"]) else "--"
  data.table(Spec = spec,
             QAIC = qaic(fit, c_hat(fit_po)),
             `beta: UKR (SE)` = b_ukr,
             `Ukrainian rho, range` = sprintf("%.2f--%.2f", min(rho[d$ukr == 1]), max(rho[d$ukr == 1])),
             `UKR 2019` = ci("2019", "1"), `UKR 2024` = ci("2024", "1"), `Non-UKR 2024` = ci("2024", "0"),
             `Failed rep.` = sum(!complete.cases(boot$t)))
}
s11_tab <- rbind(s11_row(fit_po, boot_s8, "S8"), s11_row(fit_s11, boot_s11, "S11"))
s11_gt <- s11_tab |>
  gt() |>
  tab_caption(paste0("S8 versus S11 (Poisson PMLE, ZUS reference): fit, the Ukrainian term in the exposure ",
                     "elasticity, the range of the fitted Ukrainian detection rate, and cluster-FWB ",
                     "estimates (median with 95% percentile CI, thousands) by origin.")) |>
  fmt_number(columns = QAIC, decimals = 1) |>
  cols_label(`beta: UKR (SE)` = "$\\hat\\beta_{\\text{UKR}}$ (SE)",
             `Ukrainian rho, range` = "$\\hat\\rho$, UKR cells",
             `UKR 2019` = "$\\hat\\xi$ UKR 2019", `UKR 2024` = "$\\hat\\xi$ UKR 2024",
             `Non-UKR 2024` = "$\\hat\\xi$ non-UKR 2024", `Failed rep.` = "Failed replicates") |>
  tab_source_note(source_note = paste0("Note: QAIC with $\\hat c$ from S8. SE = HC1 robust. $\\hat\\rho$ is admissible ",
                  "only if at most 1. Failed replicates = bootstrap draws (of 999) in which the fit did not ",
                  "converge. S8: $\\alpha \\sim$ year $\\times$ UKR + sex, $\\beta \\sim$ year; S11: same $\\alpha$, ",
                  "$\\beta \\sim$ year + sex + UKR.")) |>
  tab_options(table.font.size = px(10))
gt_to_tex(s11_gt, "tables/tblA-s11-identification.tex", label = "tbl-s11-identification")

## ---- Figure: S8 vs S11 cluster-FWB estimates by origin (Figure 7 styling) ----
## Same layout, palette (Set1) and CI encoding as Figure 7 in the main text: median with
## 80% (thick) and 95% (thin) percentile CIs, facets by Ukrainian origin, S8 and S11 dodged.
s11_boot_dt <- function(boot, spec) {
  ps  <- as.data.table(boot$popsize)
  q80 <- apply(boot$t, 2, quantile, probs = c(0.10, 0.90), na.rm = TRUE)
  ps[, `:=`(lower_80 = q80[1, ], upper_80 = q80[2, ], spec = spec)]
  ps[, c("year", "ukr") := tstrsplit(group, ", ", fixed = TRUE)]
  ps[]
}
s11_plot_dt <- rbind(s11_boot_dt(boot_s8, "S8"), s11_boot_dt(boot_s11, "S11"))
s11_plot_dt[, ukr_label := factor(fifelse(ukr == "1", "Ukraine", "Non-Ukraine"), c("Ukraine", "Non-Ukraine"))]
s11_plot_dt[, spec := factor(spec, c("S8", "S11"))]
pd_s11 <- position_dodge(width = 0.5)
ggplot(s11_plot_dt, aes(x = as.factor(year), y = estimate, colour = spec)) +
  geom_linerange(aes(ymin = pmax(lower, 0),    ymax = upper,    linewidth = "95%"), position = pd_s11) +
  geom_linerange(aes(ymin = pmax(lower_80, 0), ymax = upper_80, linewidth = "80%"), position = pd_s11) +
  geom_point(position = pd_s11, size = 2, shape = 21, fill = "white", stroke = 0.7) +
  facet_wrap(~ ukr_label, scales = "free_y") +
  scale_linewidth_manual(name = "CI", breaks = c("80%", "95%"), values = c("80%" = 1.6, "95%" = 0.5)) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "Year", y = expression(hat(xi)), colour = NULL) +
  scale_color_brewer(type = "qual", palette = "Set1") +
  guides(colour = guide_legend(order = 1), linewidth = guide_legend(order = 2)) +   # fixed legend order (reproducible output)
  theme(text = element_text(size = 15)) -> pA_s11
ggsave(plot = pA_s11, filename = "figs-appen/figA-s11-ukr-estimates.pdf", width = 10, height = 5)

## ---- Figure: all specifications S1-S11 by year and origin, S8 as reference ----
## One panel per specification (3 columns), Poisson PMLE cluster-FWB median with 80% (thick)
## and 95% (thin) percentile CIs by Ukrainian origin (Set1 colours as in Figure 7), with the
## S8 median (line) and 80% CI (band) repeated in every panel as the reference. Log scale.
## Bootstraps: R = 999, cluster by country, seed 2026 -- the paper's settings. S8 and S11
## reuse the caches above; the other nine are cached as results/boot-ukr-po-<spec>.rds.
all_specs <- c(lapply(specs, function(sp) list(a = sp$a, b = sp$b)),
               list(S9  = list(a = ~ year * ukr + sex, b = ~ year + sex),
                    S10 = list(a = ~ year * ukr + sex, b = ~ year + ukr),
                    S11 = list(a = ~ year * ukr + sex, b = ~ year + sex + ukr)))
spec_labels <- c(
  S1  = "S1: alpha = intercept only, beta = constant",
  S2  = "S2: alpha = year, beta = year",
  S3  = "S3: alpha = year + sex, beta = year",
  S7  = "S7: alpha = year + sex + UKR, beta = year",
  S4  = "S4: alpha = year x UKR, beta = year",
  S5  = "S5: alpha = year x simplified, beta = year",
  S6  = "S6: alpha = year + continent, beta = year",
  S8  = "S8: alpha = year x UKR + sex, beta = year",
  S9  = "S9: alpha = year x UKR + sex, beta = year + sex",
  S10 = "S10: alpha = year x UKR + sex, beta = year + UKR",
  S11 = "S11: alpha = year x UKR + sex, beta = year + sex + UKR")   # nested order: alpha chain, then beta chain
boot_all <- lapply(names(spec_labels), function(sn) {
  f <- switch(sn, S8 = "results/boot-ukr-po.rds", S11 = "results/boot-ukr-po-s11.rds",
              sprintf("results/boot-ukr-po-%s.rds", sn))
  if (!recompute_appendix && file.exists(f)) return(readRDS(f))
  fit <- fitm(model_zus, "poisson", all_specs[[sn]]$a, all_specs[[sn]]$b)
  b <- bootstrap_popsize(fit, by = ~ year + ukr, R = 999, cluster = ~ country_code, level = 0.95, seed = 2026)
  saveRDS(b, f); b
})
names(boot_all) <- names(spec_labels)
## strip labels: formulas on line 1, AIC and BIC of the Poisson PMLE fit (rounded) on line 2
fits_all <- lapply(names(spec_labels), function(sn) fitm(model_zus, "poisson", all_specs[[sn]]$a, all_specs[[sn]]$b))
ic_all   <- vapply(fits_all, function(f) sprintf("AIC = %s, BIC = %s", formatC(round(AIC(f)), format = "d", big.mark = ","),
                                                  formatC(round(BIC(f)), format = "d", big.mark = ",")), character(1))
strip_labels <- setNames(paste0(spec_labels, "\n", ic_all), names(spec_labels))
all_dt <- rbindlist(lapply(names(boot_all), function(sn) s11_boot_dt(boot_all[[sn]], sn)))
all_dt[, origin := factor(fifelse(ukr == "1", "Ukraine", "Non-Ukraine"), c("Ukraine", "Non-Ukraine"))]
all_dt[, spec := factor(spec, names(strip_labels), strip_labels)]
all_dt[, year := as.integer(year)]
all_dt[, `:=`(est_k = estimate / 1e3, l95 = pmax(lower, 1) / 1e3, u95 = upper / 1e3,
              l80 = pmax(lower_80, 1) / 1e3, u80 = upper_80 / 1e3)]
ref_s8 <- all_dt[spec == strip_labels["S8"], .(year, origin, s8 = est_k, s8_l80 = l80, s8_u80 = u80)]
all_dt <- merge(all_dt, ref_s8, by = c("year", "origin"))
pd_all <- position_dodge(width = 0.45)
ggplot(all_dt, aes(x = year)) +
  geom_ribbon(aes(ymin = s8_l80, ymax = s8_u80, group = origin), fill = "grey80", alpha = 0.5) +
  geom_line(aes(y = s8, group = origin, linetype = origin), colour = "grey40") +
  geom_linerange(aes(ymin = l95, ymax = u95, colour = origin, linewidth = "95%"), position = pd_all) +
  geom_linerange(aes(ymin = l80, ymax = u80, colour = origin, linewidth = "80%"), position = pd_all) +
  geom_point(aes(y = est_k, colour = origin), shape = 21, fill = "white", size = 1.8, stroke = 0.7, position = pd_all) +
  facet_wrap(~ spec, ncol = 3) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_x_continuous(breaks = 2019:2024) +
  scale_colour_brewer(type = "qual", palette = "Set1", name = NULL) +
  scale_linetype_manual(values = c(Ukraine = "solid", `Non-Ukraine` = "22"), name = "S8 reference (median, 80% CI)") +
  scale_linewidth_manual(name = "CI", breaks = c("80%", "95%"), values = c("80%" = 1.6, "95%" = 0.5)) +
  labs(x = NULL, y = expression(hat(xi) ~ "(thousands, log scale)")) +
  guides(colour = guide_legend(order = 1), linewidth = guide_legend(order = 2), linetype = guide_legend(order = 3)) +
  theme(text = element_text(size = 13), legend.position = "bottom",
        strip.text = element_text(size = 8.5, lineheight = 0.95)) -> pA_all_specs
ggsave(plot = pA_all_specs, filename = "figs-appen/figA-all-specs-ukr-estimates.pdf", width = 12, height = 12)

## ---- Figure: UKR vs non-UKR alpha trajectories across S4, S8, S11 and registers ----
specs_app <- list(S4  = list(a = ~ year * ukr,       b = ~ year),
                  S8  = list(a = ~ year * ukr + sex, b = ~ year),
                  S11 = list(a = ~ year * ukr + sex, b = ~ year + sex + ukr))
if (!recompute_appendix && file.exists("results/appendix-alpha-traj.rds")) {
  alpha_traj <- readRDS("results/appendix-alpha-traj.rds")
} else {
  alpha_traj <- rbindlist(lapply(names(specs_app), function(sn) {
    sp <- specs_app[[sn]]
    rbindlist(lapply(registers, function(reg) rbindlist(lapply(c("poisson", "nb"), function(meth) {
      fit <- fitm(datasets[[reg]], meth, sp$a, sp$b); if (is.null(fit)) return(NULL)
      d <- copy(fit$data); d[, alpha := alpha_vals(fit)]
      d[, .(alpha = mean(alpha)), by = .(Year = as.integer(as.character(year)), ukr)][
        , `:=`(Spec = sn, Model = toupper(meth), Register = reg,
               Group = fifelse(ukr == 1, "Ukraine", "Non-Ukraine"))]
    }))))
  }))
  saveRDS(alpha_traj, "results/appendix-alpha-traj.rds")
}
alpha_traj[Model == "POISSON", Model := "Poisson PMLE"]   # normalise labels from older caches
alpha_traj[Model == "NB",      Model := "NB-MLE"]
alpha_traj[, Spec     := factor(Spec, c("S4", "S8", "S11"))]
alpha_traj[, Register := factor(Register, registers)]
alpha_traj[, Model    := factor(Model, c("Poisson PMLE", "NB-MLE"))]   # Poisson on top row

ggplot(alpha_traj, aes(Year, alpha, colour = Group, linetype = Group)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
  facet_grid(Model ~ Spec + Register) +
  scale_colour_brewer(type = "qual", palette = "Set1") +
  labs(y = expression(alpha), colour = "", linetype = "") +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) -> pA_alpha
ggsave(pA_alpha, filename = "figs-appen/figA-alpha-ukr-appendix.pdf", width = 14, height = 8)

## ---- Table: S8 coefficient estimates (Poisson PML vs unconstrained NB) -------
## Bootstrap medians with 95% percentile CIs, computed from the per-replicate refit
## coefficients (boot_params) stored in the cached cluster-FWB S8 bootstraps
## (R = 999, clustered by country) -- the same bootstrap used for the population sizes.
boot_coef <- function(file) {
  bp  <- readRDS(file)$boot_params
  med <- apply(bp, 2, median,   na.rm = TRUE)
  lo  <- apply(bp, 2, quantile, 0.025, na.rm = TRUE)
  hi  <- apply(bp, 2, quantile, 0.975, na.rm = TRUE)
  setNames(sprintf("%.3f (%.3f, %.3f)", med, lo, hi), colnames(bp))
}
cp <- boot_coef("results/boot-ukr-po.rds")   # S8 Poisson, cluster-FWB
cn <- boot_coef("results/boot-ukr-nb.rds")   # S8 NB, cluster-FWB
amap <- c(
  "alpha:(Intercept)" = "Intercept (2019, female, non-UKR)",
  "alpha:year2020" = "Year 2020", "alpha:year2021" = "Year 2021", "alpha:year2022" = "Year 2022",
  "alpha:year2023" = "Year 2023", "alpha:year2024" = "Year 2024",
  "alpha:ukr" = "Ukrainian origin", "alpha:sexm" = "Male",
  "alpha:year2020:ukr" = "Year 2020 $\\times$ Ukrainian", "alpha:year2021:ukr" = "Year 2021 $\\times$ Ukrainian",
  "alpha:year2022:ukr" = "Year 2022 $\\times$ Ukrainian", "alpha:year2023:ukr" = "Year 2023 $\\times$ Ukrainian",
  "alpha:year2024:ukr" = "Year 2024 $\\times$ Ukrainian")
bmap <- c(
  "beta:(Intercept)" = "Intercept (2019)", "beta:year2020" = "Year 2020", "beta:year2021" = "Year 2021",
  "beta:year2022" = "Year 2022", "beta:year2023" = "Year 2023", "beta:year2024" = "Year 2024")
mk_block <- function(map, blk) data.table(block = blk, Term = unname(map),
                                          Poisson = cp[names(map)], NB = cn[names(map)])
## QAIC and the Pearson dispersion c-hat are reported for the Poisson PMLE only (c-hat of
## the S8 fit itself; NB-MLE models the dispersion).
gof <- data.table(block = "Fit statistics",
  Term = c("Num.Obs.", "AIC", "BIC", "QAIC", "Dispersion ($\\hat c$)", "Log.Lik."),
  Poisson = c(format(nobs(fit_po), big.mark = ","), sprintf("%.1f", AIC(fit_po)),
              sprintf("%.1f", BIC(fit_po)), sprintf("%.1f", qaic(fit_po, c_hat(fit_po))),
              sprintf("%.2f", c_hat(fit_po)), sprintf("%.3f", as.numeric(logLik(fit_po)))),
  NB      = c(format(nobs(fit_nb), big.mark = ","), sprintf("%.1f", AIC(fit_nb)),
              sprintf("%.1f", BIC(fit_nb)), "", "", sprintf("%.3f", as.numeric(logLik(fit_nb)))))
coef_tab <- rbind(mk_block(amap, "Community anchor ($\\hat{\\alpha}$)"),
                  mk_block(bmap, "Exposure elasticity ($\\hat{\\beta}$)"), gof)
coef_gt <- coef_tab |>
  gt(groupname_col = "block") |>
  tab_caption(paste0("S8 coefficient estimates under Poisson PML and unconstrained NB. ",
                     "Cluster fractional-weighted bootstrap medians with 95% percentile ",
                     "confidence intervals in parentheses (R = 999, clustered by country). ",
                     "ZUS reference population. Both models unconstrained.")) |>
  cols_label(Term = "") |>
  tab_options(table.font.size = px(11))
gt_to_tex(coef_gt, "tables/tblA4-coefficients.tex", label = "tbl-coefficients-appendix")

## ---- Table: alpha values per (year, origin, sex) cell, three specifications --
if (!recompute_appendix && file.exists("results/fit-nbc-s8.rds")) {
  fit_nbc <- readRDS("results/fit-nbc-s8.rds")
} else {
  fit_nbc <- estimate_hidden_pop(data = model_zus, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
    method = "nb", constrained = TRUE, cov_alpha = f_alpha, cov_beta = f_beta, countries = ~ country_code)
  saveRDS(fit_nbc, "results/fit-nbc-s8.rds")
}
get_alpha_cells <- function(fit, label) {
  eta <- as.numeric(fit$X_alpha %*% fit$alpha_coefs)
  av  <- if (isTRUE(fit$constrained)) plogis(eta) else eta
  dt  <- copy(fit$data); dt[, alpha := av]
  dt  <- unique(dt[, .(year, ukr, sex, alpha)]); dt[, Model := label]; dt
}
cells <- rbind(get_alpha_cells(fit_po, "Poisson"),
               get_alpha_cells(fit_nb, "NB"),
               get_alpha_cells(fit_nbc, "NB (constrained)"))
cells[, Origin := fifelse(ukr == 1L, "UKR", "Non-UKR")]
cells[, Sex    := factor(sex, c("f", "m"), c("Females", "Males"))]
cells_wide <- dcast(cells, year + Origin + Sex ~ Model, value.var = "alpha")
setorder(cells_wide, year, -Origin, Sex)
cells_wide[, Year := as.character(year)][, year := NULL]
setcolorder(cells_wide, c("Year", "Origin", "Sex", "NB", "NB (constrained)", "Poisson"))

cells_gt <- cells_wide |>
  gt() |>
  tab_caption(paste0("Alpha values per (year, Ukrainian origin, sex) cell under three ",
                     "specifications. Constrained NB back-transformed via logit-inverse.")) |>
  fmt_number(columns = c(NB, `NB (constrained)`, Poisson), decimals = 4) |>
  tab_options(table.font.size = px(11))
gt_to_tex(cells_gt, "tables/tblA5-alpha-cells.tex", label = "tbl-alpha-cells")


## ---- Figures: population size and share of reference population, by register --------
## S8 popsize by year x ukr for each register and method, via cluster fractional-weighted
## bootstrap (R = 199, clustered by country): point estimate = bootstrap median, 95% percentile
## CI. Cached -- set recompute_appendix <- TRUE (or delete the .rds) to regenerate.
if (!recompute_appendix && file.exists("results/appendix-ps-register-boot.rds")) {
  ps_register <- readRDS("results/appendix-ps-register-boot.rds")
} else {
  ps_register <- rbindlist(lapply(registers, function(reg) rbindlist(lapply(c("poisson", "nb"), function(meth) {
    fit <- fitm(datasets[[reg]], meth, f_alpha, f_beta); if (is.null(fit)) return(NULL)
    b   <- bootstrap_popsize(fit, by = ~ year + ukr, R = 199, cluster = ~ country_code,
                             level = 0.95, seed = 2026)
    ps  <- as.data.table(b$popsize)
    ps[, c("year", "ukr") := tstrsplit(group, ", ", fixed = TRUE)]
    ps[, .(year = as.integer(year),
           ukr_label = fifelse(ukr == "1", "Ukrainian", "Non-Ukrainian"),
           register = reg,
           method = c(poisson = "Poisson PMLE", nb = "NB-MLE")[[meth]],
           estimate, lower, upper)]
  }))))
  saveRDS(ps_register, "results/appendix-ps-register-boot.rds")
}

## reference-population totals per register x year x origin (iterate `registers`, which is unnamed)
ref_totals <- rbindlist(lapply(registers, function(reg) {
  d <- datasets[[reg]]
  totals <- d[, .(N_total = sum(N)), by = .(year, ukr)]
  totals[, `:=`(register = reg, ukr_label = fifelse(ukr == 1L, "Ukrainian", "Non-Ukrainian"))]
  totals
}))
ref_totals[, year := as.integer(as.character(year))]

## share of the reference population (merge on plain character keys, then factor for display)
ps_share <- merge(copy(ps_register), ref_totals[, .(year, register, ukr_label, N_total)],
                  by = c("year", "register", "ukr_label"), all.x = TRUE)
ps_share[, `:=`(share    = 100 * estimate / N_total,
                share_lo = 100 * pmax(lower, 0) / N_total,
                share_hi = 100 * upper / N_total)]

lev_reg <- registers; lev_meth <- c("Poisson PMLE", "NB-MLE"); lev_ukr <- c("Ukrainian", "Non-Ukrainian")
ps_register[, `:=`(register = factor(register, lev_reg), method = factor(method, lev_meth),
                   ukr_label = factor(ukr_label, lev_ukr))]
ps_share[, `:=`(register = factor(register, lev_reg), method = factor(method, lev_meth),
                ukr_label = factor(ukr_label, lev_ukr), year = factor(year, sort(unique(year))))]

pd_reg <- position_dodge(width = 0.4)

## Figure 1: estimated unauthorised population (thousands)
ggplot(ps_register, aes(x = factor(year), y = estimate / 1e3, colour = register, shape = register)) +
  geom_pointrange(aes(ymin = pmax(lower, 0) / 1e3, ymax = upper / 1e3),
                  position = pd_reg, size = 0.4, linewidth = 0.5) +
  facet_grid(ukr_label ~ method) +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_color_brewer(type = "qual", palette = "Set1") +
  coord_cartesian(ylim = c(0, 1000)) +                 # truncate at 1 million
  labs(x = "Year", y = expression(hat(xi) ~ "(thousands)"),
       colour = "Register", shape = "Register") +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))  -> pA_popsize
ggsave(pA_popsize, filename = "figs-appen/figA-popsize-register.pdf", width = 10, height = 5)

## Figure 2: estimate as a share of the reference population
ggplot(ps_share, aes(x = year, y = share, colour = register, shape = register)) +
  geom_pointrange(aes(ymin = share_lo, ymax = share_hi),
                  position = pd_reg, size = 0.4, linewidth = 0.5) +
  facet_grid(ukr_label ~ method) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_color_brewer(type = "qual", palette = "Set1") +
  coord_cartesian(ylim = c(0, 100)) +                  # truncate at 100% of reference pop
  labs(x = "Year", y = expression(hat(xi) / N ~ "(%)"),
       colour = "Register", shape = "Register") +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) -> pA_share
ggsave(pA_share, filename = "figs-appen/figA-popsize-share.pdf", width = 10, height = 5)


## ============================================================================
## Sensitivity of S8 to the auxiliary detection count n
## (police total / police PESEL-ID / prison)
## ----------------------------------------------------------------------------
## The three datasets share identical m (border), N (ZUS) and community rows and
## differ ONLY in the auxiliary count n. Refitting the SAME S8 specification on each
## isolates the effect of the n source; because m and the rows are identical, AIC/BIC
## are directly comparable and any change in xi-hat is attributable to the n source.
## ============================================================================
## The bootstrap is expensive -> cached to results/ (set recompute_nsource <- TRUE).
recompute_nsource <- FALSE

## three ZUS-reference modelling datasets that differ only in the auxiliary count n
## (model_zus, n = police total, is defined in 3-main-paper-analysis.R; add the two variants)
model_zus_police <- full_database_processed[pop_insured > 0,
  .(year = as.factor(year), sex, country_code, m = border, n = police_id_yes,
    N = pop_insured, ukr, simplified_proc, continent)]
model_zus_prison <- full_database_processed[pop_insured > 0,
  .(year = as.factor(year), sex, country_code, m = border, n = prison,
    N = pop_insured, ukr, simplified_proc, continent)]

n_sources <- list("Police (total)"    = model_zus,         # n = police (police_id_yes + police_id_no)
                  "Police (PESEL ID)" = model_zus_police,  # n = police_id_yes (person had a PESEL number)
                  "Prison"            = model_zus_prison)  # n = prison

## S8 main-paper specification (see codes/3-main-paper-analysis.R)
s8_alpha <- ~ year * ukr + sex
s8_beta  <- ~ year

## fit Poisson PMLE + NB-MLE for every n source (NULL on failure)
fits_n <- lapply(n_sources, function(d)
  lapply(c(poisson = "poisson", nb = "nb"), function(meth) tryCatch(
    estimate_hidden_pop(data = d, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
                        method = meth, cov_alpha = s8_alpha, cov_beta = s8_beta,
                        countries = ~ country_code),
    error = function(e) NULL)))

## fit-quality comparison (AIC / BIC / QAIC) per n source and method. QAIC uses one c-hat
## for all three sources: the Pearson dispersion of the primary fit (police total).
chat_nsrc <- if (!is.null(fits_n[["Police (total)"]]$poisson))
  c_hat(fits_n[["Police (total)"]]$poisson) else NA_real_
cmp_ic <- rbindlist(lapply(names(fits_n), function(nm) {
  f <- fits_n[[nm]]
  data.table(n_source = nm,
    po_AIC  = if (!is.null(f$poisson)) AIC(f$poisson) else NA_real_,
    po_BIC  = if (!is.null(f$poisson)) BIC(f$poisson) else NA_real_,
    po_QAIC = if (!is.null(f$poisson)) qaic(f$poisson, chat_nsrc) else NA_real_,
    po_chat = if (!is.null(f$poisson)) c_hat(f$poisson) else NA_real_,
    nb_AIC  = if (!is.null(f$nb))      AIC(f$nb)      else NA_real_,
    nb_BIC  = if (!is.null(f$nb))      BIC(f$nb)      else NA_real_)
}))
## yearly point estimates and summed-yearly totals per n source and method
cmp_ps <- rbindlist(lapply(names(fits_n), function(nm)
  rbindlist(lapply(names(fits_n[[nm]]), function(meth) {
    fit <- fits_n[[nm]][[meth]]; if (is.null(fit)) return(NULL)
    ps <- tryCatch(as.data.table(popsize(fit, by = ~ year)), error = function(e) NULL)
    if (is.null(ps)) return(NULL)
    ps[, .(n_source = nm, method = c(poisson = "Poisson PMLE", nb = "NB-MLE")[meth],
           year = as.integer(as.character(group)), estimate)]
  }))))
cmp_tot <- dcast(cmp_ps[, .(total = sum(estimate)), by = .(n_source, method)],
                 n_source ~ method, value.var = "total")

## ---- bootstrap (199 cluster reps): median and 80% / 95% percentile CIs per group ----
## For each fit we draw R = 199 cluster (by country) fractional-weight bootstrap replicates
## and summarise the population-size distribution per group:
##   point  = bootstrap median (b$popsize$estimate),
##   80% CI = 10th / 90th percentiles of the replicates (b$t),
##   95% CI = 2.5th / 97.5th percentiles of the replicates (b$t).
## boot_grid() runs boot_tidy() over every (n source x method) fit for a grouping `by`.
## boot_grid_cached() memoises the result to results/ (recompute_nsource controls refresh).
boot_tidy <- function(fit, by) {
  b  <- bootstrap_popsize(fit, by = by, R = 199, cluster = ~ country_code, seed = 2026)
  ps <- as.data.table(b$popsize)
  qs <- apply(b$t, 2, quantile, probs = c(0.025, 0.10, 0.90, 0.975), na.rm = TRUE)
  ps[, `:=`(lo95 = qs[1, ], lo80 = qs[2, ], hi80 = qs[3, ], hi95 = qs[4, ])]
  ps[, .(group = as.character(group), estimate, lo80, hi80, lo95, hi95)]
}
boot_grid <- function(by) {
  out <- rbindlist(lapply(names(fits_n), function(nm)
    rbindlist(lapply(names(fits_n[[nm]]), function(meth) {
      fit <- fits_n[[nm]][[meth]]; if (is.null(fit)) return(NULL)
      r <- tryCatch(boot_tidy(fit, by), error = function(e) NULL); if (is.null(r)) return(NULL)
      r[, `:=`(n_source = nm, method = c(poisson = "Poisson PMLE", nb = "NB-MLE")[meth])]
    }))), use.names = TRUE)
  out[, n_source := factor(n_source, names(n_sources))]
  out[, method   := factor(method, c("Poisson PMLE", "NB-MLE"))]
  out[]
}
boot_grid_cached <- function(by, file) {
  if (!recompute_nsource && file.exists(file)) return(readRDS(file))
  g <- boot_grid(by); saveRDS(g, file); g
}

boot_year <- boot_grid_cached(~ year,       "results/nsource-boot-year.rds")
boot_ukr  <- boot_grid_cached(~ year + ukr, "results/nsource-boot-ukr.rds")
boot_sex  <- boot_grid_cached(~ year + sex, "results/nsource-boot-sex.rds")

## attach year and split the "year, ukr" / "year, sex" group labels for faceting
boot_year[, year := as.integer(group)]
boot_ukr[, c("year", "ukr") := tstrsplit(group, ", ", fixed = TRUE)]
boot_ukr[, `:=`(year = as.integer(year),
                ukr_label = factor(fifelse(ukr == "1", "Ukraine", "Non-Ukraine"),
                                   c("Ukraine", "Non-Ukraine")))]
boot_sex[, c("year", "sex") := tstrsplit(group, ", ", fixed = TRUE)]
boot_sex[, `:=`(year = as.integer(year),
                sex_label = factor(fifelse(sex == "f", "Females", "Males"),
                                   c("Females", "Males")))]

## ---- figures: thick bar = 80% CI, thin bar = 95% CI (percentile); colour = n source ----
pd_cmp <- position_dodge(width = 0.5)
cmp_layers <- list(
  geom_linerange(aes(ymin = pmax(lo95, 0), ymax = hi95, linewidth = "95% CI"), position = pd_cmp),
  geom_linerange(aes(ymin = pmax(lo80, 0), ymax = hi80, linewidth = "80% CI"), position = pd_cmp),
  geom_point(position = pd_cmp, size = 2, shape = 21, fill = "white", stroke = 0.7),
  scale_linewidth_manual(name = "CI", breaks = c("80% CI", "95% CI"),
                         values = c("80% CI" = 1.6, "95% CI" = 0.5)),
  scale_y_continuous(labels = scales::label_comma()),
  scale_color_brewer(type = "qual", palette = "Set1"),
  labs(x = "Year", y = expression(hat(xi)), colour = "Auxiliary count (n)"),
  theme(text = element_text(size = 15)))

## (1) by year, faceted by model
ggplot(boot_year, aes(as.factor(year), estimate, colour = n_source, group = n_source)) +
  cmp_layers + facet_wrap(~ method) -> pA_nsource_year
ggsave(pA_nsource_year, filename = "figs-appen/figA-nsource-year.pdf", width = 10, height = 5)

## (2) by year x Ukrainian origin (S8 alpha already interacts year * ukr)
ggplot(boot_ukr, aes(as.factor(year), estimate, colour = n_source, group = n_source)) +
  cmp_layers + facet_grid(rows = vars(method), cols = vars(ukr_label)) -> pA_nsource_ukr
ggsave(pA_nsource_ukr, filename = "figs-appen/figA-nsource-ukr.pdf", width = 10, height = 5)

## (3) by year x sex (sex is additive in the S8 alpha). Prison has many female zeros, so
## the female detection signal is weak -> expect wider / unstable female CIs, esp. for Prison.
ggplot(boot_sex, aes(as.factor(year), estimate, colour = n_source, group = n_source)) +
  cmp_layers + facet_grid(rows = vars(method), cols = vars(sex_label)) -> pA_nsource_sex
ggsave(pA_nsource_sex, filename = "figs-appen/figA-nsource-sex.pdf", width = 10, height = 5)

## ---- diagnostics: correlation of each n with m, and zero shares by sex ----
## (a) correlation of each auxiliary count with the response m, across communities
aux_dt <- full_database_processed[pop_insured > 0,
  .(m = border, police, police_id_yes, prison, sex)]
cor_sp <- cor(aux_dt[, .(m, police, police_id_yes, prison)], method = "spearman")
cor_pe <- cor(aux_dt[, .(m, police, police_id_yes, prison)], method = "pearson")
## (b) share of structural zeros in each auxiliary count, by sex (the prison/female issue)
zero_share <- aux_dt[, .(police = mean(police == 0), police_id_yes = mean(police_id_yes == 0),
                         prison = mean(prison == 0)), by = sex]

## ---- Table: auxiliary-count characteristics (corr with m; zero shares by sex) ----
nv <- function(dt, col) setNames(dt[[col]], as.character(dt$n_source))
zf <- function(x) sprintf("%.0f", 100 * x)
aux_tab <- data.table(
  src = c("Police (total)", "Police (PESEL ID)", "Prison"),
  r   = c(cor_pe["m", "police"], cor_pe["m", "police_id_yes"], cor_pe["m", "prison"]),
  rho = c(cor_sp["m", "police"], cor_sp["m", "police_id_yes"], cor_sp["m", "prison"]),
  z_f = c(zero_share[sex == "f", police], zero_share[sex == "f", police_id_yes], zero_share[sex == "f", prison]),
  z_m = c(zero_share[sex == "m", police], zero_share[sex == "m", police_id_yes], zero_share[sex == "m", prison]))
aux_tex <- c(
  "\\begin{table}[H]", "\\centering",
  paste0("\\caption{\\label{tbl-nsource-aux}Auxiliary detection counts used in the S8 ",
         "sensitivity check: Pearson and Spearman correlation with the apprehension count $m$ ",
         "across community--year--sex cells, and the share of cells with a zero auxiliary count by sex.}"),
  "\\begin{tabular}{lrrrr} \\toprule",
  "Auxiliary count ($n$) & Pearson $r(m,n)$ & Spearman $\\rho(m,n)$ & Zeros, Females (\\%) & Zeros, Males (\\%) \\\\ \\midrule",
  aux_tab[, paste0(src, " & ", sprintf("%.2f", r), " & ", sprintf("%.2f", rho), " & ", zf(z_f), " & ", zf(z_m), " \\\\")],
  "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(aux_tex, "tables/tblA-nsource-aux.tex")

## ---- Table: fit (AIC/BIC) and summed-yearly point total xi-hat per source ----
order_src <- c("Police (total)", "Police (PESEL ID)", "Prison")
poA <- nv(cmp_ic, "po_AIC"); poB <- nv(cmp_ic, "po_BIC")
poQ <- nv(cmp_ic, "po_QAIC"); poC <- nv(cmp_ic, "po_chat")
nbA <- nv(cmp_ic, "nb_AIC"); nbB <- nv(cmp_ic, "nb_BIC")
totPo <- nv(cmp_tot, "Poisson PMLE"); totNb <- nv(cmp_tot, "NB-MLE")
fmt_ic <- function(v) formatC(round(v), format = "d", big.mark = ",")
fmt_k  <- function(v) formatC(round(v / 1e3), format = "d", big.mark = ",")
fit_tex <- c(
  "\\begin{table}[H]", "\\centering",
  paste0("\\caption{\\label{tbl-nsource-fit}Model fit and estimated population size under S8 ",
         "with three auxiliary detection counts. AIC, BIC and QAIC for Poisson PMLE, AIC and BIC ",
         "for NB-MLE; $\\hat c$ is the Pearson dispersion of each Poisson PMLE fit and QAIC ",
         "$= -2\\ell/\\hat c + 2(K+1)$ uses the dispersion of the primary fit, police total ",
         sprintf("($\\hat c = %.1f$), for all three sources. ", chat_nsrc),
         "$\\hat\\xi$ is the sum of the yearly point estimates (in thousands). Lower NB-MLE ",
         "AIC/BIC reflects over-dispersion accommodation and is not grounds for preferring it ",
         "(Section~\\ref{sec-cov-selection}).}"),
  "\\begin{tabular}{lrrrrrr|rr} \\toprule",
  "& \\multicolumn{4}{c}{Poisson PMLE} & \\multicolumn{2}{c|}{NB-MLE} & \\multicolumn{2}{c}{$\\hat\\xi$ (000s)} \\\\",
  "Auxiliary count ($n$) & AIC & BIC & QAIC & $\\hat c$ & AIC & BIC & Po & NB \\\\ \\midrule",
  vapply(order_src, function(s) paste0(
    s, " & ",
    fmt_ic(poA[s]), " & ", fmt_ic(poB[s]), " & ", fmt_ic(poQ[s]), " & ", sprintf("%.1f", poC[s]), " & ",
    fmt_ic(nbA[s]), " & ", fmt_ic(nbB[s]), " & ",
    fmt_k(totPo[s]), " & ", fmt_k(totNb[s]), " \\\\"), character(1)),
  "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(fit_tex, "tables/tblA-nsource-fit.tex")


#Leave-one-out sensitivity by country of
#origin under the S8 specification: change in \(\hat\xi\) when each
#country is dropped, top 20 countries shown. Both Poisson and NB are
#reported; the influence patterns are similar in rank but the absolute
#magnitudes differ owing to NB's quadratic variance function (see the
#observation-level analysis below).

## ---- LOO sensitivity under S8, re-plotted in ggplot2 (Set1, same theme as the rest) ----
## compare_loo() returns $table (ordered by influence) with the percentage change in xi-hat
## when each country / observation is dropped, for Poisson and NB; we re-plot those values so
## the LOO figures share the colour scheme and theme of the rest of the paper.
loo_bar_gg <- function(cmp, n_top) {
  tab <- as.data.table(cmp$table)
  pc1 <- names(tab)[4]; pc2 <- names(tab)[6]            # pct_<model1>, pct_<model2>
  top <- head(tab, n_top)
  long <- rbind(data.table(label = top$label, Model = cmp$labels[1], pct = top[[pc1]]),
                data.table(label = top$label, Model = cmp$labels[2], pct = top[[pc2]]))
  long[, label := factor(label, rev(top$label))]       # most influential at the top
  long[, Model := factor(Model, cmp$labels)]
  ggplot(long, aes(pct, label, fill = Model)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    labs(x = expression("Change in " * hat(xi) * " when dropped (%)"), y = NULL, fill = NULL) +
    theme(text = element_text(size = 15), legend.position = "top")
}
loo_scatter_gg <- function(cmp, n_top) {
  tab <- as.data.table(cmp$table)
  pc1 <- names(tab)[4]; pc2 <- names(tab)[6]
  tab[, lab := fifelse(seq_len(.N) <= n_top, label, NA_character_)]
  ggplot(tab, aes(.data[[pc1]], .data[[pc2]])) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(colour = "#377EB8", alpha = 0.6) +
    geom_text(aes(label = lab), size = 3, vjust = -0.6, check_overlap = TRUE, na.rm = TRUE) +
    labs(x = paste0(cmp$labels[1], " influence (%)"),
         y = paste0(cmp$labels[2], " influence (%)")) +
    theme(text = element_text(size = 15))
}

## LOO is expensive (drops each country / observation and refits) -> cached to results/
if (!recompute_appendix && file.exists("results/appendix-loo.rds")) {
  loo_cmp <- readRDS("results/appendix-loo.rds")
  loo_comp_s8 <- loo_cmp$cty; loo_obs_comp <- loo_cmp$obs
} else {
  loo_comp_s8 <- compare_loo(loo(fit_po, by = "country", verbose = FALSE),
                             loo(fit_nb, by = "country", verbose = FALSE),
                             labels = c("Poisson S8", "NB S8"))
  loo_obs_comp <- compare_loo(loo(fit_po, by = "obs", verbose = FALSE),
                              loo(fit_nb, by = "obs", verbose = FALSE),
                              labels = c("Poisson S8", "NB S8"),
                              data = model_zus, label_vars = c("country_code", "year", "sex"))
  saveRDS(list(cty = loo_comp_s8, obs = loo_obs_comp), "results/appendix-loo.rds")
}

ggsave(plot = loo_bar_gg(loo_comp_s8, 20),      filename = "figs-appen/figA-loo-s8.pdf",          width = 10, height = 10)
ggsave(plot = loo_bar_gg(loo_obs_comp, 30),     filename = "figs-appen/figA-loo-obs.pdf",         width = 10, height = 10)
ggsave(plot = loo_scatter_gg(loo_obs_comp, 30), filename = "figs-appen/figA-loo-obs-scatter.pdf", width = 10, height = 5)


#Leave-one-out sensitivity by country of
#origin under the S8 specification: change in \(\hat\xi\) when each
#country is dropped, top 20 countries shown. Both Poisson and NB are
#reported; the influence patterns are similar in rank but the absolute
#magnitudes differ owing to NB's quadratic variance function (see the
#observation-level analysis below) 

## (observation-level LOO is computed in the cached block above and plotted via loo_bar_gg /
##  loo_scatter_gg; the country-level and observation-level figures are written there.)

## ---- Anscombe residuals vs fitted, in ggplot2 (paper theme; faceted by method) ----
r_po <- residuals(fit_po, type = "anscombe")
r_nb <- residuals(fit_nb, type = "anscombe")
x_po <- sqrt(fit_po$fitted.values)
x_nb <- sqrt(fit_nb$fitted.values)
res_dt <- rbind(data.table(Model = "Poisson S8", x = x_po, r = r_po),
                data.table(Model = "NB S8",      x = x_nb, r = r_nb))
res_dt[, Model := factor(Model, c("Poisson S8", "NB S8"))]
ggplot(res_dt, aes(x, r)) +
  geom_point(shape = 1, colour = adjustcolor("gray30", 0.4)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_smooth(method = "loess", se = FALSE, colour = "#377EB8", linewidth = 1) +
  facet_wrap(~ Model) +
  coord_cartesian(ylim = c(-25, 20)) +
  labs(x = expression(sqrt(hat(mu))), y = "Anscombe residual") +
  theme(text = element_text(size = 15)) -> p_resid
ggsave(plot = p_resid, filename = "figs-appen/figA-residuals.pdf", width = 10, height = 5)

## ---- Gamma profiles for S8, in ggplot2 (paper theme) ----
## profile_gamma(..., plot = FALSE) returns gamma, xi and loglik over the grid; the dashed
## red line marks the estimated gamma.
gamma_profile_gg <- function(fit) {
  pr <- as.data.table(profile_gamma(fit, gamma_grid = seq(1e-4, 0.3, length.out = 30), plot = FALSE))
  long <- melt(pr, id.vars = "gamma", measure.vars = c("xi", "loglik"))
  long[, variable := factor(variable, c("xi", "loglik"), c("Population size", "Log-likelihood"))]
  ggplot(long, aes(gamma, value)) +
    geom_line(colour = "#377EB8") +
    geom_point(colour = "#377EB8", size = 1) +
    geom_vline(xintercept = fit$gamma, linetype = 2, colour = "#E41A1C") +
    facet_wrap(~ variable, scales = "free_y") +
    labs(x = expression(gamma), y = NULL) +
    theme(text = element_text(size = 15))
}
ggsave(plot = gamma_profile_gg(fit_po), filename = "figs-appen/figA-fig-gamma-profile-po.pdf", width = 10, height = 5)
ggsave(plot = gamma_profile_gg(fit_nb), filename = "figs-appen/figA-fig-gamma-profile-nb.pdf", width = 10, height = 5)


## ---- Table: sensitivity to covariate-varying gamma --------------------------
## S8 (alpha = f_alpha, beta = f_beta) with gamma varying by different covariate sets,
## Poisson and NB, ZUS reference population. Expensive (14 fits) -> cached.
gamma_specs_cg <- list(
  "constant"      = NULL,
  "~year"         = ~ year,
  "~sex"          = ~ sex,
  "~ukr"          = ~ ukr,
  "~year+ukr"     = ~ year + ukr,
  "~year+sex"     = ~ year + sex,
  "~year+sex+ukr" = ~ year + sex + ukr)
gamma_labels <- c(
  "constant"      = "constant",
  "~year"         = "$\\gamma \\sim \\text{year}$",
  "~sex"          = "$\\gamma \\sim \\text{sex}$",
  "~ukr"          = "$\\gamma \\sim \\text{UKR}$",
  "~year+ukr"     = "$\\gamma \\sim \\text{year} + \\text{UKR}$",
  "~year+sex"     = "$\\gamma \\sim \\text{year} + \\text{sex}$",
  "~year+sex+ukr" = "$\\gamma \\sim \\text{year} + \\text{sex} + \\text{UKR}$")

## keep the fitted models (fits_cg) -- needed for the figures below
if (!recompute_appendix && file.exists("results/appendix-cov-gamma-fits.rds")) {
  fits_cg <- readRDS("results/appendix-cov-gamma-fits.rds")
} else {
  fits_cg <- list()
  for (gname in names(gamma_specs_cg)) for (meth in c("poisson", "nb")) {
    fit <- tryCatch(suppressWarnings(estimate_hidden_pop(
      data = model_zus, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
      method = meth, vcov = "HC1", cov_alpha = f_alpha, cov_beta = f_beta,
      gamma = "estimate", cov_gamma = gamma_specs_cg[[gname]],
      countries = ~ country_code)), error = function(e) NULL)
    if (!is.null(fit)) fits_cg[[paste(gname, meth, sep = "_")]] <- fit
  }
  saveRDS(fits_cg, "results/appendix-cov-gamma-fits.rds")
}

## summary stats derived from the fitted models. QAIC (Poisson only) uses the Pearson
## dispersion of the adopted specification (constant gamma = S8), the same yardstick as in
## Table 3 and the beta-comparison table.
stats_cg_dt <- rbindlist(lapply(names(fits_cg), function(key) {
  fit <- fits_cg[[key]]
  data.table(cov_gamma = sub("_(poisson|nb)$", "", key),
             method    = toupper(sub(".*_(poisson|nb)$", "\\1", key)),
             AIC = round(AIC(fit), 1), BIC = round(BIC(fit), 1),
             k_gamma = fit$p_gamma,   # number of gamma parameters (not distinct fitted values)
             c_hat = if (grepl("_poisson$", key)) c_hat(fit) else NA_real_,
             theta = if (!is.null(fit$theta)) fit$theta else NA_real_)
}))
chat_cg <- stats_cg_dt[method == "POISSON" & cov_gamma == "constant", c_hat]
stats_cg_dt[method == "POISSON",
            QAIC := round(sapply(paste0(cov_gamma, "_poisson"), function(key) qaic(fits_cg[[key]], chat_cg)), 1)]

## Poisson | NB side by side, ordered by the gamma-specification list
cg_po <- stats_cg_dt[method == "POISSON", .(cov_gamma, AIC_PO = AIC, BIC_PO = BIC, QAIC_PO = QAIC,
                                            chat_PO = round(c_hat, 2), k_PO = k_gamma)]
cg_nb <- stats_cg_dt[method == "NB",      .(cov_gamma, AIC_NB = AIC, BIC_NB = BIC, k_NB = k_gamma)]
cg_display <- merge(cg_po, cg_nb, by = "cov_gamma")
cg_display[, Spec := gamma_labels[cov_gamma]]
cg_display[, cov_gamma := factor(cov_gamma, names(gamma_specs_cg))]
setorder(cg_display, cov_gamma)
cg_display <- cg_display[, .(Spec, AIC_PO, BIC_PO, QAIC_PO, chat_PO, k_PO, AIC_NB, BIC_NB, k_NB)]

cg_gt <- cg_display |>
  gt() |>
  tab_caption(paste0("AIC, BIC and QAIC by gamma specification and estimation method. All models use ",
                     "the S8 specification for alpha and beta with the ZUS reference population.")) |>
  tab_spanner(label = "Poisson PMLE", columns = c(AIC_PO, BIC_PO, QAIC_PO, chat_PO, k_PO)) |>
  tab_spanner(label = "NB-MLE",       columns = c(AIC_NB, BIC_NB, k_NB)) |>
  cols_label(Spec = "$\\gamma$ specification",
             AIC_PO = "AIC", BIC_PO = "BIC", QAIC_PO = "QAIC", chat_PO = "$\\hat c$", k_PO = "$k_\\gamma$",
             AIC_NB = "AIC", BIC_NB = "BIC", k_NB = "$k_\\gamma$") |>
  fmt_number(columns = c(AIC_PO, BIC_PO, QAIC_PO, AIC_NB, BIC_NB), decimals = 1) |>
  fmt_number(columns = chat_PO, decimals = 2) |>
  tab_source_note(source_note = paste0("Note: S8 specification: $\\alpha \\sim \\text{year} ",
                  "\\times \\text{UKR} + \\text{sex}$, $\\beta \\sim \\text{year}$. ZUS reference ",
                  "population. $k_\\gamma$ = number of gamma parameters; $\\hat c$ = Pearson dispersion ",
                  "of the Poisson PMLE fit; QAIC $= -2\\ell/\\hat c_{S8} + 2(K+1)$ with $\\hat c_{S8}$ from the ",
                  sprintf("adopted specification, constant $\\gamma$ ($\\hat c_{S8} = %.1f$). Lower AIC/BIC/QAIC is better.", chat_cg))) |>
  tab_options(table.font.size = px(11))
gt_to_tex(cg_gt, "tables/tblA6-cov-gamma-aic.tex", label = "tbl-cov-gamma-aic")

## ---- Figure: profile log-likelihood of constant gamma (Poisson vs NB) -------
## profile_gamma(plot = FALSE) returns the grid (gamma, xi, loglik); we plot loglik.
prof_cg <- rbindlist(lapply(c("poisson", "nb"), function(meth) {
  d <- as.data.table(profile_gamma(fits_cg[[paste0("constant_", meth)]],
                                   gamma_grid = seq(0.001, 0.05, length.out = 50), plot = FALSE))
  d[, Method := c(poisson = "Poisson", nb = "NB")[[meth]]]
}))
prof_cg[, Method := factor(Method, c("Poisson", "NB"))]

ggplot(prof_cg, aes(gamma, loglik)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1) +
  facet_wrap(~ Method, scales = "free_y") +
  labs(x = expression(gamma), y = "Profile log-likelihood") +
  theme(text = element_text(size = 15)) -> pA_cg_profile
ggsave(pA_cg_profile, filename = "figs-appen/figA-cov-gamma-profile.pdf", width = 10, height = 5)

## ---- Figure: gamma coefficients (log scale) with 95% Wald CIs ----------------
gamma_coef_cg <- rbindlist(lapply(c("~year", "~ukr", "~year+ukr", "~year+sex+ukr"), function(gname)
  rbindlist(lapply(c("poisson", "nb"), function(meth) {
    fit <- fits_cg[[paste(gname, meth, sep = "_")]]
    if (is.null(fit) || is.null(fit$gamma_coefs)) return(NULL)
    gc <- fit$gamma_coefs
    p_ab <- fit$p_alpha + fit$p_beta
    gamma_idx <- (p_ab + 1):(p_ab + length(gc))
    se <- sqrt(diag(fit$vcov_full[gamma_idx, gamma_idx, drop = FALSE]))
    data.table(cov_gamma = gname, method = toupper(meth), term = names(gc),
               estimate = gc, lower = gc - 1.96 * se, upper = gc + 1.96 * se)
  }))))

ggplot(gamma_coef_cg, aes(x = term, y = estimate, colour = method)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.5), size = 0.4, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  facet_wrap(~ cov_gamma, scales = "free_x") +
  scale_colour_brewer(type = "qual", palette = "Set1") +
  labs(x = NULL, y = "Gamma coefficient (log scale)", colour = "Method") +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1)) -> pA_cg_coefs
ggsave(pA_cg_coefs, filename = "figs-appen/figA-cov-gamma-coefs.pdf", width = 12, height = 8)

## ---- Figure: estimated population by year, constant vs covariate-varying gamma ----
ps_cg_year <- rbindlist(lapply(c("constant", "~year", "~year+sex+ukr"), function(gname)
  rbindlist(lapply(c("poisson", "nb"), function(meth) {
    ps <- as.data.table(popsize(fits_cg[[paste(gname, meth, sep = "_")]], by = ~ year))
    ps[, .(year = as.integer(group), estimate, lower, upper,
           method = c(poisson = "Poisson", nb = "NB")[[meth]], cov_gamma = gname)]
  }))))
ps_cg_year[, cov_gamma := factor(cov_gamma, c("constant", "~year", "~year+sex+ukr"))]
ps_cg_year[, method    := factor(method, c("Poisson", "NB"))]
pd_cg <- position_dodge(width = 0.4)

ggplot(ps_cg_year, aes(x = factor(year), y = estimate, colour = method)) +
  geom_pointrange(aes(ymin = pmax(lower, 0), ymax = upper),
                  position = pd_cg, size = 0.4, linewidth = 0.5) +
  facet_wrap(~ cov_gamma) +
  scale_colour_brewer(type = "qual", palette = "Set1") +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "Year", y = expression(hat(xi)), colour = "Method") +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) -> pA_cg_popsize
ggsave(pA_cg_popsize, filename = "figs-appen/figA-cov-gamma-popsize.pdf", width = 12, height = 6)

## ---- Figure: estimated population by year x UKR, constant vs varying gamma ----
ps_cg_ukr <- rbindlist(lapply(c("constant", "~year+sex+ukr"), function(gname)
  rbindlist(lapply(c("poisson", "nb"), function(meth) {
    ps <- as.data.table(popsize(fits_cg[[paste(gname, meth, sep = "_")]], by = ~ year + ukr))
    ps[, .(year = as.integer(sub(",.*", "", group)),
           ukr_label = fifelse(grepl(",\\s*1$", group), "Ukrainian", "Non-Ukrainian"),
           estimate, lower, upper,
           method = c(poisson = "Poisson", nb = "NB")[[meth]], cov_gamma = gname)]
  }))))
ps_cg_ukr[, cov_gamma := factor(cov_gamma, c("constant", "~year+sex+ukr"))]
ps_cg_ukr[, method    := factor(method, c("Poisson", "NB"))]

ggplot(ps_cg_ukr, aes(x = factor(year), y = estimate, colour = method)) +
  geom_pointrange(aes(ymin = pmax(lower, 0), ymax = upper),
                  position = pd_cg, size = 0.35, linewidth = 0.5) +
  facet_grid(ukr_label ~ cov_gamma, scales = "free_y") +
  scale_colour_brewer(type = "qual", palette = "Set1") +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "Year", y = expression(hat(xi)), colour = "Method") +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) -> pA_cg_popsize_ukr
ggsave(pA_cg_popsize_ukr, filename = "figs-appen/figA-cov-gamma-popsize-ukr.pdf", width = 12, height = 8)


## ---- Figure: profile log-likelihood of the alpha and beta intercepts --------
## profile_alpha()/profile_beta() (plot = FALSE) return value, xi, loglik; we plot the
## relative log-likelihood of each intercept for the S8 Poisson and (unconstrained) NB fits.
prof_a_po <- as.data.table(profile_alpha(fit_po, plot = FALSE))
prof_a_nb <- as.data.table(profile_alpha(fit_nb, plot = FALSE))
prof_b_po <- as.data.table(profile_beta(fit_po,  plot = FALSE))
prof_b_nb <- as.data.table(profile_beta(fit_nb,  plot = FALSE))
prof_all <- rbind(
  prof_a_po[, .(param = "alpha intercept", model = "Poisson S8", value, loglik)],
  prof_a_nb[, .(param = "alpha intercept", model = "NB S8",      value, loglik)],
  prof_b_po[, .(param = "beta intercept",  model = "Poisson S8", value, loglik)],
  prof_b_nb[, .(param = "beta intercept",  model = "NB S8",      value, loglik)])
prof_all[, rel_ll := loglik - max(loglik, na.rm = TRUE), by = .(param, model)]

ggplot(prof_all, aes(x = value, y = rel_ll, colour = model)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ param, scales = "free_x") +
  coord_cartesian(ylim = c(-15, 1)) +
  scale_colour_brewer(type = "qual", palette = "Set1") +
  labs(x = "Coefficient value", y = expression(log * L - max * " " * log * L), colour = NULL) +
  theme(text = element_text(size = 15), legend.position = "bottom") -> pA_profile_ab
ggsave(pA_profile_ab, filename = "figs-appen/figA-profile-alpha-beta.pdf", width = 10, height = 4.5)


## ============================================================================
## Simulation-study tables (from results/simulation-results/, written by 4-simulation-study.R)
## ============================================================================
sim_dir <- "results/simulation-results"
sim_ml  <- c(nb = "NB", nb_c = "NB (c)", nls = "NLS", ols = "OLS",
             poisson = "Poisson", poisson_c = "Poisson (c)")
sim_f1  <- function(x) sprintf("%.1f", round(x, 1))
sim_rf  <- function(x) formatC(round(x), format = "d", big.mark = ",")
sim_f4  <- function(x) sprintf("%.4f", x)

## ---- Table: simulation summary (Setup-I and Setup-II, Po-PMLE vs NB-MLE) ----
## Simulation results, 2000 replications. Bias relative to $\xi = \sum_{i=1}^{100} \xi_i$. RMSE, root mean squared error. Coverage of $\xi$. Confidence interval (CI) of $\xi$ by \eqref{eq-xi-confint}, CI width relative to $\xi$.
## Built from the analytical-CI metrics written by 4-simulation-study.R:
##   Setup-I  (gamma = 0, regrouped 'rest' community) -> sim1_metrics_analytical.csv
##   Setup-II (gamma > 0, m_i n_i = 0 allowed)        -> sim2_metrics_analytical.csv
## Shown: Po(n)-Po(m) and NB(n)-NB(m) scenarios; estimators poisson (Po-PMLE) and nb (NB-MLE).
sim1_m   <- fread(file.path(sim_dir, "sim1_metrics_analytical.csv"))
sim2_m   <- fread(file.path(sim_dir, "sim2_metrics_analytical.csv"))

## one body row per (scenario, estimator), ordered Po(n)-Po(m) then NB(n)-NB(m)
sim_body <- function(dt) {
  d <- dt[scenario %in% c("Po(n)-Po(m)", "NB(n)-NB(m)") & model %in% c("poisson", "nb")]
  d[, Estimator := c(poisson = "Po-PMLE", nb = "NB-MLE")[model]]
  d[, scen_ord := match(scenario, c("Po(n)-Po(m)", "NB(n)-NB(m)"))]
  d[, est_ord  := match(model, c("poisson", "nb"))]
  setorder(d, scen_ord, est_ord)   # Po(n)-Po(m) before NB(n)-NB(m); Po-PMLE before NB-MLE
  d[, paste0(scenario, " & ", Estimator, " & ",
             sprintf("%.1f", round(`Bias (%)`, 1)), " & ",
             formatC(round(RMSE), format = "d", big.mark = ","), " & ",
             sprintf("%.1f", round(`Coverage (%)`, 1)), " & ",
             sprintf("%.1f", round(`CI width (%)`, 1)), " \\\\")]
}

sim_hdr <- "Data & Estimator & Bias (\\%) & RMSE & Coverage (\\%) & CI width (\\%)\\\\ \\midrule"
sim_tex <- c(
  "\\begin{table}[H]",
  "\\centering",
  paste0("\\caption{\\label{tbl-sim-summary}Simulation results, 2000 replications. Bias relative ",
         "to $\\xi = \\sum_{i=1}^{100} \\xi_i$. RMSE, root mean squared error. Coverage of $\\xi$. ",
         "Confidence interval (CI) of $\\xi$ by \\eqref{eq-xi-confint}, CI width relative to $\\xi$.}"),
  "\\begin{tabular}{llrrrr} \\toprule",
  "& \\multicolumn{5}{c}{Setup-I: $\\gamma =0$, regrouped `rest' community to void $m_i n_i =0$} \\\\ \\cline{2-6}",
  sim_hdr,
  sim_body(sim1_m),
  "\\bottomrule",
  "& \\multicolumn{5}{c}{Setup-II: $\\gamma >0$, communities with $m_i n_i =0$ allowed}\\\\ \\cline{2-6}",
  sim_hdr,
  sim_body(sim2_m),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(sim_tex, "tables/tblA-sim-summary.tex")

## ---- Table: Simulation 1 (aggregation, gamma = 0) ----
s1 <- fread(file.path(sim_dir, "sim1_metrics_analytical.csv"))[, Model := sim_ml[model]]
setorder(s1, scenario, model)
s1[, .(Scenario = scenario, Model, Bias = sim_f1(`Bias (%)`), BiasBC = sim_f1(`Bias BC (%)`),
       RMSE = sim_rf(RMSE), Cov = sim_f1(`Coverage (%)`), CIw = sim_f1(`CI width (%)`),
       Rok = as.character(R_ok))] |>
  gt(groupname_col = "Scenario") |>
  tab_caption(paste0("Simulation 1: bias, RMSE, coverage and CI width across four DGP scenarios under ",
                     "aggregation (all $m>0$, $n>0$). Analytical (delta-method) CIs. $R = 2{,}000$ replications.")) |>
  cols_label(Bias = "Bias (%)", BiasBC = "Bias BC (%)", RMSE = "RMSE", Cov = "Coverage (%)",
             CIw = "CI width (%)", Rok = "$R_{\\text{ok}}$") |>
  cols_align("right", columns = c(Bias, BiasBC, RMSE, Cov, CIw, Rok)) |>
  tab_source_note(paste0("Note: bias and CI width as a percentage of true $\\xi$; RMSE on the original ",
                         "scale; coverage is the share of replications whose 95% analytical CI contains ",
                         "$\\xi$; $R_{\\text{ok}}$ is the number of finite replications out of 2,000. ",
                         "NLS bias, bias BC and RMSE use median-based summaries (root-median-squared error) ",
                         "because a few replicates are numerically extreme under NB-generated data; the ",
                         "other estimators use mean-based summaries.")) |>
  tab_options(table.font.size = px(10), source_notes.font.size = px(9)) |>
  gt_to_tex("tables/tblA-sim1.tex", label = "tbl-sim1")

## ---- Table: Simulation 2, analytical CI (zeros, gamma estimated) ----
s2 <- fread(file.path(sim_dir, "sim2_metrics_analytical.csv"))[, Model := sim_ml[model]]
setorder(s2, scenario, model)
s2[, .(Scenario = scenario, Model, Bias = sim_f1(`Bias (%)`), BiasBC = sim_f1(`Bias BC (%)`),
       Biasmed = sim_f1(`Bias median (%)`), RMSE = sim_rf(RMSE), Cov = sim_f1(`Coverage (%)`),
       CIw = sim_f1(`CI width (%)`), Rok = as.character(R_ok))] |>
  gt(groupname_col = "Scenario") |>
  tab_caption(paste0("Simulation 2: bias, RMSE, coverage and CI width with zeros allowed ($\\gamma$ ",
                     "estimated). Analytical CIs. $R = 2{,}000$ replications.")) |>
  cols_label(Bias = "Bias (%)", BiasBC = "Bias BC (%)", Biasmed = "Bias med (%)", RMSE = "RMSE",
             Cov = "Coverage (%)", CIw = "CI width (%)", Rok = "$R_{\\text{ok}}$") |>
  cols_align("right", columns = c(Bias, BiasBC, Biasmed, RMSE, Cov, CIw, Rok)) |>
  tab_source_note(paste0("Note: bias and CI width as a percentage of true $\\xi$; RMSE on the original ",
                         "scale. NLS bias, bias BC, bias med and RMSE use median-based summaries ",
                         "(root-median-squared error) under NB-generated data; other estimators use ",
                         "mean-based summaries.")) |>
  tab_options(table.font.size = px(9), source_notes.font.size = px(8)) |>
  gt_to_tex("tables/tblA-sim2-analytical.tex", label = "tbl-sim2-analytical")

## ---- Table: Simulation 2, analytical vs bootstrap CI ----
ci <- fread(file.path(sim_dir, "ci_comparison.csv"))[, Model := sim_ml[model]]
ci[, ci_ord := match(ci, c("Analytical", "Bootstrap"))]
setorder(ci, scenario, model, ci_ord)
ci[, .(Scenario = scenario, CI = ci, Model, Cov = sim_f1(`Coverage (%)`), CIw = sim_f1(`CI width (%)`))] |>
  gt(groupname_col = "Scenario") |>
  tab_caption("Simulation 2: analytical vs bootstrap CI coverage and width, by scenario and estimator.") |>
  cols_label(Cov = "Coverage (%)", CIw = "CI width (%)") |>
  cols_align("right", columns = c(Cov, CIw)) |>
  tab_source_note(paste0("Note: coverage is the share of replications whose 95% CI contains $\\xi$; ",
                         "CI width is the median width as a percentage of $\\xi$. Bootstrap CIs use ",
                         "$R_{\\text{boot}} = 499$ per replicate without clustering. NLS CI widths can be ",
                         "very large under NB-generated scenarios.")) |>
  tab_options(table.font.size = px(10), source_notes.font.size = px(9)) |>
  gt_to_tex("tables/tblA-sim2-ci.tex", label = "tbl-sim2-ci-comparison")

## ---- Table: Simulation 2, gamma recovery ----
g <- fread(file.path(sim_dir, "sim2_gamma_summary.csv"))[, Model := sim_ml[model]]
setorder(g, scenario, model)
g[, .(Scenario = scenario, Model, Median = sim_f4(gamma_median), Mean = sim_f4(gamma_mean),
      Q75 = sim_f4(gamma_q75), boundary = sim_f1(pct_boundary), gt01 = sim_f1(pct_gt_01))] |>
  gt(groupname_col = "Scenario") |>
  tab_caption(paste0("Simulation 2: recovery of $\\gamma$ across scenarios and estimators. Median, mean, ",
                     "and the share of replications where $\\hat\\gamma$ hits the upper boundary ",
                     "($\\geq 0.49$). True $\\gamma = 0.005$.")) |>
  cols_label(boundary = "% boundary", gt01 = "% $>0.1$") |>
  cols_align("right", columns = c(Median, Mean, Q75, boundary, gt01)) |>
  tab_source_note(paste0("Note: boundary = $\\hat\\gamma \\geq 0.49$ (upper optimisation bound). The ",
                         "unconstrained NB hits this boundary in 16% of Po(n)-Po(m) replications, ",
                         "inflating $\\hat\\gamma$ by a factor of 24 relative to the true value.")) |>
  tab_options(table.font.size = px(10), source_notes.font.size = px(9)) |>
  gt_to_tex("tables/tblA-sim2-gamma.tex", label = "tbl-sim2-gamma")

