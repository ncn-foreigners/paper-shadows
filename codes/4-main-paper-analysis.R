## Main-paper analysis: modelling plus every main-text figure and table.
## Uses full_database_processed (2-prepare-data.R) and estimate_hidden_pop()/popsize()
## from the uncounted package, and the helpers in 1-functions.R.
##
## Part A -- Modelling (specifications, fits, bootstraps)
## Part B -- Main-text figures and tables

## ZUS-based modelling dataset (reference population N = social-insurance register)
model_zus <- full_database_processed[pop_insured > 0,
  .(year = as.factor(year), sex, country_code,
    m = border, n = police, N = pop_insured,
    ukr, simplified_proc, continent)]


## ---- Estimates for Figure 5 (fig5-init-estimates.pdf; plotted in Part B below) ----
## Separate yearly estimation of the unauthorised population total by NB-MLE without offset,
## Poisson PMLE with offset, and NB-MLE with offset. Two CIs (95% and 80%) reported per spec.
## Result (baseline_cmp) is kept in memory and saved to results/baseline-cmp.rds.
years_vec <- sort(unique(as.character(model_zus$year)))

## --- 1. Separate yearly NB (Zhang 2008, aggregated, no gamma) ---
baseline_yearly <- rbindlist(lapply(years_vec, function(y) {
  d_y <- model_zus[as.character(year) == y]
  d_y[, flag := (n > 0 & m > 0 & n / N < 1)]
  d_y[, country_agg := fifelse(flag, as.character(country_code),
                                paste0("_rest_", sex))]
  rest_n <- d_y[startsWith(country_agg, "_rest_"), .(n_rest = sum(n)), by = sex]
  for (s in rest_n[n_rest == 0, sex]) {
    cands <- d_y[sex == s & flag & n > 0]
    if (nrow(cands) > 0) {
      donor <- cands[which.min(N), country_agg]
      d_y[sex == s & country_agg == donor,
          country_agg := paste0("_rest_", s)]
    }
  }
  d_agg <- d_y[, .(m = sum(m), n = sum(n), N = sum(N), sex = first(sex)),
               by = country_agg]
  setnames(d_agg, "country_agg", "country_code")
  d_agg <- d_agg[n > 0 & m > 0]
  n_kept <- sum(!startsWith(d_agg$country_code, "_rest_"))
  fit_nb <- tryCatch(
    suppressWarnings(estimate_hidden_pop(
      data = d_agg, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
      method = "nb", gamma = NULL, vcov = "HC1",
      countries = ~ country_code)),
    error = function(e) NULL)
  if (is.null(fit_nb)) return(NULL)

  ps <- popsize(fit_nb)
  ps2 <- popsize(fit_nb, level = 0.80)

  data.table(year = as.integer(y),
             estimate = sum(ps$estimate),
             lower = sum(ps$lower), upper = sum(ps$upper),
             lower_80 = sum(ps2$lower), upper_80 = sum(ps2$upper),
             n_countries = n_kept,
             Spec = "NB-MLE without offset")
}))

## --- 2. Pooled Poisson with year effects, gamma estimated ---
fit_pool_po <- estimate_hidden_pop(
  data = model_zus, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
  method = "poisson",
  cov_alpha = ~ year, cov_beta = ~ year,
  gamma = "estimate", vcov = "HC1",
  countries = ~ country_code)
ps_po <- as.data.table(popsize(fit_pool_po, by = ~ year))
ps_po2 <- as.data.table(popsize(fit_pool_po, by = ~ year, level = 0.80))

ps_po[, `:=`(year = as.integer(as.character(group)),
             Spec = "Poisson PMLE with offset",
             lower_80 = ps_po2$lower,
             upper_80 = ps_po2$upper)]


## --- 3. Pooled NB with year effects, gamma estimated ---
fit_pool_nb <- estimate_hidden_pop(
  data = model_zus, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
  method = "nb",
  cov_alpha = ~ year, cov_beta = ~ year,
  gamma = "estimate", vcov = "HC1",
  countries = ~ country_code)

ps_nb <- as.data.table(popsize(fit_pool_nb, by = ~ year))
ps_nb2 <- as.data.table(popsize(fit_pool_nb, by = ~ year, level = 0.80))

ps_nb[, `:=`(year = as.integer(as.character(group)),
             Spec = "NB-MLE with offset",
             lower_80 = ps_nb2$lower,
             upper_80 = ps_nb2$upper)]

baseline_cmp <- rbind(
  baseline_yearly[, .(year, estimate, lower, upper, lower_80, upper_80, Spec)],
  ps_po[, .(year, estimate, lower, upper, lower_80, upper_80, Spec)],
  ps_nb[, .(year, estimate, lower, upper, lower_80, upper_80, Spec)]
)

baseline_cmp[, Spec := factor(Spec, levels = c(
  "NB-MLE without offset",
  "Poisson PMLE with offset",
  "NB-MLE with offset"))]

saveRDS(baseline_cmp, "results/baseline-cmp.rds")

# Section -- Parameter specification and estimation

## Fit specifications for model comparison
specs <- list(
 S1 = list(cov_alpha = ~ 1,                        cov_beta = ~ 1,     label = "Intercept only"),
 S2 = list(cov_alpha = ~ year,                     cov_beta = ~ year,  label = "Year"),
 S3 = list(cov_alpha = ~ year + sex,               cov_beta = ~ year,  label = "Year + sex"),
 S4 = list(cov_alpha = ~ year * ukr,               cov_beta = ~ year,  label = "Year x UKR"),
 S5 = list(cov_alpha = ~ year * simplified_proc,   cov_beta = ~ year,  label = "Year x simplified"),
 S6 = list(cov_alpha = ~ year + continent,         cov_beta = ~ year,  label = "Year + continent"),
 S7 = list(cov_alpha = ~ year + sex + ukr,         cov_beta = ~ year,  label = "Year + sex + UKR"),
 S8 = list(cov_alpha = ~ year * ukr + sex,         cov_beta = ~ year,  label = "Year x UKR + sex")
)

## Fit every specification (Poisson + NB). Expensive -> cached to results/ and reused on
## later runs (e.g. by 5-supplement.R); set recompute_fits <- TRUE to refit from scratch.
recompute_fits <- FALSE
if (!recompute_fits && file.exists("results/fit-all.rds")) {
  fit_all <- readRDS("results/fit-all.rds")
} else {
  fit_all <- lapply(specs, function(sp) {
    po <- estimate_hidden_pop(
      data = model_zus, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
      method = "poisson", cov_alpha = sp$cov_alpha, cov_beta = sp$cov_beta,
      countries = ~ country_code
    )
    nb <- estimate_hidden_pop(
      data = model_zus, observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
      method = "nb", cov_alpha = sp$cov_alpha, cov_beta = sp$cov_beta,
      countries = ~ country_code
    )
    list(po = po, nb = nb)
  })
  saveRDS(fit_all, "results/fit-all.rds")
}

## S8 is the main-paper specification; also persisted individually for reuse in the appendix
fit_po <- fit_all$S8$po
fit_nb <- fit_all$S8$nb
saveRDS(fit_po, "results/fit-po-s8.rds")
saveRDS(fit_nb, "results/fit-nb-s8.rds")

## Cluster fractional-weighted bootstrap (FWB) of population size for S8, by group.
## R = 999 replicates, cluster by country; each call is seeded (seed = 2026) so it is
## reproducible independent of run order. Results are cached under results/ and reused on
## later runs -- set recompute_boot <- TRUE (or delete the .rds files) to regenerate.
recompute_boot <- FALSE
boot_or_load <- function(file, fit, by) {
  if (!recompute_boot && file.exists(file)) return(readRDS(file))
  b <- bootstrap_popsize(fit, by = by, R = 999, cluster = ~ country_code,
                         level = 0.95, seed = 2026)
  saveRDS(b, file)
  b
}

ps_ukr_po_bootstrap <- boot_or_load("results/boot-ukr-po.rds", fit_po, ~ year + ukr)
ps_ukr_nb_bootstrap <- boot_or_load("results/boot-ukr-nb.rds", fit_nb, ~ year + ukr)
ps_sex_po_bootstrap <- boot_or_load("results/boot-sex-po.rds", fit_po, ~ year + sex)
ps_sex_nb_bootstrap <- boot_or_load("results/boot-sex-nb.rds", fit_nb, ~ year + sex)




## ============================================================================
## Part B -- Main-text figures and tables
## ============================================================================

# Table 1 Simulation results, 2000 replications. Bias relative to $\xi = \sum_{i=1}^{100} \xi_i$. RMSE, root mean squared error. Coverage of $\xi$. Confidence interval (CI) of $\xi$ by \eqref{eq-xi-confint}, CI width relative to $\xi$.
## Built from the analytical-CI metrics written by 3-simulation-study.R:
##   Setup-I  (gamma = 0, regrouped 'rest' community) -> sim1_metrics_analytical.csv
##   Setup-II (gamma > 0, m_i n_i = 0 allowed)        -> sim2_metrics_analytical.csv
## Shown: Po(n)-Po(m) and NB(n)-NB(m) scenarios; estimators poisson (Po-PMLE) and nb (NB-MLE).
sim_dir  <- "results/simulation-results"
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
writeLines(sim_tex, "tables/tbl1-sim-summary.tex")




# Figure 1: Number of unauthorised migrants aged 18+ apprehended by the Polish Border Guard, by place of apprehension, 2014--2024, aggregated by half-year periods. These data refer to individuals detained or detected by the Border Guard for illegal stay. Identification at the border refers to foreigners identified as unauthorised when leaving Poland (not on legal or illegal entry to Poland).

border_country[, .(m=sum(border)), 
               keyby=.(where=fcase(where %in% c("external air", "external sea", "internal air"), 
               "air/sea",
                                   default = where), period)] |>
  transform(where = factor(where, 
                           c("country", "Ukraine", "Belarus", "Russia", "air/sea"),
                           c("Interior", 
                             "Polish-Ukraine\nborder (on exit)", 
                             "Polish-Belarus\nborder (on exit)",
                             "Polish-Russian\nborder (on exit)", 
                             "Air or sea border"))) |>
  transform(period = stri_replace_all_fixed(period, "_", " H")) |>
  transform(period = stri_replace_all_fixed(period, "II", "2")) |>
  transform(period = stri_replace_all_fixed(period, "I", "1")) |>
  ggplot(data = _, aes(x = period, y = m, group = where, color = where)) +
  geom_line() +
  geom_point() + 
  labs(x = "Period",  
      color = "Place of identification\nby Polish Border Guard", 
      y = "Observed unauthorized population") +
  theme(axis.text = element_text(angle = 45, vjust = 1, hjust = 1),
        text = element_text(size=15),
        legend.key.height = unit(1, "cm")) +
  scale_color_brewer(type = "qual", palette = "Paired") -> p1

ggsave(plot=p1, filename = "figs/fig1-border.pdf", width = 10, height = 5)

# Figure 2 - Number of 18+ foreign nationals by register, country of citizenship (Schengen, non-Schengen) and sex. All reference populations are as of the end of these years. (Tax register panel omitted in the public replication -- restricted data.)

full_database[year %in% 2017:2024,.(pesel =sum(pop_register,na.rm=T),
                 social=sum(pop_insured, na.rm=T)),
              .(year, group = paste0(schengen, "\n(", ifelse(sex=="f", "Females", "Males"), ")"))] |>
  melt(id.vars = 1:2) |>
  transform(value = ifelse(value == 0, NA, value)) |>
  transform(variable = factor(variable,
                           c("pesel", "social"),
                           c("PESEL register", "Social Insurance register"))) |>
  ggplot(data = _, aes(x = as.factor(year), y = value/1e6, fill = group)) +
  geom_col(color = "black") + 
  facet_wrap(~variable) +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  labs(x = "Year",  y = "Population (in mln)", fill = "Foreigners") +
  theme(text = element_text(size=15),
        axis.text.x=element_text(angle=45, vjust=1,hjust=1),
        legend.key.height = unit(1, "cm")) -> p2
  
ggsave(plot=p2, filename = "figs/fig2-pop-registers.pdf", width = 10, height = 5)

## Table 2 Yearly information on the observed unauthorised, auxiliary and reference 18+ populations for non-Schengen countries between 2019 and 2024 by sex (in thousands)

## aggregate non-Schengen, 2019-2024, by year and sex; convert counts to thousands
tab2_sex <- full_database_processed[,
                          .(Border = sum(border,       na.rm = TRUE) / 1e3,   # m  (Border Guard)
                            Police = sum(police,       na.rm = TRUE) / 1e3,   # n  (auxiliary, rho)
                            PESEL  = sum(pop_register,  na.rm = TRUE) / 1e3,   # N  (PESEL register)
                            ZUS    = sum(pop_insured,   na.rm = TRUE) / 1e3),  # N  (Social Insurance)
                          keyby = .(year, sex)]

tab2_tot <- tab2_sex[, lapply(.SD, sum), by = year, .SDcols = Border:ZUS]
tab2_tot[, sex := "total"]

## explicit row order so Total sits *after* Females/Males within each year
tab2_sex[, ord := fifelse(sex == "f", 1L, 2L)]
tab2_tot[, ord := 3L]
tab2 <- rbind(tab2_sex, tab2_tot)
setorder(tab2, year, ord)

tab2[, Sex  := factor(ord, 1:3, c("Females", "Males", "Total"))]
tab2[, Year := fifelse(ord == 1L, as.character(year), "")]   # year printed once per block
setcolorder(tab2, c("Year", "Sex", "Border", "Police", "PESEL", "ZUS"))
tab2[, c("year", "sex", "ord") := NULL]


tab2_gt <- tab2 |>
  gt() |>
  tab_caption(paste0("Yearly information on the observed unauthorised, auxiliary and reference ",
                     "18+ populations for non-Schengen countries between 2019 and 2024 by sex ",
                     "(in thousands)")) |>
  tab_spanner(label = "Register", columns = c(PESEL, ZUS)) |>
  fmt_number(columns = Border:ZUS, decimals = 1) |>          # 1 decimal + thousands separators
  cols_align(align = "right", columns = everything()) |>
  ## column footnotes -> marks a, b, c (PESEL/ZUS share the same note -> single mark c)
  tab_footnote(footnote = paste0("Observed unauthorised population according to Polish Border Guard ",
                                 "inspections within the country (excluding airports); this is $m$ in our model."),
               locations = cells_column_labels(columns = Border)) |>
  tab_footnote(footnote = paste0("Auxiliary population used to model the detection rate $\\rho$; ",
                                 "this is $n$ in our model."),
               locations = cells_column_labels(columns = Police)) |>
  tab_footnote(footnote = paste0("Reference population registers; serve as $N$ in our model. ",
                                 "Data recorded as of 31 December of each study year."),
               locations = cells_column_labels(columns = c(PESEL, ZUS))) |>
  ## general (unmarked) notes
  tab_source_note(source_note = "Note: All numbers are in thousands.") |>
  tab_source_note(source_note = paste0("Data from Border Guard and Police are counts of persons ",
                                       "within the year, while reference populations are stocks as of year-end.")) |>
  ## italicise the whole Total row
  tab_style(style = cell_text(style = "italic"),
            locations = cells_body(rows = Sex == "Total")) |>
  opt_footnote_marks(marks = "letters") |>
  tab_options(table.font.size        = px(16),   # gt: pt = px * 0.75 -> 16px = 12pt, line 14pt
              footnotes.font.size     = px(12),   # smaller notes (~9pt, footnotesize)
              source_notes.font.size  = px(12))

gt_to_tex(tab2_gt, "tables/tbl2-yearly-data.tex", pos = "H", label = "tbl-yearly-data")   # gt_to_tex() defined in 1-functions.R
tab2_gt

## Figure 3: Relation between standardised positive detection count $m/N$ and reference (ZUS)
## population size $N$ on log scale (alpha diagnostic).
## Figure 4: Relation between $m/N$ and detection feature $n/N$ on log scale (beta diagnostic).
## Both use register_diag_plots() (defined in 1-functions.R); ZUS is the reference register N.
pos <- register_diag_plots("pop_insured",
                           "ZUS",
                           "figs/fig3-zus-alpha.pdf",
                           "figs/fig4-zus-beta.pdf",
                           facet = facet_wrap(~year, nrow = 2),
                           extra = scale_color_brewer(type = "qual", palette = "Set1"))

## Per year-sex OLS slopes: alpha = slope(log_mN ~ log_N) + 1 ; beta = slope(log_mN ~ log_nN)
alpha_slopes <- pos[, .(alpha_hat = round(coef(lm(log_mN ~ log_N))[2] + 1, 2)), by = .(year, sex)]
beta_slopes  <- pos[, .(beta_hat  = round(coef(lm(log_mN ~ log_nN))[2],    2)), by = .(year, sex)]

# Figure 5 Separate yearly estimation of unauthorised population total by NB-MLE without offset, Poisson PMLE with offset, and NB-MLE with offset. Two confidence intervals by \eqref{eq-xi-confint} reported
## baseline_cmp is built in Part A above (also saved to results/baseline-cmp.rds).

pd <- position_dodge(width = 0.5)   # one shared dodge so all layers line up

ggplot(baseline_cmp, aes(x = as.factor(year), 
        y = estimate, colour = Spec)) +
  geom_linerange(aes(ymin = pmax(lower, 0),    ymax = pmin(upper, 450000), linewidth = "95%"),
                 position = pd) +
  geom_linerange(aes(ymin = pmax(lower_80, 0), ymax = pmin(upper_80, 450000), linewidth = "80%"),
                 position = pd) +
  geom_point(position = pd, size = 2, shape = 21, fill = "white", stroke = 0.7) +
  scale_linewidth_manual(name = "CI", breaks = c("80%", "95%"),
                         values = c("80%" = 1.6, "95%" = 0.5)) +
  coord_cartesian(ylim = c(0, 450000)) +
  scale_y_continuous(labels = scales::label_comma(),
                    breaks = seq(0, 450000, 50000)) +
  labs(x = "Year", y = expression(hat(xi)), colour = NULL,
       caption = "Point: estimate. Intervals truncated at 450,000.") +
  scale_color_brewer(type = "qual", palette = "Set1") +
  theme(text = element_text(size = 15)) -> p5

ggsave(plot=p5, filename = "figs/fig5-init-estimates.pdf", width = 10, height = 5)


## Table 3: AIC and BIC by community anchor parameter alpha_it, given yearly beta_t and constant gamma.
## Poisson PMLE (po) and NB-MLE (nb) fits per specification S1-S8 from fit_all (Part A).

## one row per specification with the four information criteria
tab3 <- rbindlist(lapply(names(fit_all), function(s) {
  x <- fit_all[[s]]
  data.table(
    Specification = s,
    Label         = specs[[s]]$label,
    po_AIC = AIC(x$po), po_BIC = BIC(x$po),
    nb_AIC = AIC(x$nb), nb_BIC = BIC(x$nb)
  )
}))

## round to whole numbers and format with thousands separators
fmt_ic <- function(v) formatC(round(v), format = "d", big.mark = ",")

body_rows <- tab3[, paste0(
  Specification, " & ", Label, " & ",
  fmt_ic(po_AIC), " & ", fmt_ic(po_BIC), " & ",
  fmt_ic(nb_AIC), " & ", fmt_ic(nb_BIC), " \\\\")]

tab3_tex <- c(
  "\\begin{table}[H]",
  "\\centering",
  paste0("\\caption{\\label{tbl-aic-specs}AIC and BIC by community anchor parameter $\\alpha_{it}$, ",
         "given yearly $\\beta_t$ and constant $\\gamma$.}"),
  "\\begin{tabular}{llrr|rr} \\toprule",
  "& & \\multicolumn{2}{c|}{Poisson PMLE} & \\multicolumn{2}{c}{NB-MLE} \\\\",
  "Specification & Label & AIC & BIC & AIC & BIC \\\\ \\midrule",
  body_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(tab3_tex, "tables/tbl3-aic-specs.tex")
tab3

# Figure - Estimated trajectory of community anchor (\(\alpha\)) by sex and Ukrainian origin, based on Poisson PMLE with S8 specification. NB. make the figure without NB-MLE 

get_alpha_by_year_ukr_sex <- function(fit) {
  X <- fit$X_alpha
  a_coefs <- fit$alpha_coefs
  d <- fit$data
  d$alpha_hat <- as.numeric(X %*% a_coefs)
  unique(d[, .(year, ukr, sex, alpha_hat)])
}

fit_po <- fit_all$S8$po
fit_nb <- fit_all$S8$nb

alpha_po <- get_alpha_by_year_ukr_sex(fit_po)
alpha_nb <- get_alpha_by_year_ukr_sex(fit_nb)
alpha_po[, method := "Poisson"]
alpha_nb[, method := "NB"]
alpha_dt <- rbind(alpha_po, alpha_nb)
alpha_dt[, ukr_label := fifelse(ukr == 1, "Ukraine", "Non-Ukraine")]
alpha_dt[, year_num := as.integer(as.character(year))]
alpha_dt[, sex := factor(sex, c("f", "m"), c("Females", "Males"))]

alpha_dt |> 
  subset(method == "Poisson") |>
  ggplot(data=_, aes(x = year_num, y = alpha_hat,
                     colour = ukr_label, linetype = sex,
                     group = interaction(ukr_label, sex))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_linetype_manual(values = c(Females = "dotted", Males = "solid"), 
                        name = "Sex") +
  labs(x = "Year", y = expression(hat(alpha)),
       colour = "Origin") +
  scale_color_brewer(type = "qual", palette = "Set1") +
  scale_y_continuous(limits = c(0.6,1.0)) +
  theme(text = element_text(size=15)) -> p6

ggsave(plot=p6, filename = "figs/fig6-alpha-trajectories.pdf", width = 10, height = 5)


# Figure - Estimated trajectory of community anchor (alpha) by sex and Ukrainian origin,
# Poisson PMLE with S8, with 80% and 95% bootstrap uncertainty bands.
# alpha_hat = X_alpha %*% alpha_coefs is linear in the coefficients, so we reuse the
# per-replicate alpha coefficients ($boot_params) from the cluster (by country) FWB already
# computed for Figure 7 / Tables 4-5 (ps_ukr_po_bootstrap; the draws are independent of the
# `by` grouping), form alpha_hat for each (year, ukr, sex) in every replicate; the central
# line is the bootstrap median and the bands are the 80%/95% percentiles. Facet by sex
# (Females dashed, Males solid); colour by origin.
an_a   <- names(fit_po$alpha_coefs)
keep_a <- !duplicated(as.data.table(fit_po$data)[, .(year, ukr, sex)])      # unique design rows
Xu_a   <- fit_po$X_alpha[keep_a, , drop = FALSE]
Ab     <- ps_ukr_po_bootstrap$boot_params[, an_a, drop = FALSE]             # R x p_alpha
Ab     <- Ab[is.finite(rowSums(Ab)), , drop = FALSE]                        # drop non-converged
qa     <- t(apply(Xu_a %*% t(Ab), 1, quantile, c(0.025, 0.10, 0.50, 0.90, 0.975), na.rm = TRUE))

alpha_bnd <- as.data.table(fit_po$data)[keep_a, .(year, ukr, sex)]
alpha_bnd[, `:=`(alpha_hat = qa[, 3],                                       # bootstrap median
                 lo95 = qa[, 1], lo80 = qa[, 2], hi80 = qa[, 4], hi95 = qa[, 5],
                 year_num  = as.integer(as.character(year)),
                 ukr_label = fifelse(ukr == 1, "Ukraine", "Non-Ukraine"),
                 sex       = factor(sex, c("f", "m"), c("Females", "Males")))]

ggplot(alpha_bnd, aes(year_num, alpha_hat, colour = ukr_label, fill = ukr_label, linetype = sex)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95, alpha = "95%"), colour = NA) +
  geom_ribbon(aes(ymin = lo80, ymax = hi80, alpha = "80%"), colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ sex) +
  scale_alpha_manual(name = "CI", breaks = c("80%", "95%"), values = c("80%" = 0.45, "95%" = 0.2)) +
  scale_linetype_manual(values = c(Females = "solid", Males = "solid"), guide = "none") +
  scale_color_brewer(name = "Origin", type = "qual", palette = "Set1") +
  scale_fill_brewer(name = "Origin", type = "qual", palette = "Set1") +
  guides(alpha    = guide_legend(override.aes = list(fill = "#E41A1C"))) + # 
  labs(x = "Year", y = expression(hat(alpha))) +
  scale_y_continuous(limits = c(0.5,1.0)) +
  theme(text = element_text(size = 15)) -> p6_bounds

ggsave(plot = p6_bounds, filename = "figs/fig6-alpha-trajectories-bounds.pdf", width = 10, height = 5)


# Figure 7 Estimated unauthorised population by origin (Ukrainian vs non-Ukrainian) and estimation method. Cluster bootstrap (by country) median with 95\% percentile CI.

## ps_ukr_{po,nb}_bootstrap (year + ukr, cluster by country) come from Part A.
## $popsize holds the median estimate and 95% percentile CI; $t holds the raw R x group
## replicate matrix (column order matches $popsize$group), from which we recover the 80% CI.
boot_to_dt <- function(boot, method) {
  ps <- as.data.table(boot$popsize)                 # group, estimate(=median), lower/upper (95%)
  q80 <- apply(boot$t, 2, quantile, probs = c(0.10, 0.90), na.rm = TRUE)
  ps[, `:=`(lower_80 = q80[1, ], upper_80 = q80[2, ], method = method)]
  ps[, c("year", "ukr") := tstrsplit(group, ", ", fixed = TRUE)]   # "year, ukr" interaction label
  ps[]
}

ukr_boot <- rbind(
  boot_to_dt(ps_ukr_po_bootstrap, "Poisson"),
  boot_to_dt(ps_ukr_nb_bootstrap, "NB")
)
ukr_boot[, ukr_label := fifelse(ukr == "1", "Ukraine", "Non-Ukraine")]
ukr_boot[, year_num  := as.integer(year)]
ukr_boot[, method    := factor(method, c("Poisson", "NB"))]

pd7 <- position_dodge(width = 0.5)   # shared dodge so both models line up

ggplot(ukr_boot, aes(x = as.factor(year_num), y = estimate, colour = method)) +
  geom_linerange(aes(ymin = pmax(lower, 0),    ymax = upper,    linewidth = "95%"),
                 position = pd7) +
  geom_linerange(aes(ymin = pmax(lower_80, 0), ymax = upper_80, linewidth = "80%"),
                 position = pd7) +
  geom_point(position = pd7, size = 2, shape = 21, fill = "white", stroke = 0.7) +
  facet_wrap(~ ukr_label, scales = "free_y") +
  scale_linewidth_manual(name = "CI", breaks = c("80%", "95%"),
                         values = c("80%" = 1.6, "95%" = 0.5)) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "Year", y = expression(hat(xi)), colour = NULL) +
  scale_color_brewer(type = "qual", palette = "Set1") +
  theme(text = element_text(size = 15)) -> p7

ggsave(plot = p7, filename = "figs/fig7-ukr-estimates.pdf", width = 10, height = 5)


## Table: estimated size of unauthorised 18+ non-Schengen population by year and sex.
## Median + 80% and 95% cluster-FWB percentile CIs for Poisson PMLE and NB-MLE (S8).
## ps_sex_{po,nb}_bootstrap (year + sex, cluster by country) come from Part A.

## observed Border (m) and ZUS reference population (N) by year x sex, plus yearly totals
obs_sex <- full_database_processed[, .(Border = sum(border), Ref = sum(pop_insured)),
                                   keyby = .(year, sex)]
obs_xi <- rbind(obs_sex,
                obs_sex[, .(sex = "total", Border = sum(Border), Ref = sum(Ref)), by = year])
obs_xi[, year := as.integer(year)]

## tidy one uncounted_boot (grouped by year + sex) into median + 80%/95% CI per group,
## then add yearly totals by summing the replicate columns *within each replicate* (boot$t)
## before taking quantiles -- CI bounds cannot simply be added across sexes.
boot_sex_dt <- function(boot, m) {
  g  <- tstrsplit(boot$popsize$group, ", ", fixed = TRUE)   # "year, sex"
  yr <- as.integer(g[[1]])
  tt <- boot$t                                              # R x group replicate matrix
  q80 <- apply(tt, 2, quantile, probs = c(0.10, 0.90), na.rm = TRUE)
  grp <- data.table(year = yr, sex = g[[2]], median = boot$popsize$estimate,
                    lo80 = q80[1, ], hi80 = q80[2, ],
                    lo95 = boot$popsize$lower, hi95 = boot$popsize$upper)
  tot <- rbindlist(lapply(unique(yr), function(y) {
    s <- rowSums(tt[, which(yr == y), drop = FALSE]); s <- s[is.finite(s)]
    data.table(year = y, sex = "total", median = median(s),
               lo80 = quantile(s, 0.10),  hi80 = quantile(s, 0.90),
               lo95 = quantile(s, 0.025), hi95 = quantile(s, 0.975))
  }))
  out <- rbind(grp, tot)
  setnames(out, c("median", "lo80", "hi80", "lo95", "hi95"),
           paste0(m, "_", c("median", "lo80", "hi80", "lo95", "hi95")))
  out[]
}

xi <- Reduce(function(a, b) merge(a, b, by = c("year", "sex")),
             list(obs_xi,
                  boot_sex_dt(ps_sex_po_bootstrap, "po"),
                  boot_sex_dt(ps_sex_nb_bootstrap, "nb")))
xi[, ord := fcase(sex == "f", 1L, sex == "m", 2L, default = 3L)]
setorder(xi, year, ord)
xi[, Sex := factor(ord, 1:3, c("Female", "Male", "Total"))]

## formatters (all figures expressed in thousands):
## estimates/reference pop -> whole thousands (10,000 -> "10"); Border -> 1 decimal (1,500 -> "1.5")
r_xi  <- function(x) formatC(round(x / 1e3), format = "d", big.mark = ",")
b_xi  <- function(x) formatC(x / 1e3, format = "f", digits = 1, big.mark = ",")
ci_xi <- function(l, h) paste0("(", r_xi(l), ", ", r_xi(h), ")")

## italicise Total rows to match Table 2: wrap each cell in {\itshape ...} (a brace group
## cannot span the & column separators, so the font switch is applied cell-by-cell)
cell <- function(v, ital) fifelse(ital, paste0("{\\itshape ", v, "}"), v)

## one body line per row; Year printed once per block, \addlinespace between years
xi[, body := paste0(
  cell(fifelse(ord == 1L, as.character(year), ""), ord == 3L), " & ",
  cell(as.character(Sex), ord == 3L), " & ",
  cell(b_xi(Border), ord == 3L), " & ", cell(r_xi(Ref), ord == 3L), " & ",
  cell(r_xi(po_median), ord == 3L), " & ",
  cell(ci_xi(po_lo80, po_hi80), ord == 3L), " & ", cell(ci_xi(po_lo95, po_hi95), ord == 3L), " & ",
  cell(r_xi(nb_median), ord == 3L), " & ",
  cell(ci_xi(nb_lo80, nb_hi80), ord == 3L), " & ", cell(ci_xi(nb_lo95, nb_hi95), ord == 3L), " \\\\",
  fifelse(ord == 3L & year != max(year), "\n\\addlinespace", ""))]

xi_tex <- c(
  "\\begin{table}[H]",
  "\\centering",
  paste0("\\caption{\\label{tbl-xi-estimates}Estimated size of unauthorised 18+ non-Schengen ",
         "population by year and sex. Observed counts (Border). Reference population (Ref. pop.) ",
         "by ZUS register. Model with S8 specification. Cluster FWB median with 80\\% and 95\\% ",
         "percentile CIs for Poisson PMLE or NB-MLE. All figures in thousands.}"),
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{llrr rcc rcc} \\toprule",
  paste0("& & & & \\multicolumn{3}{c}{Poisson PMLE} & \\multicolumn{3}{c}{NB-MLE} \\\\ ",
         "\\cline{5-7} \\cline{8-10}"),
  paste0("Year & Sex & Border & Ref. pop. & Median & 80\\% CI & 95\\% CI & ",
         "Median & 80\\% CI & 95\\% CI \\\\ \\midrule"),
  xi$body,
  "\\bottomrule",
  "\\end{tabular}",
  "}",
  "\\end{table}"
)

writeLines(xi_tex, "tables/tbl4-xi-estimates.tex")


## Table: estimated unauthorised population (Poisson PMLE cluster-FWB median) compared with
## administrative refused-residence (Refusal) and return-obligation (Return) counts, by year and
## Ukrainian origin. Estimate = cluster-FWB median from ps_ukr_po_bootstrap (year + ukr; Part A)
## -- the same object and column plotted as the Poisson series in Figure 7, so the two agree.
## Refusal and Return are external administrative figures entered as given.

est_ukr <- as.data.table(ps_ukr_po_bootstrap$popsize)
est_ukr[, c("year", "ukr") := tstrsplit(group, ", ", fixed = TRUE)]
est_ukr <- est_ukr[, .(year = as.integer(year), ukr = as.integer(ukr), Estimate = estimate)]

## administrative counts (ukr = 1 Ukrainian, ukr = 0 non-Ukrainian)
admin_ukr <- data.table(
  year    = rep(2019:2024, times = 2),
  ukr     = rep(c(1L, 0L), each = 6),
  Refusal = c(19700, 24500, 25400, 15100,  2300,  1900,    # Ukrainian
              13200, 14500, 12300, 18800, 21600, 27500),   # non-Ukrainian
  Return  = c(21700,  8700,  5400,  1200,   500,   600,
               7400,  3300,  4700,  7200, 10000, 11600))

cmp_ukr <- merge(admin_ukr, est_ukr, by = c("year", "ukr"))

## one wide row per year: Ukrainian block (ukr 1) then non-Ukrainian block (ukr 0)
r_h <- function(x) formatC(round(x / 100) * 100, format = "d", big.mark = ",")   # nearest hundred
blk <- function(d) d[order(year), paste0(r_h(Estimate), " & ", r_h(Refusal), " & ", r_h(Return))]

body_ukr <- paste0(2019:2024, " & ",
                   blk(cmp_ukr[ukr == 1]), " & ", blk(cmp_ukr[ukr == 0]), " \\\\")

ukr_tex <- c(
  "\\begin{table}[H]",
  "\\centering",
  paste0("\\caption{\\label{tbl-comparison-ukr}Estimated unauthorised population size (Estimate), ",
         "number of refused temporary residence decisions (Refusal) or foreigner return obligations ",
         "(Return), by origin (Ukrainian vs non-Ukrainian).}"),
  "\\begin{tabular}{lrrr|rrr} \\toprule",
  "& \\multicolumn{3}{c|}{Ukrainian} & \\multicolumn{3}{c}{Non-Ukrainian} \\\\",
  "Year & Estimate & Refusal & Return & Estimate & Refusal & Return \\\\ \\midrule",
  body_ukr,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(ukr_tex, "tables/tbl5-comparison-ukr.tex")

##