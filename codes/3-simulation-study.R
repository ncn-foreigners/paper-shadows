#!/usr/bin/env Rscript
## =============================================================
## Simulation study for unauthorized population estimation
## Run on server: Rscript run-simulation.R
## =============================================================

## -- Dependencies ---------------------------------------------
required_pkgs <- c("data.table", "MASS", "future", "future.apply",
                   "progressr", "ggplot2", "fwb", "minpack.lm", "sandwich")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

## Install uncounted from GitHub if not already installed
if (!requireNamespace("uncounted", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", repos = "https://cloud.r-project.org")
  remotes::install_github("ncn-foreigners/uncounted")
}

library(data.table)
library(MASS)
library(future)
library(future.apply)
library(ggplot2)
library(uncounted)
library(progressr)
handlers(handler_txtprogressbar(style = 3))

cat("== Simulation study ==\n")
cat("uncounted version:", as.character(packageVersion("uncounted")), "\n")
cat("Workers:", max(1, parallel::detectCores() - 1), "\n")

## -- Output directory (relative to script location) -----------
script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
out_dir <- file.path(script_dir, "results", "simulation-results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## =============================================================
## DGP parameters (calibrated to Polish ZUS data)
## =============================================================
p_alpha <- 0.70
p_beta  <- 0.50
p_gamma <- 0.005
C <- 100  # countries per simulation

N_meanlog <- 5.5
N_sdlog   <- 2.75
N_min     <- 3

## n | N models
n_po_intercept <- -3.8073
n_po_slope     <- 0.9226
n_nb_intercept <- -4.0507
n_nb_slope     <- 0.9397
n_nb_phi       <- 1.887

## m | N, n overdispersion
m_nb_phi <- 1.9229

## Simulation settings
R_sim  <- 2000
R_boot <- 499

## =============================================================
## Helper: generate one dataset
## =============================================================
generate_sim_data <- function(seed, alpha, beta, gamma_true,
                              n_type = "poisson", m_type = "poisson",
                              aggregate = FALSE) {
  set.seed(seed)

  N <- pmax(round(rlnorm(C, N_meanlog, N_sdlog)), N_min)

  if (n_type == "poisson") {
    lambda_n <- exp(n_po_intercept + n_po_slope * log(N))
    n <- rpois(C, lambda_n)
  } else {
    mu_n <- exp(n_nb_intercept + n_nb_slope * log(N))
    n <- rnbinom(C, mu = mu_n, size = n_nb_phi)
  }

  mu_m <- N^alpha * (gamma_true + n / N)^beta
  if (m_type == "poisson") {
    m <- rpois(C, mu_m)
  } else {
    m <- rnbinom(C, mu = mu_m, size = m_nb_phi)
  }

  dt <- data.table(country = seq_len(C), m = m, n = n, N = N)
  xi_true <- sum(N^alpha)

  if (aggregate) {
    dt[, flag := n > 0 & m > 0 & n / N < 1]
    dt[, country_agg := ifelse(flag, as.character(country), "_rest_")]
    dt <- dt[, .(m = sum(m), n = sum(n), N = sum(N)), by = country_agg]
    setnames(dt, "country_agg", "country")
    dt[m == 0, m := 1L]
  }

  list(data = dt, xi_true = xi_true)
}

## =============================================================
## Helper: fit all models and extract results
## =============================================================
fit_and_extract <- function(dt, xi_true, use_gamma = NULL,
                            do_bootstrap = FALSE, R_boot = 99) {

  methods <- c("ols", "nls", "poisson", "nb", "iols")
  constrained_opts <- c(FALSE, TRUE)

  results <- list()

  ## For iOLS two-stage: get gamma from Poisson first
  gamma_from_poisson <- NULL

  for (meth in methods) {
    for (constr in constrained_opts) {
      if (meth %in% c("ols", "nls", "iols") && constr) next

      label <- if (constr) paste0(meth, "_c") else meth

      ## iOLS: use fixed gamma from Poisson (gamma="estimate" not supported)
      gamma_for_fit <- use_gamma
      if (meth == "iols") {
        if (!is.null(gamma_from_poisson)) {
          gamma_for_fit <- gamma_from_poisson
        } else {
          gamma_for_fit <- p_gamma  # fallback to true gamma
        }
      }

      fit <- tryCatch(
        suppressWarnings(estimate_hidden_pop(
          data = dt, observed = ~m, auxiliary = ~n, reference_pop = ~N,
          method = meth, gamma = gamma_for_fit, vcov = "HC3",
          constrained = constr, countries = ~country
        )),
        error = function(e) NULL
      )

      ## Store Poisson gamma for iOLS two-stage
      if (meth == "poisson" && !constr && !is.null(fit) && !is.null(fit$gamma)) {
        gamma_from_poisson <- fit$gamma
      }

      if (is.null(fit)) {
        results[[label]] <- data.table(
          model = label, alpha = NA_real_, gamma = NA_real_,
          xi = NA_real_, xi_bc = NA_real_, xi_boot_median = NA_real_,
          xi_low_a = NA_real_, xi_upp_a = NA_real_,
          xi_low_b = NA_real_, xi_upp_b = NA_real_,
          xi_true = xi_true
        )
        next
      }

      ## Analytical CI
      ps <- tryCatch(popsize(fit, bias_correction = TRUE), error = function(e) NULL)
      if (is.null(ps)) {
        results[[label]] <- data.table(
          model = label, alpha = NA_real_, gamma = NA_real_,
          xi = NA_real_, xi_bc = NA_real_, xi_boot_median = NA_real_,
          xi_low_a = NA_real_, xi_upp_a = NA_real_,
          xi_low_b = NA_real_, xi_upp_b = NA_real_,
          xi_true = xi_true
        )
        next
      }
      xi_est <- sum(ps$estimate)
      xi_bc  <- sum(ps$estimate_bc)
      xi_low_a <- sum(ps$lower)
      xi_upp_a <- sum(ps$upper)

      ## Bootstrap CI + median
      xi_boot_median <- NA_real_
      xi_low_b <- NA_real_
      xi_upp_b <- NA_real_
      gamma_est <- fit$gamma

      if (do_bootstrap) {
        boot_res <- tryCatch({
          b <- bootstrap_popsize(fit, R = R_boot, cluster = NULL,
                                 seed = NULL, verbose = FALSE)
          list(median = sum(b$popsize$estimate),
               lower  = sum(b$popsize$lower),
               upper  = sum(b$popsize$upper))
        }, error = function(e) NULL)

        if (!is.null(boot_res)) {
          xi_boot_median <- boot_res$median
          xi_low_b <- boot_res$lower
          xi_upp_b <- boot_res$upper
        }
      }

      results[[label]] <- data.table(
        model = label,
        alpha = if (constr) fit$alpha_values[1] else fit$alpha_coefs[1],
        gamma = if (!is.null(gamma_est)) gamma_est else NA_real_,
        xi = xi_est, xi_bc = xi_bc, xi_boot_median = xi_boot_median,
        xi_low_a = xi_low_a, xi_upp_a = xi_upp_a,
        xi_low_b = xi_low_b, xi_upp_b = xi_upp_b,
        xi_true = xi_true
      )
    }
  }

  rbindlist(results)
}

## =============================================================
## Helper: run one replicate
## =============================================================
run_one_replicate <- function(seed, alpha, beta, gamma_true,
                              aggregate = FALSE, do_bootstrap = FALSE,
                              R_boot = 99) {
  scenarios <- list(
    `Po(n)-Po(m)` = c("poisson", "poisson"),
    `Po(n)-NB(m)` = c("poisson", "negbin"),
    `NB(n)-Po(m)` = c("negbin", "poisson"),
    `NB(n)-NB(m)` = c("negbin", "negbin")
  )

  use_gamma <- if (aggregate) NULL else "estimate"

  results <- lapply(names(scenarios), function(sc_name) {
    sc <- scenarios[[sc_name]]
    sim <- generate_sim_data(seed, alpha, beta, gamma_true,
                             n_type = sc[1], m_type = sc[2],
                             aggregate = aggregate)

    res <- fit_and_extract(sim$data, sim$xi_true, use_gamma = use_gamma,
                           do_bootstrap = do_bootstrap, R_boot = R_boot)
    res[, scenario := sc_name]
    res
  })

  rbindlist(results)
}

## =============================================================
## Helper: compute metrics
## =============================================================
compute_metrics <- function(results_df, ci_type = "analytical") {
  if (ci_type == "analytical") {
    lo <- "xi_low_a"; up <- "xi_upp_a"
  } else {
    lo <- "xi_low_b"; up <- "xi_upp_b"
  }
  results_df[, .(
    `Bias (%)` = mean((xi - xi_true) / xi_true, na.rm = TRUE) * 100,
    `Bias BC (%)` = mean((xi_bc - xi_true) / xi_true, na.rm = TRUE) * 100,
    `Bias median (%)` = mean((xi_boot_median - xi_true) / xi_true, na.rm = TRUE) * 100,
    RMSE = sqrt(mean((xi - xi_true)^2, na.rm = TRUE)),
    `Coverage (%)` = mean(xi_true > get(lo) & xi_true < get(up), na.rm = TRUE) * 100,
    `CI width (%)` = median((get(up) - get(lo)) / xi_true, na.rm = TRUE) * 100,
    R_ok = sum(!is.na(xi))
  ), keyby = .(scenario, model)]
}

## =============================================================
## Simulation 1: with aggregation (no zeros, no gamma)
## Analytical CI only
## =============================================================
cat("\n========================================\n")
cat("Simulation 1: with aggregation (no gamma)\n")
cat("R_sim =", R_sim, "\n")
cat("========================================\n")

plan(multisession, workers = max(1, parallel::detectCores() - 1))

t1 <- system.time(with_progress({
  p <- progressor(R_sim)
  sim1_results <- future_lapply(seq_len(R_sim), function(i) {
    res <- run_one_replicate(seed = 2025 + i, alpha = p_alpha, beta = p_beta,
                             gamma_true = p_gamma, aggregate = TRUE,
                             do_bootstrap = FALSE)
    p(sprintf("Sim1: %d/%d", i, R_sim))
    res
  }, future.seed = 2025)
}))

plan(sequential)
sim1_df <- rbindlist(sim1_results, idcol = "sim")
cat("Simulation 1 done in", round(t1[3] / 60, 1), "minutes\n")

## Save raw results
fwrite(sim1_df, file.path(out_dir, "sim1_results.csv"))

## Metrics
sim1_metrics <- compute_metrics(sim1_df, ci_type = "analytical")
cat("\n--- Simulation 1: Analytical CI metrics ---\n")
print(sim1_metrics)
fwrite(sim1_metrics, file.path(out_dir, "sim1_metrics_analytical.csv"))

## =============================================================
## Simulation 2: with zeros (gamma estimated)
## Both analytical and bootstrap CI
## =============================================================
cat("\n========================================\n")
cat("Simulation 2: gamma estimated, bootstrap R =", R_boot, "\n")
cat("R_sim =", R_sim, "\n")
cat("========================================\n")

plan(multisession, workers = max(1, parallel::detectCores() - 1))

t2 <- system.time(with_progress({
  p <- progressor(R_sim)
  sim2_results <- future_lapply(seq_len(R_sim), function(i) {
    res <- run_one_replicate(seed = 2025 + i, alpha = p_alpha, beta = p_beta,
                             gamma_true = p_gamma, aggregate = FALSE,
                             do_bootstrap = TRUE, R_boot = R_boot)
    p(sprintf("Sim2: %d/%d", i, R_sim))
    res
  }, future.seed = 2025)
}))

plan(sequential)
sim2_df <- rbindlist(sim2_results, idcol = "sim")
cat("Simulation 2 done in", round(t2[3] / 60, 1), "minutes\n")

## Save raw results
fwrite(sim2_df, file.path(out_dir, "sim2_results.csv"))

## Metrics — analytical
sim2_metrics_a <- compute_metrics(sim2_df, ci_type = "analytical")
cat("\n--- Simulation 2: Analytical CI metrics ---\n")
print(sim2_metrics_a)
fwrite(sim2_metrics_a, file.path(out_dir, "sim2_metrics_analytical.csv"))

## Metrics — bootstrap
sim2_metrics_b <- compute_metrics(sim2_df, ci_type = "bootstrap")
cat("\n--- Simulation 2: Bootstrap CI metrics ---\n")
print(sim2_metrics_b)
fwrite(sim2_metrics_b, file.path(out_dir, "sim2_metrics_bootstrap.csv"))

## =============================================================
## Plots
## =============================================================
cat("\n== Generating plots ==\n")

## Bias boxplots — Sim 1
p1 <- ggplot(sim1_df[!is.na(xi)],
             aes(x = model, y = (xi - xi_true) / xi_true * 100)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~scenario, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  coord_cartesian(ylim = c(-100, 100)) +
  labs(title = "Simulation 1: Relative bias (aggregated, no gamma)",
       y = "Relative bias (%)", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "sim1_bias.pdf"), p1, width = 10, height = 6)

## Bias boxplots — Sim 2
p2 <- ggplot(sim2_df[!is.na(xi)],
             aes(x = model, y = (xi - xi_true) / xi_true * 100)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~scenario, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  coord_cartesian(ylim = c(-100, 100)) +
  labs(title = "Simulation 2: Relative bias (gamma estimated)",
       y = "Relative bias (%)", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "sim2_bias.pdf"), p2, width = 10, height = 6)

## Alpha estimation — Sim 2
p3 <- ggplot(sim2_df[!is.na(alpha)],
             aes(x = model, y = alpha)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~scenario) +
  geom_hline(yintercept = p_alpha, linetype = 2, color = "red") +
  coord_cartesian(ylim = c(0, 1.5)) +
  labs(title = "Simulation 2: Estimated alpha",
       subtitle = paste("True alpha =", p_alpha),
       y = expression(hat(alpha)), x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "sim2_alpha.pdf"), p3, width = 10, height = 6)

## Gamma recovery — Sim 2
sim2_gamma <- sim2_df[!is.na(gamma) & is.finite(gamma) & gamma < 1]
p4 <- ggplot(sim2_gamma, aes(x = model, y = gamma)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~scenario) +
  geom_hline(yintercept = p_gamma, linetype = 2, color = "red") +
  labs(title = "Simulation 2: Estimated gamma",
       subtitle = paste("True gamma =", p_gamma),
       y = expression(hat(gamma)), x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "sim2_gamma.pdf"), p4, width = 10, height = 6)

## Coverage comparison
ci_comp <- rbind(
  cbind(ci = "Analytical", sim2_metrics_a[, .(scenario, model, `Coverage (%)`, `CI width (%)`)]),
  cbind(ci = "Bootstrap",  sim2_metrics_b[, .(scenario, model, `Coverage (%)`, `CI width (%)`)])
)
fwrite(ci_comp, file.path(out_dir, "ci_comparison.csv"))

p5 <- ggplot(ci_comp, aes(x = model, y = `Coverage (%)`, fill = ci)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  facet_wrap(~scenario) +
  geom_hline(yintercept = 95, linetype = 2, color = "red") +
  labs(title = "Coverage: analytical vs bootstrap CI",
       y = "Coverage (%)", x = "", fill = "CI type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "coverage_comparison.pdf"), p5, width = 10, height = 5)

## Population size ratio: xi_hat / xi_true — Sim 1
p6 <- ggplot(sim1_df[!is.na(xi)],
             aes(x = model, y = xi / xi_true)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~scenario) +
  geom_hline(yintercept = 1, linetype = 2, color = "red") +
  coord_cartesian(ylim = c(0, 3)) +
  labs(title = "Simulation 1: Population size ratio (aggregated, no gamma)",
       y = expression(hat(xi) / xi), x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "sim1_popsize_ratio.pdf"), p6, width = 10, height = 6)

## Population size ratio: xi_hat / xi_true — Sim 2
p7 <- ggplot(sim2_df[!is.na(xi)],
             aes(x = model, y = xi / xi_true)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~scenario) +
  geom_hline(yintercept = 1, linetype = 2, color = "red") +
  coord_cartesian(ylim = c(0, 3)) +
  labs(title = "Simulation 2: Population size ratio (gamma estimated)",
       y = expression(hat(xi) / xi), x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "sim2_popsize_ratio.pdf"), p7, width = 10, height = 6)

## Point estimator comparison: plug-in vs BC vs bootstrap median — Sim 2
sim2_est_long <- melt(
  sim2_df[!is.na(xi) & !is.na(xi_boot_median),
          .(scenario, model, sim,
            `Plug-in` = xi / xi_true,
            `Bias-corrected` = xi_bc / xi_true,
            `Boot median` = xi_boot_median / xi_true)],
  id.vars = c("scenario", "model", "sim"),
  variable.name = "Estimator", value.name = "ratio"
)

p8 <- ggplot(sim2_est_long,
             aes(x = model, y = ratio, fill = Estimator)) +
  geom_boxplot(outlier.size = 0.3, position = position_dodge(0.8)) +
  facet_wrap(~scenario) +
  geom_hline(yintercept = 1, linetype = 2, color = "red") +
  coord_cartesian(ylim = c(0, 3)) +
  labs(title = "Simulation 2: Point estimator comparison",
       subtitle = "Plug-in vs bias-corrected vs bootstrap median",
       y = expression(hat(xi) / xi), x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "sim2_estimator_comparison.pdf"), p8, width = 12, height = 6)

## Bias correction comparison
bc_comp <- rbind(
  cbind(sim = "Sim1 (aggreg.)", sim1_metrics[, .(scenario, model, `Bias (%)`, `Bias BC (%)`)]),
  cbind(sim = "Sim2 (zeros)",   sim2_metrics_a[, .(scenario, model, `Bias (%)`, `Bias BC (%)`)])
)
bc_comp[, improvement_pp := round(`Bias (%)` - `Bias BC (%)`, 2)]
fwrite(bc_comp, file.path(out_dir, "bias_correction_comparison.csv"))

## =============================================================
## Done
## =============================================================
cat("\n========================================\n")
cat("All done!\n")
cat("Results saved in:", normalizePath(out_dir), "\n")
cat("Total time:", round((t1[3] + t2[3]) / 60, 1), "minutes\n")
cat("========================================\n")

sessionInfo()
