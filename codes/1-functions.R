## Shared helper functions for the pipeline (sourced by 0-run-all.R before the data stage).

## gt -> stand-alone LaTeX, forcing a float position and restoring raw math ($...$).
## label (optional): injected as \label{...} inside the caption (\caption{\label{..}..})
## so cross-references work; gt's as_latex does not emit a label on its own.
gt_to_tex <- function(x, file, pos = "ht!", label = NULL) {
  tex <- as.character(as_latex(x))
  tex <- sub("\\\\begin\\{table\\}(\\[[^]]*\\])?", paste0("\\\\begin{table}[", pos, "]"), tex)
  ## undo gt's text-escaping so embedded LaTeX math survives (e.g. row-group headers
  ## and captions written as "$\hat{\alpha}$"): \$ -> $, \textbackslash{} -> \, \{ \} -> { }
  tex <- gsub("\\\\\\$", "$", tex)
  tex <- gsub("\\\\textbackslash\\{\\}", "\\\\", tex)
  tex <- gsub("\\\\\\{", "{", tex)
  tex <- gsub("\\\\\\}", "}", tex)
  tex <- gsub("\\\\_", "_", tex)   # restore math subscripts (e.g. $k_\gamma$)
  if (!is.null(label))
    tex <- sub("\\\\caption\\{", paste0("\\\\caption{\\\\label{", label, "}"), tex)
  writeLines(tex, file)
}

## alpha/beta detection diagnostics for a given reference register N. Saves the two plots
## and (invisibly) returns the positive-cell data so callers can fit the slopes.
## color_var: column mapped to colour (e.g. "sex", "year", "continent").
## facet: a ggplot2 facet layer added to the plot (e.g. facet_grid(cols = vars(year)),
##        facet_wrap(vars(continent)), facet_grid(rows = vars(sex), cols = vars(year))).
## subset: optional unquoted filter on the plotted data, e.g. subset = continent == "Africa"
##         or subset = sex == "Males" & year %in% 2022:2024 (default keeps all rows).
## extra: extra ggplot layer(s) added to both plots, e.g. a colour scale or theme
##        (extra = scale_color_brewer(type = "qual", palette = "Paired")); pass several as a list().
## pos carries year, sex, country_code, country, continent, wbd_region7, schengen for use here.
register_diag_plots <- function(reg_col, reg_label, file_alpha, file_beta,
                                color_var = "sex", facet = facet_grid(cols = vars(year)),
                                subset = TRUE, extra = NULL) {
  keep <- substitute(subset)
  pos <- full_database_processed[, .(year = as.factor(year), sex, country_code, country,
                                     continent, wbd_region7, schengen,
                                     m = border, n = police, N = get(reg_col))][
                                   N > 0 & m > 0 & n > 0]
  pos[, `:=`(log_mN = log(m / N), log_N = log(N), log_nN = log(n / N))]
  pos[, sex := factor(sex, c("f", "m"), c("Females", "Males"))]
  pos <- pos[eval(keep, pos, parent.frame())]   # apply `subset` (NSE) on the final columns

  ## alpha diagnostic: log(m/N) vs log(N)
  ggplot(pos, aes(x = log_N, y = log_mN, label = country_code, color = .data[[color_var]])) +
    geom_text(alpha = 0.6, show.legend = TRUE, aes(size = sqrt(N))) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
    facet +
    scale_size_continuous(range = c(1.5, 4)) +
    coord_cartesian(clip = "off") +
    labs(x = expression(log(N)), y = expression(log(m / N)),
         title = bquote(.(reg_label) ~ "register diagnostic: " * alpha < 1 ~ "(negative slope expected)"),
         size = expression(sqrt(N)), color = color_var) +
    theme(text = element_text(size = 15)) +
    extra -> p_alpha
  ggsave(plot = p_alpha, filename = file_alpha, width = 10, height = 5)

  ## beta diagnostic: log(m/N) vs log(n/N)
  ggplot(pos, aes(x = log_nN, y = log_mN, label = country_code, color = .data[[color_var]])) +
    geom_text(alpha = 0.6, show.legend = TRUE, aes(size = sqrt(N))) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
    facet +
    scale_size_continuous(range = c(1.5, 5)) +
    labs(x = expression(log(n / N)), y = expression(log(m / N)),
         title = bquote(.(reg_label) ~ "register diagnostic: " * beta > 0 ~ "(positive slope expected)"),
         size = expression(sqrt(N)), color = color_var) +
    theme(text = element_text(size = 15)) +
    extra -> p_beta
  ggsave(plot = p_beta, filename = file_beta, width = 10, height = 5)

  invisible(pos)
}


## ---- dispersion and quasi-AIC for `uncounted` fits (Burnham & Anderson 2002, Sec. 2.5) ----
## c_hat(): Pearson chi-square / residual df of a Poisson PMLE fit -- the over-dispersion
##   factor (1 under a correctly specified Poisson).
## qaic(): -2 logLik / c_hat + 2 (K + 1), where K is the parameter count AIC() uses for the
##   fit (coefficients plus the offset gamma) and the +1 counts c_hat itself. c_hat is taken
##   from the most general model of the set being compared and applied to every member of
##   the set. Defined for Poisson PMLE only (NB-MLE models the dispersion itself).
c_hat <- function(fit) {
  mu <- pmax(fit$fitted.values, 1e-8)
  sum((fit$m - mu)^2 / mu) / fit$df.residual
}
qaic <- function(fit, chat) {
  ll <- logLik(fit)
  -2 * as.numeric(ll) / chat + 2 * (attr(ll, "df") + 1)
}
