###############################################################################
#  Regenerate ADF table with proper methodology:
#    - urca::ur.df() instead of tseries::adf.test()
#    - Deterministic component explicitly specified
#    - Lag order selected by AIC, reported
#    - Level tests for ALL variables (to establish integration order)
#    - Transformed-variable tests (to confirm stationarity for analysis)
#
#  Deterministic specification rationale:
#    Level tests:
#      "drift" (constant only) for spreads, log ratios, quantities, prices
#      — spreads/ratios mean-revert; prices/quantities follow random walk
#         with possible drift; no economic reason for deterministic trend
#    Transformed-variable tests:
#      "none" for first differences and log returns
#      — returns/changes have approximately zero mean, no trend
#
#  Reference: MacKinnon (1996) for critical values.
#  Lag selection: AIC over {0, ..., floor((N-1)^(1/3))}
###############################################################################

library(readr)
library(dplyr)
library(urca)
library(xtable)

OUT_DIR <- "/Users/jordi/Downloads/University/TFG/Data/Thesis/Empirical Analysis/tables"

# ── Load data ───────────────────────────────────────────────────────────────

df <- read_csv("/Users/jordi/Downloads/University/TFG/Data/Thesis/data_master.csv",
               show_col_types = FALSE)
df$Date <- as.Date(df$Date)

# Post-2018 sample
df_post <- df %>% filter(Date >= as.Date("2018-04-03"))

# Dollar sample (all ETFs available)
df_dollar <- df %>%
  filter(!is.na(DXY), !is.na(HYG), !is.na(LQD), !is.na(SHV), !is.na(EMB))

cat(sprintf("Post-2018: %d obs (%s to %s)\n",
            nrow(df_post), min(df_post$Date), max(df_post$Date)))
cat(sprintf("Dollar:    %d obs (%s to %s)\n",
            nrow(df_dollar), min(df_dollar$Date), max(df_dollar$Date)))

# ── ADF test wrapper ───────────────────────────────────────────────────────

run_adf <- function(x, type = "drift") {
  x <- na.omit(x)
  n <- length(x)
  if (n < 30) return(list(stat = NA, lags = NA, n = NA,
                           cv1 = NA, cv5 = NA, cv10 = NA))
  max_lag <- trunc((n - 1)^(1/3))
  if (max_lag < 1) max_lag <- 1

  test <- ur.df(x, type = type, lags = max_lag, selectlags = "AIC")

  # tau statistic is first element of teststat
  tau <- test@teststat[1]
  cv  <- test@cval[1, ]  # critical values for tau

  # Actual lags used (after AIC selection)
  # testreg is summary.lm; coefficients is a matrix (rows x 4)
  n_coef <- nrow(test@testreg$coefficients)
  if (type == "none")  lags_used <- n_coef - 1       # y_{t-1} + lags
  if (type == "drift") lags_used <- n_coef - 2       # intercept + y_{t-1} + lags
  if (type == "trend") lags_used <- n_coef - 3       # intercept + trend + y_{t-1} + lags

  n_used <- length(residuals(test@testreg))

  list(stat = round(tau, 2),
       lags = lags_used,
       n = n_used,
       cv1 = round(cv[1], 2),
       cv5 = round(cv[2], 2),
       cv10 = round(cv[3], 2))
}

# ── Define all tests ──────────────────────────────────────────────────────

# Panel A: Level tests (establish integration order)
level_tests <- list(
  # Post-2018 spreads & quantities
  list(label = "SOFR--EFFR",      var = "SOFR_EFFR",  data = "post", type = "drift",
       form = "Level", sample = "Post-2018"),
  list(label = "Baa spread",      var = "Baa_spread",  data = "post", type = "drift",
       form = "Level", sample = "Post-2018"),
  list(label = "Aaa spread",      var = "Aaa_spread",  data = "post", type = "drift",
       form = "Level", sample = "Post-2018"),
  list(label = "Baa--Aaa",        var = "Baa_Aaa",     data = "post", type = "drift",
       form = "Level", sample = "Post-2018"),
  list(label = "Reserves",        var = "Reserves",    data = "post", type = "drift",
       form = "Level", sample = "Post-2018"),
  list(label = "TGA",             var = "TGA",         data = "post", type = "drift",
       form = "Level", sample = "Post-2018"),
  list(label = "ON RRP",          var = "ON_RRP",      data = "post", type = "drift",
       form = "Level", sample = "Post-2018"),

  # Dollar sample levels
  list(label = "DXY",             var = "DXY",         data = "dollar", type = "drift",
       form = "Level", sample = "2008--2026"),
  list(label = "log(HYG/LQD)",    var = "log_HYG_LQD", data = "dollar", type = "drift",
       form = "Level", sample = "2008--2026"),
  list(label = "log(HYG/SHV)",    var = "log_HYG_SHV", data = "dollar", type = "drift",
       form = "Level", sample = "2008--2026"),
  list(label = "EMB",             var = "EMB",         data = "dollar", type = "drift",
       form = "Level", sample = "2008--2026")
)

# Panel B: Transformed variables (confirm stationarity for analysis)
trans_tests <- list(
  # Post-2018 first differences
  list(label = "Baa spread",      var = "d_Baa",       data = "post", type = "none",
       form = "First diff.", sample = "Post-2018"),
  list(label = "Aaa spread",      var = "d_Aaa",       data = "post", type = "none",
       form = "First diff.", sample = "Post-2018"),
  list(label = "Baa--Aaa",        var = "d_Baa_Aaa",   data = "post", type = "none",
       form = "First diff.", sample = "Post-2018"),
  list(label = "Reserves",        var = "d_Reserves",  data = "post", type = "none",
       form = "First diff.", sample = "Post-2018"),
  list(label = "TGA",             var = "d_TGA",       data = "post", type = "none",
       form = "First diff.", sample = "Post-2018"),
  list(label = "ON RRP",          var = "d_ON_RRP",    data = "post", type = "none",
       form = "First diff.", sample = "Post-2018"),

  # Dollar sample returns / diffs
  list(label = "DXY",             var = "lr_DXY",          data = "dollar", type = "none",
       form = "Log return", sample = "2008--2026"),
  list(label = "HYG/LQD",         var = "d_log_HYG_LQD",  data = "dollar", type = "none",
       form = "$\\Delta$log ratio", sample = "2008--2026"),
  list(label = "HYG/SHV",         var = "d_log_HYG_SHV",  data = "dollar", type = "none",
       form = "$\\Delta$log ratio", sample = "2008--2026"),
  list(label = "EMB",             var = "lr_EMB",          data = "dollar", type = "none",
       form = "Log return", sample = "2008--2026")
)

# ── Run all tests ──────────────────────────────────────────────────────────

run_all <- function(test_list) {
  do.call(rbind, lapply(test_list, function(t) {
    d <- if (t$data == "post") df_post else df_dollar
    res <- run_adf(d[[t$var]], type = t$type)
    data.frame(
      Variable = t$label,
      Form     = t$form,
      Sample   = t$sample,
      Lags     = res$lags,
      N        = res$n,
      ADF_tau  = res$stat,
      CV_1     = res$cv1,
      CV_5     = res$cv5,
      CV_10    = res$cv10,
      Stationary = ifelse(res$stat < res$cv5, "Yes", "No"),
      stringsAsFactors = FALSE
    )
  }))
}

panel_a <- run_all(level_tests)
panel_b <- run_all(trans_tests)

cat("\n── Panel A: Level Tests ──\n")
print(panel_a, row.names = FALSE)
cat("\n── Panel B: Transformed Variables ──\n")
print(panel_b, row.names = FALSE)

# ── Generate LaTeX table ──────────────────────────────────────────────────

# Combine with a separator marker
all_results <- rbind(panel_a, panel_b)

# Format for LaTeX
tbl <- all_results %>%
  select(Variable, Form, Sample, Lags, N, ADF_tau, CV_1, CV_5, CV_10, Stationary)

colnames(tbl) <- c("Variable", "Form", "Sample", "Lags", "$N$",
                    "ADF $\\tau$", "1\\%", "5\\%", "10\\%", "I(0)?")

xt <- xtable(tbl,
  caption = paste0("Augmented Dickey--Fuller Unit Root Tests. ",
    "Panel~A tests variables in levels to establish integration order. ",
    "Panel~B tests transformed variables entering the analysis. ",
    "Level tests use the ``drift'' specification (constant, no trend); ",
    "tests on returns and first differences use ``none'' (no deterministic components). ",
    "Lag order selected by AIC over $\\{0, \\ldots, \\lfloor(N{-}1)^{1/3}\\rfloor\\}$. ",
    "Critical values from MacKinnon (1996)."),
  label = "tab:adf_all",
  digits = c(0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 0),
  align = c("l", "l", "l", "l", "r", "r", "r", "r", "r", "r", "l")
)

# Custom print with panel headers
sink(file.path(OUT_DIR, "table_adf_all.tex"))
cat("\\begin{table}[H]\n")
cat("\\centering\n")
cat("\\scalebox{0.72}{\n")
cat("\\begin{tabular}{lllrrrrrrl}\n")
cat("  \\toprule\n")
cat("  Variable & Form & Sample & Lags & $N$ & ADF $\\tau$ & 1\\% c.v. & 5\\% c.v. & 10\\% c.v. & I(0)? \\\\ \n")
cat("  \\midrule\n")
cat("  \\multicolumn{10}{l}{\\textit{Panel A: Level tests}} \\\\ \n")
cat("  \\addlinespace\n")

for (i in 1:nrow(panel_a)) {
  r <- panel_a[i, ]
  cat(sprintf("  %s & %s & %s & %d & %d & %.2f & %.2f & %.2f & %.2f & %s \\\\ \n",
              r$Variable, r$Form, r$Sample, r$Lags, r$N,
              r$ADF_tau, r$CV_1, r$CV_5, r$CV_10, r$Stationary))
}

cat("  \\addlinespace\n")
cat("  \\multicolumn{10}{l}{\\textit{Panel B: Analysis variables (transformed)}} \\\\ \n")
cat("  \\addlinespace\n")

for (i in 1:nrow(panel_b)) {
  r <- panel_b[i, ]
  cat(sprintf("  %s & %s & %s & %d & %d & %.2f & %.2f & %.2f & %.2f & %s \\\\ \n",
              r$Variable, r$Form, r$Sample, r$Lags, r$N,
              r$ADF_tau, r$CV_1, r$CV_5, r$CV_10, r$Stationary))
}

cat("  \\bottomrule\n")
cat("\\end{tabular}\n")
cat("}\n")
cat(paste0("\\caption{Augmented Dickey--Fuller unit root tests. ",
    "Panel~A tests variables in levels to establish integration order. ",
    "Panel~B confirms stationarity of transformed variables entering the analysis. ",
    "Level tests use the ``drift'' specification (constant, no trend); ",
    "tests on returns and first differences use ``none'' (no deterministic components). ",
    "Lag order selected by AIC over $\\{0, \\ldots, \\lfloor(N{-}1)^{1/3}\\rfloor\\}$. ",
    "Critical values from MacKinnon (1996).}\n"))
cat("\\label{tab:adf_all}\n")
cat("\\end{table}\n")
sink()

cat(sprintf("\nWritten: %s/table_adf_all.tex\n", OUT_DIR))
cat("Done.\n")
