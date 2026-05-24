##############################################################################
##  Step 3 – Dollar Analysis: Robustness & Diagnostics
##  Lag selection, Breusch-Godfrey, ARCH-LM, Toda-Yamamoto
##  Jordi Camps Triay – TFG 2026
##############################################################################

library(tidyverse)
library(lmtest)
library(sandwich)
library(tseries)
library(xtable)
library(vars)

select <- dplyr::select

WD <- "/Users/jordi/Downloads/University/TFG/Data/Thesis/step3_dollar_analysis"
setwd(WD)

# ── 1. Load data (same sample as main script) ──────────────────────────────

raw <- read_csv("../data_master.csv", show_col_types = FALSE)
raw$Date <- as.Date(raw$Date)

df <- raw %>%
  filter(!is.na(DXY), !is.na(HYG), !is.na(LQD), !is.na(SHV), !is.na(EMB)) %>%
  arrange(Date)

cat("Robustness sample:", as.character(range(df$Date)), "  N =", nrow(df), "\n")

# ── 2. Information Criterion Lag Selection (AIC / BIC) ─────────────────────

# All bivariate pairs from Batteries 1-3
lag_pairs <- list(
  # Battery 1: DXY → Credit
  c("lr_DXY", "d_log_HYG_LQD"),
  c("lr_DXY", "d_log_HYG_SHV"),
  c("lr_DXY", "lr_EMB"),
  # Battery 2: Credit → DXY
  c("d_log_HYG_LQD", "lr_DXY"),
  c("d_log_HYG_SHV", "lr_DXY"),
  c("lr_EMB", "lr_DXY"),
  # Battery 3: DXY ↔ Quantities
  c("d_Reserves", "lr_DXY"),
  c("d_TGA", "lr_DXY"),
  c("lr_DXY", "d_Reserves"),
  c("lr_DXY", "d_TGA")
)

clean_name <- function(x) {
  x <- gsub("d_log_HYG_LQD", "$\\Delta\\log$(HYG/LQD)", x, fixed = TRUE)
  x <- gsub("d_log_HYG_SHV", "$\\Delta\\log$(HYG/SHV)", x, fixed = TRUE)
  x <- gsub("lr_EMB", "$\\Delta\\log$ EMB", x, fixed = TRUE)
  x <- gsub("lr_DXY", "$\\Delta\\log$ DXY", x, fixed = TRUE)
  x <- gsub("d_Reserves", "$\\Delta$Reserves", x, fixed = TRUE)
  x <- gsub("d_TGA", "$\\Delta$TGA", x, fixed = TRUE)
  x
}

cat("\n=== Lag Selection (AIC / BIC, max lag = 15) ===\n")

lag_results <- map_dfr(lag_pairs, function(pair) {
  d2 <- df[, pair] %>% drop_na()
  # Use VARselect for bivariate system
  vs <- tryCatch({
    VARselect(d2, lag.max = 15, type = "const")
  }, error = function(e) NULL)

  if (is.null(vs)) return(NULL)

  aic_lag <- vs$selection["AIC(n)"]
  bic_lag <- vs$selection["SC(n)"]   # SC = Schwarz (BIC)

  tibble(X = clean_name(pair[1]),
         Y = clean_name(pair[2]),
         AIC_lag = aic_lag,
         BIC_lag = bic_lag,
         N = nrow(d2))
})

print(lag_results)

lag_xt <- xtable(lag_results,
                 caption = "Information Criterion Lag Selection for Bivariate Granger Pairs (Dollar Sample)",
                 label = "tab:lag_selection_dollar",
                 digits = c(0, 0, 0, 0, 0, 0))
print(lag_xt, file = "table_lag_selection_dollar.tex",
      include.rownames = FALSE, sanitize.text.function = identity,
      booktabs = TRUE, floating = TRUE, table.placement = "H",
      comment = FALSE)

# ── 3. Residual Diagnostics: Breusch-Godfrey & ARCH-LM ────────────────────

cat("\n=== Residual Diagnostics ===\n")

# Key specifications from Battery 1 (forward) that showed significance
diag_specs <- list(
  list(x = "lr_DXY", y = "d_log_HYG_LQD", lags = c(1, 5)),
  list(x = "lr_DXY", y = "d_log_HYG_SHV", lags = c(1, 5)),
  list(x = "lr_DXY", y = "lr_EMB",         lags = c(1, 5)),
  list(x = "d_Reserves", y = "lr_DXY",     lags = c(1, 5)),
  list(x = "d_TGA", y = "lr_DXY",          lags = c(1, 5)),
  list(x = "lr_DXY", y = "d_TGA",          lags = c(1, 5))
)

run_diagnostics <- function(data, x_var, y_var, lag_order) {
  d2 <- data[, c(x_var, y_var)] %>% drop_na()
  n <- nrow(d2)
  maxlag <- lag_order

  Y  <- d2[(maxlag + 1):n, y_var, drop = TRUE]
  Xl <- matrix(NA, nrow = n - maxlag, ncol = lag_order)
  Yl <- matrix(NA, nrow = n - maxlag, ncol = lag_order)
  for (k in 1:lag_order) {
    Xl[, k] <- d2[(maxlag + 1 - k):(n - k), x_var, drop = TRUE]
    Yl[, k] <- d2[(maxlag + 1 - k):(n - k), y_var, drop = TRUE]
  }
  colnames(Xl) <- paste0("x_lag", 1:lag_order)
  colnames(Yl) <- paste0("y_lag", 1:lag_order)
  dat_reg <- data.frame(Y = Y, Yl, Xl)

  mod_u <- lm(Y ~ ., data = dat_reg)
  resids <- residuals(mod_u)

  # Breusch-Godfrey (order 5)
  bg <- bgtest(mod_u, order = 5)

  # ARCH-LM (order 5)
  resid2 <- resids^2
  n_r <- length(resid2)
  if (n_r > 10) {
    arch_dat <- data.frame(r2 = resid2[6:n_r])
    for (k in 1:5) arch_dat[[paste0("r2_lag", k)]] <- resid2[(6 - k):(n_r - k)]
    arch_mod <- lm(r2 ~ ., data = arch_dat)
    arch_r2 <- summary(arch_mod)$r.squared
    arch_lm <- n_r * arch_r2
    arch_p  <- pchisq(arch_lm, df = 5, lower.tail = FALSE)
  } else {
    arch_lm <- NA; arch_p <- NA
  }

  tibble(
    Direction = paste0(clean_name(x_var), " $\\rightarrow$ ", clean_name(y_var)),
    Lag = lag_order,
    BG_stat = round(bg$statistic, 3),
    BG_p    = round(bg$p.value, 4),
    BG_reject = ifelse(bg$p.value < 0.05, "Yes", "No"),
    ARCH_LM = round(arch_lm, 3),
    ARCH_p  = round(arch_p, 4),
    ARCH_reject = ifelse(!is.na(arch_p) & arch_p < 0.05, "Yes", "No")
  )
}

diag_results <- map_dfr(diag_specs, function(spec) {
  map_dfr(spec$lags, function(lag) {
    run_diagnostics(df, spec$x, spec$y, lag)
  })
})

print(diag_results)

colnames(diag_results) <- c("Direction", "Lag", "BG stat", "BG $p$",
                             "BG reject?", "ARCH-LM", "ARCH $p$", "ARCH reject?")

diag_xt <- xtable(diag_results,
                  caption = "Residual Diagnostics: Breusch-Godfrey (order 5) and ARCH-LM (order 5) Tests",
                  label = "tab:diagnostics_dollar",
                  digits = c(0, 0, 0, 3, 4, 0, 3, 4, 0))
print(diag_xt, file = "table_diagnostics_dollar.tex",
      include.rownames = FALSE, sanitize.text.function = identity,
      booktabs = TRUE, floating = TRUE, table.placement = "H",
      scalebox = 0.72, comment = FALSE)

# ── 4. Toda-Yamamoto Procedure ─────────────────────────────────────────────

cat("\n=== Toda-Yamamoto Procedure ===\n")

# All series are I(0) (log returns / first differences), so d_max = 1
# We estimate VAR(p + d_max) and test only the first p lags

run_toda_yamamoto <- function(data, x_var, y_var, p, d_max = 1) {
  d2 <- data[, c(x_var, y_var)] %>% drop_na()
  n <- nrow(d2)
  total_lags <- p + d_max

  Y  <- d2[(total_lags + 1):n, y_var, drop = TRUE]
  Xl <- matrix(NA, nrow = n - total_lags, ncol = total_lags)
  Yl <- matrix(NA, nrow = n - total_lags, ncol = total_lags)
  for (k in 1:total_lags) {
    Xl[, k] <- d2[(total_lags + 1 - k):(n - k), x_var, drop = TRUE]
    Yl[, k] <- d2[(total_lags + 1 - k):(n - k), y_var, drop = TRUE]
  }
  colnames(Xl) <- paste0("x_lag", 1:total_lags)
  colnames(Yl) <- paste0("y_lag", 1:total_lags)
  dat_reg <- data.frame(Y = Y, Yl, Xl)

  # Full model: all p + d_max lags
  mod_full <- lm(Y ~ ., data = dat_reg)

  # Restricted model: drop only first p lags of X (keep d_max extra lags)
  # We test H0: gamma_1 = ... = gamma_p = 0
  # The "extra" lags (p+1 to p+d_max) remain in both models
  x_test_cols <- paste0("x_lag", 1:p)
  x_extra_cols <- if (d_max > 0) paste0("x_lag", (p + 1):total_lags) else character(0)
  y_all_cols <- colnames(Yl)

  keep_cols <- c("Y", y_all_cols, x_extra_cols)
  mod_rest <- lm(Y ~ ., data = dat_reg[, keep_cols])

  # Wald test (chi-squared)
  wt <- waldtest(mod_rest, mod_full, test = "Chisq")
  chi2 <- wt$Chisq[2]
  chi2_p <- wt$`Pr(>Chisq)`[2]

  # Also get F-stat version
  wt_f <- waldtest(mod_rest, mod_full)
  f_stat <- wt_f$F[2]
  f_p <- wt_f$`Pr(>F)`[2]

  nobs <- nrow(dat_reg)

  tibble(Direction = paste0(clean_name(x_var), " $\\rightarrow$ ", clean_name(y_var)),
         p = p, d_max = d_max, N = nobs,
         chi2 = chi2, chi2_p = chi2_p,
         F_stat = f_stat, F_p = f_p)
}

# Run Toda-Yamamoto for key significant specifications
ty_specs <- list(
  # Battery 1 key results
  list(x = "lr_DXY", y = "lr_EMB", p = 1),
  list(x = "lr_DXY", y = "lr_EMB", p = 5),
  list(x = "lr_DXY", y = "lr_EMB", p = 10),
  list(x = "lr_DXY", y = "d_log_HYG_SHV", p = 1),
  list(x = "lr_DXY", y = "d_log_HYG_SHV", p = 5),
  list(x = "lr_DXY", y = "d_log_HYG_LQD", p = 5),
  # Battery 3 key results
  list(x = "d_Reserves", y = "lr_DXY", p = 1),
  list(x = "lr_DXY", y = "d_TGA", p = 1)
)

ty_results <- map_dfr(ty_specs, function(spec) {
  run_toda_yamamoto(df, spec$x, spec$y, spec$p)
})

print(ty_results)

sig_stars <- function(p) {
  case_when(p < 0.001 ~ "***",
            p < 0.01  ~ "**",
            p < 0.05  ~ "*",
            p < 0.10  ~ ".",
            TRUE      ~ "")
}

ty_tab <- ty_results %>%
  mutate(Sig_chi2 = sig_stars(chi2_p),
         Sig_F    = sig_stars(F_p)) %>%
  select(Direction, p, d_max, N, chi2, chi2_p, F_stat, F_p, Sig_F)

colnames(ty_tab) <- c("Direction", "$p$", "$d_{\\max}$", "$N$",
                       "$\\chi^2$", "$\\chi^2$ $p$",
                       "$F$", "$F$ $p$", "")

ty_xt <- xtable(ty_tab,
                caption = "Toda-Yamamoto Granger Non-Causality Tests ($d_{\\max} = 1$)",
                label = "tab:toda_yamamoto_dollar",
                digits = c(0, 0, 0, 0, 0, 3, 4, 3, 4, 0))
print(ty_xt, file = "table_toda_yamamoto_dollar.tex",
      include.rownames = FALSE, sanitize.text.function = identity,
      booktabs = TRUE, floating = TRUE, table.placement = "H",
      scalebox = 0.75, comment = FALSE)

# ── 5. Save all results ────────────────────────────────────────────────────

write_csv(lag_results, "results_lag_selection_dollar.csv")
write_csv(diag_results, "results_diagnostics_dollar.csv")
write_csv(ty_results, "results_toda_yamamoto_dollar.csv")

cat("\n=== Robustness diagnostics complete ===\n")
cat("Tables: table_lag_selection_dollar.tex\n")
cat("        table_diagnostics_dollar.tex\n")
cat("        table_toda_yamamoto_dollar.tex\n")
