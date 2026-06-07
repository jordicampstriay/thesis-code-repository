library(tidyverse)
library(lmtest)
library(sandwich)
library(tseries)
library(xtable)
library(vars)

select <- dplyr::select

WD <- "/Users/jordi/Downloads/University/TFG/Data/Thesis/step3_dollar_analysis"
setwd(WD)

raw <- read_csv("../data_master.csv", show_col_types = FALSE)
raw$Date <- as.Date(raw$Date)

df <- raw %>%
  filter(!is.na(DXY), !is.na(HYG), !is.na(LQD), !is.na(SHV), !is.na(EMB)) %>%
  arrange(Date)

cat("Robustness sample:", as.character(range(df$Date)), "  N =", nrow(df), "\n")

lag_pairs <- list(

  c("lr_DXY", "d_log_HYG_LQD"),
  c("lr_DXY", "d_log_HYG_SHV"),
  c("lr_DXY", "lr_EMB"),
  c("d_log_HYG_LQD", "lr_DXY"),
  c("d_log_HYG_SHV", "lr_DXY"),
  c("lr_EMB", "lr_DXY"),
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
  vs <- tryCatch({
    VARselect(d2, lag.max = 15, type = "const")
  }, error = function(e) NULL)

  if (is.null(vs)) return(NULL)

  aic_lag <- vs$selection["AIC(n)"]
  bic_lag <- vs$selection["SC(n)"]

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

cat("\n=== Residual Diagnostics ===\n")

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

  bg <- bgtest(mod_u, order = 5)

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



# ── 5. Save all results ────────────────────────────────────────────────────

write_csv(lag_results, "results_lag_selection_dollar.csv")
write_csv(diag_results, "results_diagnostics_dollar.csv")

cat("\n=== Robustness diagnostics complete ===\n")
cat("Tables: table_lag_selection_dollar.tex\n")
cat("        table_diagnostics_dollar.tex\n")
