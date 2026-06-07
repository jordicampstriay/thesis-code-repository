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

cat("Dollar analysis sample:", as.character(range(df$Date)), "  N =", nrow(df), "\n")

adf_vars <- c("lr_DXY", "d_log_HYG_LQD", "d_log_HYG_SHV", "lr_EMB",
              "d_Reserves", "d_TGA")
adf_labels <- c("$\\Delta\\log$ DXY", "$\\Delta\\log$(HYG/LQD)",
                "$\\Delta\\log$(HYG/SHV)", "$\\Delta\\log$ EMB",
                "$\\Delta$Reserves", "$\\Delta$TGA")

adf_results <- map_dfr(seq_along(adf_vars), function(i) {
  x <- na.omit(df[[adf_vars[i]]])
  tt <- adf.test(x, alternative = "stationary")
  tibble(Variable  = adf_labels[i],
         N         = length(x),
         `ADF stat` = round(tt$statistic, 2),
         `$p$-value` = round(tt$p.value, 4))
})

print(adf_results)

adf_xt <- xtable(adf_results,
                 caption = "ADF Unit Root Tests (Dollar Analysis Sample)",
                 label   = "tab:adf_dollar",
                 digits  = c(0, 0, 0, 2, 4))
print(adf_xt, file = "table_adf_dollar.tex",
      include.rownames = FALSE, sanitize.text.function = identity,
      booktabs = TRUE, floating = TRUE, table.placement = "H",
      comment = FALSE)

sum_vars <- c("lr_DXY", "d_log_HYG_LQD", "d_log_HYG_SHV", "lr_EMB")
sum_labels <- c("$\\Delta\\log$ DXY (\\%)", "$\\Delta\\log$(HYG/LQD)",
                "$\\Delta\\log$(HYG/SHV)", "$\\Delta\\log$ EMB (\\%)")

sumstats <- map_dfr(seq_along(sum_vars), function(i) {
  x <- na.omit(df[[sum_vars[i]]])
  tibble(Variable = sum_labels[i],
         N     = length(x),
         Mean  = round(mean(x), 5),
         SD    = round(sd(x), 5),
         Skew  = round(moments::skewness(x), 2),
         Kurt  = round(moments::kurtosis(x), 2),
         Min   = round(min(x), 4),
         Max   = round(max(x), 4))
})

print(sumstats)

ss_xt <- xtable(sumstats,
                caption = "Summary Statistics: Daily Log Returns / Changes (2008--2026)",
                label   = "tab:sumstats_dollar",
                digits  = c(0, 0, 0, 5, 5, 2, 2, 4, 4))
print(ss_xt, file = "table_sumstats_dollar.tex",
      include.rownames = FALSE, sanitize.text.function = identity,
      booktabs = TRUE, floating = TRUE, table.placement = "H",
      comment = FALSE)

cor_vars   <- c("lr_DXY", "d_log_HYG_LQD", "d_log_HYG_SHV", "lr_EMB")
cor_labels <- c("DXY", "HYG/LQD", "HYG/SHV", "EMB")

cor_data <- df[, cor_vars] %>% drop_na()
cor_mat  <- cor(cor_data, use = "complete.obs")
rownames(cor_mat) <- colnames(cor_mat) <- cor_labels

cor_xt <- xtable(cor_mat,
                 caption = "Contemporaneous Correlation Matrix (Daily Log Returns)",
                 label   = "tab:corr_dollar",
                 digits  = 3)
print(cor_xt, file = "table_corr_dollar.tex",
      sanitize.text.function = identity,
      booktabs = TRUE, floating = TRUE, table.placement = "H",
      comment = FALSE)

run_granger <- function(data, x_var, y_var, lag_order) {
  d2 <- data[, c(x_var, y_var)] %>% drop_na()
  n  <- nrow(d2)
  fmla_u <- paste0(y_var, " ~ ",
                   paste0("lag(", y_var, ", ", 1:lag_order, ")", collapse = " + "),
                   " + ",
                   paste0("lag(", x_var, ", ", 1:lag_order, ")", collapse = " + "))
  fmla_r <- paste0(y_var, " ~ ",
                   paste0("lag(", y_var, ", ", 1:lag_order, ")", collapse = " + "))

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
  mod_r <- lm(Y ~ ., data = dat_reg[, c("Y", colnames(Yl))])

  # Classical F-test
  wt <- waldtest(mod_r, mod_u)
  f_cl <- wt$F[2]
  p_cl <- wt$`Pr(>F)`[2]

  wt_hac <- waldtest(mod_r, mod_u, vcov = vcovHAC)
  f_hac <- wt_hac$F[2]
  p_hac <- wt_hac$`Pr(>F)`[2]

  beta1 <- coef(mod_u)[paste0("x_lag1")]

  nobs <- nrow(dat_reg)

  tibble(X = x_var, Y = y_var, Lag = lag_order, N = nobs,
         F_classical = f_cl, p_classical = p_cl,
         F_HAC = f_hac, p_HAC = p_hac,
         beta1 = beta1)
}

sig_stars <- function(p) {
  case_when(p < 0.001 ~ "***",
            p < 0.01  ~ "**",
            p < 0.05  ~ "*",
            p < 0.10  ~ ".",
            TRUE      ~ "")
}

cat("\n=== Battery 1: DXY → Credit measures ===\n")

b1_pairs <- list(
  c("lr_DXY", "d_log_HYG_LQD"),
  c("lr_DXY", "d_log_HYG_SHV"),
  c("lr_DXY", "lr_EMB")
)

b1 <- map_dfr(b1_pairs, function(pair) {
  map_dfr(c(1, 5, 10), function(lag) {
    run_granger(df, pair[1], pair[2], lag)
  })
})

print(b1 %>% select(X, Y, Lag, N, F_classical, p_classical, F_HAC, p_HAC))

cat("\n=== Battery 2: Credit → DXY (reverse) ===\n")

b2_pairs <- list(
  c("d_log_HYG_LQD", "lr_DXY"),
  c("d_log_HYG_SHV", "lr_DXY"),
  c("lr_EMB", "lr_DXY")
)

b2 <- map_dfr(b2_pairs, function(pair) {
  map_dfr(c(1, 5, 10), function(lag) {
    run_granger(df, pair[1], pair[2], lag)
  })
})

print(b2 %>% select(X, Y, Lag, N, F_classical, p_classical, F_HAC, p_HAC))

cat("\n=== Battery 3: DXY ↔ Quantities ===\n")

b3_pairs <- list(
  c("d_Reserves", "lr_DXY"),
  c("d_TGA", "lr_DXY"),
  c("lr_DXY", "d_Reserves"),
  c("lr_DXY", "d_TGA")
)

b3 <- map_dfr(b3_pairs, function(pair) {
  map_dfr(c(1, 5, 10), function(lag) {
    run_granger(df, pair[1], pair[2], lag)
  })
})

print(b3 %>% select(X, Y, Lag, N, F_classical, p_classical, F_HAC, p_HAC))

all_gc <- bind_rows(
  b1 %>% mutate(Battery = "1_DXY_to_Credit"),
  b2 %>% mutate(Battery = "2_Credit_to_DXY"),
  b3 %>% mutate(Battery = "3_DXY_Quantities")
)

write_csv(all_gc, "results_granger_dollar.csv")

clean_name <- function(x) {
  x <- gsub("d_log_HYG_LQD", "$\\Delta\\log$(HYG/LQD)", x, fixed = TRUE)
  x <- gsub("d_log_HYG_SHV", "$\\Delta\\log$(HYG/SHV)", x, fixed = TRUE)
  x <- gsub("lr_EMB", "$\\Delta\\log$ EMB", x, fixed = TRUE)
  x <- gsub("lr_DXY", "$\\Delta\\log$ DXY", x, fixed = TRUE)
  x <- gsub("d_Reserves", "$\\Delta$Reserves", x, fixed = TRUE)
  x <- gsub("d_TGA", "$\\Delta$TGA", x, fixed = TRUE)
  x
}

make_gc_table <- function(data, cap, lab, fname) {
  tab <- data %>%
    mutate(Direction = paste0(clean_name(X), " $\\rightarrow$ ", clean_name(Y)),
           Sig_cl  = sig_stars(p_classical),
           Sig_hac = sig_stars(p_HAC)) %>%
    select(Direction, Lag, N, F_classical, p_classical, Sig_cl,
           F_HAC, p_HAC, Sig_hac)

  colnames(tab) <- c("Direction", "Lag", "$N$",
                      "$F$ (class.)", "$p$ (class.)", "",
                      "$F$ (HAC)", "$p$ (HAC)", "")

  xt <- xtable(tab, caption = cap, label = lab,
               digits = c(0,0,0,0, 3,4,0, 3,4,0))
  print(xt, file = fname,
        include.rownames = FALSE, sanitize.text.function = identity,
        booktabs = TRUE, floating = TRUE, table.placement = "H",
        scalebox = 0.72, comment = FALSE)
}

make_gc_table(b1,
              "Granger Causality: DXY $\\rightarrow$ Credit/Risk Measures (2008--2026)",
              "tab:gc_dxy_credit", "table_granger_dxy_credit.tex")

make_gc_table(b2,
              "Granger Causality: Credit/Risk Measures $\\rightarrow$ DXY (Reverse, 2008--2026)",
              "tab:gc_credit_dxy", "table_granger_credit_dxy.tex")

make_gc_table(b3,
              "Granger Causality: DXY $\\leftrightarrow$ Reserves, TGA (2008--2026)",
              "tab:gc_dxy_qty", "table_granger_dxy_qty.tex")

roll_cor <- function(x, y, w = 250) {
  n <- length(x)
  rc <- rep(NA_real_, n)
  for (i in w:n) {
    idx <- (i - w + 1):i
    ok  <- complete.cases(x[idx], y[idx])
    if (sum(ok) > 30) rc[i] <- cor(x[idx][ok], y[idx][ok])
  }
  rc
}

df_rc <- df %>%
  mutate(
    rc_HYG_LQD = roll_cor(lr_DXY, d_log_HYG_LQD, 250),
    rc_HYG_SHV = roll_cor(lr_DXY, d_log_HYG_SHV, 250),
    rc_EMB     = roll_cor(lr_DXY, lr_EMB, 250)
  )

pdf("fig_rolling_corr_dxy.pdf", width = 9, height = 5)
par(mar = c(4, 4, 2.5, 1), mgp = c(2.5, 0.8, 0))
plot(df_rc$Date, df_rc$rc_HYG_LQD, type = "n",
     ylim = c(-0.7, 0.3), xlab = "", ylab = "Correlation",
     main = "250-Day Rolling Correlation: DXY Returns vs Credit Measures")
abline(h = 0, col = "grey60", lty = 2)
lines(df_rc$Date, df_rc$rc_HYG_LQD, col = "#2166AC", lwd = 1.5)
lines(df_rc$Date, df_rc$rc_HYG_SHV, col = "#D6604D", lwd = 1.5)
lines(df_rc$Date, df_rc$rc_EMB,     col = "#1B7837", lwd = 1.5)
abline(v = as.Date("2020-03-01"), col = "red", lty = 3)
text(as.Date("2020-03-01"), 0.25, "COVID", col = "red", cex = 0.8, pos = 4)
abline(v = as.Date("2008-09-15"), col = "darkred", lty = 3)
text(as.Date("2008-09-15"), 0.25, "Lehman", col = "darkred", cex = 0.8, pos = 4)
legend("bottomleft",
       legend = c("DXY vs HYG/LQD", "DXY vs HYG/SHV", "DXY vs EMB"),
       col = c("#2166AC", "#D6604D", "#1B7837"), lwd = 2, cex = 0.85,
       bg = "white")
dev.off()


roll_granger <- function(data, x_var, y_var, lag_order = 5, window = 400) {
  d2 <- data[, c("Date", x_var, y_var)] %>% drop_na()
  n <- nrow(d2)
  pvals <- rep(NA_real_, n)
  dates <- d2$Date

  for (i in window:n) {
    idx <- (i - window + 1):i
    sub <- d2[idx, ]
    tryCatch({
      res <- run_granger(sub, x_var, y_var, lag_order)
      pvals[i] <- res$p_HAC 
    }, error = function(e) {})
  }
  tibble(Date = dates, p_value = pvals)
}

cat("\nComputing rolling Granger (400-day, lag 5, HAC-robust)...\n")

rg_hyglqd <- roll_granger(df, "lr_DXY", "d_log_HYG_LQD") %>%
  mutate(Direction = "DXY → HYG/LQD")
rg_hygshv <- roll_granger(df, "lr_DXY", "d_log_HYG_SHV") %>%
  mutate(Direction = "DXY → HYG/SHV")
rg_emb <- roll_granger(df, "lr_DXY", "lr_EMB") %>%
  mutate(Direction = "DXY → EMB")

rg_all <- bind_rows(rg_hyglqd, rg_hygshv, rg_emb)
write_csv(rg_all, "results_rolling_granger_dollar.csv")

pdf("fig_rolling_granger_dxy.pdf", width = 9, height = 5)
par(mar = c(4, 4, 2.5, 1), mgp = c(2.5, 0.8, 0))
plot(rg_hyglqd$Date, rg_hyglqd$p_value, type = "n",
     ylim = c(0, 1), xlab = "", ylab = "HAC-robust p-value",
     main = "Rolling Granger p-values: DXY → Credit Measures\n(400-day window, lag 5)")
abline(h = 0.05, col = "black", lty = 2)
lines(rg_hyglqd$Date, rg_hyglqd$p_value, col = "#2166AC", lwd = 1.3)
lines(rg_hygshv$Date, rg_hygshv$p_value, col = "#D6604D", lwd = 1.3)
lines(rg_emb$Date,    rg_emb$p_value,    col = "#1B7837", lwd = 1.3)
abline(v = as.Date("2020-03-01"), col = "red", lty = 3)
text(as.Date("2020-03-01"), 0.95, "COVID", col = "red", cex = 0.8, pos = 4)
abline(v = as.Date("2008-09-15"), col = "darkred", lty = 3)
text(as.Date("2008-09-15"), 0.95, "Lehman", col = "darkred", cex = 0.8, pos = 4)
legend("topright",
       legend = c("DXY → HYG/LQD", "DXY → HYG/SHV", "DXY → EMB",
                  "5% significance"),
       col = c("#2166AC", "#D6604D", "#1B7837", "black"),
       lty = c(1, 1, 1, 2), lwd = c(2, 2, 2, 1), cex = 0.8, bg = "white")
dev.off()


pdf("fig_scatter_dxy.pdf", width = 10, height = 3.5)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1), mgp = c(2.5, 0.8, 0))

scatter_pair <- function(x, y, xlab, ylab, main_title) {
  ok <- complete.cases(x, y)
  plot(x[ok], y[ok], pch = 16, cex = 0.3, col = rgb(0.2, 0.2, 0.6, 0.3),
       xlab = xlab, ylab = ylab, main = main_title)
  abline(h = 0, v = 0, col = "grey60", lty = 2)
  fit <- lm(y[ok] ~ x[ok])
  abline(fit, col = "red", lwd = 2)
  r2 <- summary(fit)$r.squared
  legend("topright",
         legend = bquote(R^2 == .(sprintf("%.3f", r2))),
         bty = "n", cex = 0.9)
}

scatter_pair(df$lr_DXY, df$d_log_HYG_LQD,
             expression(Delta*log~DXY~("%")),
             expression(Delta*log~(HYG/LQD)),
             "DXY vs Risk Aversion")

scatter_pair(df$lr_DXY, df$d_log_HYG_SHV,
             expression(Delta*log~DXY~("%")),
             expression(Delta*log~(HYG/SHV)),
             "DXY vs Flight to Quality")

scatter_pair(df$lr_DXY, df$lr_EMB,
             expression(Delta*log~DXY~("%")),
             expression(Delta*log~EMB~("%")),
             "DXY vs EM Bonds")

dev.off()

t0_hyglqd <- df$HYG[1] / df$LQD[1]
t0_hygshv <- df$HYG_SHV[1]
t0_emb    <- df$EMB[1]

pdf("fig_ts_dollar_credit.pdf", width = 10, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.5, 0.8, 0))

plot(df$Date, df$DXY, type = "l", col = "#2166AC", lwd = 1.2,
     xlab = "", ylab = "Index level", main = "(a) US Dollar Index (DXY)")
abline(v = as.Date("2008-09-15"), col = "darkred", lty = 3)
abline(v = as.Date("2020-03-01"), col = "red", lty = 3)

hyg_lqd <- df$HYG / df$LQD
plot(df$Date, hyg_lqd, type = "l", col = "#D6604D", lwd = 1.2,
     xlab = "", ylab = "Ratio", main = "(b) HYG/LQD (Risk Aversion)")
abline(v = as.Date("2008-09-15"), col = "darkred", lty = 3)
abline(v = as.Date("2020-03-01"), col = "red", lty = 3)

plot(df$Date, df$HYG_SHV, type = "l", col = "#B2182B", lwd = 1.2,
     xlab = "", ylab = "Ratio", main = "(c) HYG/SHV (Flight to Quality)")
abline(v = as.Date("2008-09-15"), col = "darkred", lty = 3)
abline(v = as.Date("2020-03-01"), col = "red", lty = 3)

plot(df$Date, df$EMB, type = "l", col = "#1B7837", lwd = 1.2,
     xlab = "", ylab = "USD", main = "(d) EMB (EM Bond ETF)")
abline(v = as.Date("2008-09-15"), col = "darkred", lty = 3)
abline(v = as.Date("2020-03-01"), col = "red", lty = 3)

dev.off()

summ <- all_gc %>%
  group_by(Battery) %>%
  summarise(
    Tests   = n(),
    Sig_5   = sum(p_classical < 0.05),
    Sig_10  = sum(p_classical < 0.10),
    Sig_HAC_5  = sum(p_HAC < 0.05),
    Sig_HAC_10 = sum(p_HAC < 0.10),
    .groups = "drop"
  )

summ_labels <- c("1: DXY $\\rightarrow$ Credit",
                  "2: Credit $\\rightarrow$ DXY",
                  "3: DXY $\\leftrightarrow$ Quantities")
summ$Battery <- summ_labels

colnames(summ) <- c("Battery", "Tests", "Sig. 5\\% (cl.)", "Sig. 10\\% (cl.)",
                     "Sig. 5\\% (HAC)", "Sig. 10\\% (HAC)")

ss_xt <- xtable(summ,
                caption = "Summary: Significant Dollar Granger Results by Battery",
                label   = "tab:summary_dollar",
                digits  = c(0, 0, 0, 0, 0, 0, 0))
print(ss_xt, file = "table_summary_dollar.tex",
      include.rownames = FALSE, sanitize.text.function = identity,
      booktabs = TRUE, floating = TRUE, table.placement = "H",
      comment = FALSE)


df_nocovid <- df %>%
  filter(Date < as.Date("2020-02-01") | Date > as.Date("2020-05-31"))

cat("\nCOVID sensitivity sample: N =", nrow(df_nocovid), "\n")

covid_pairs <- list(
  c("lr_DXY", "d_log_HYG_LQD"),
  c("lr_DXY", "d_log_HYG_SHV"),
  c("lr_DXY", "lr_EMB")
)

covid_full <- map_dfr(covid_pairs, function(pair) {
  map_dfr(c(1, 5, 10), function(lag) {
    run_granger(df, pair[1], pair[2], lag) %>% mutate(Sample = "Full")
  })
})

covid_excl <- map_dfr(covid_pairs, function(pair) {
  map_dfr(c(1, 5, 10), function(lag) {
    run_granger(df_nocovid, pair[1], pair[2], lag) %>% mutate(Sample = "Excl. COVID")
  })
})

covid_comp <- bind_rows(covid_full, covid_excl)
write_csv(covid_comp, "results_covid_dollar.csv")

covid_tab <- covid_comp %>%
  mutate(Direction = paste0(clean_name(X), " $\\rightarrow$ ", clean_name(Y)),
         Sig_hac = sig_stars(p_HAC)) %>%
  select(Sample, Direction, Lag, N, F_HAC, p_HAC, Sig_hac)

colnames(covid_tab) <- c("Sample", "Direction", "Lag", "$N$",
                          "$F$ (HAC)", "$p$ (HAC)", "")

ct_xt <- xtable(covid_tab,
                caption = "COVID Sensitivity: DXY $\\rightarrow$ Credit (HAC-Robust)",
                label   = "tab:covid_dollar",
                digits  = c(0,0,0,0,0, 3,4,0))
print(ct_xt, file = "table_covid_dollar.tex",
      include.rownames = FALSE, sanitize.text.function = identity,
      booktabs = TRUE, floating = TRUE, table.placement = "H",
      scalebox = 0.75, comment = FALSE)

cat("\n=== Dollar analysis complete ===\n")
cat("Tables: table_adf_dollar.tex, table_sumstats_dollar.tex, table_corr_dollar.tex\n")
cat("        table_granger_dxy_credit.tex, table_granger_credit_dxy.tex\n")
cat("        table_granger_dxy_qty.tex, table_summary_dollar.tex\n")
cat("        table_covid_dollar.tex\n")
cat("Figures: fig_ts_dollar_credit.pdf, fig_scatter_dxy.pdf\n")
cat("         fig_rolling_corr_dxy.pdf, fig_rolling_granger_dxy.pdf\n")
cat("CSV:     results_granger_dollar.csv, results_rolling_granger_dollar.csv\n")
cat("         results_covid_dollar.csv\n")
