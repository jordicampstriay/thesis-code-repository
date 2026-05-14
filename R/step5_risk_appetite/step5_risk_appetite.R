###############################################################################
#  STEP 5 — RISK APPETITE EXTENSION (Section 5.6)
#  Thesis: "The Plumbing of Liquidity"
#  Author: Jordi Camps Triay
#
#  This script extends the funding-to-credit transmission analysis using:
#    (A) ETF-based credit measures as alternatives to Moody's OAS spreads
#    (B) JPY/USD 3-month cross-currency basis as a global dollar funding proxy
#
#  Input : ../data_master.csv, ../cross_currency_basis.csv
#  Output: figures, LaTeX tables, CSV results — all in this folder
###############################################################################

# ── PACKAGES ─────────────────────────────────────────────────────────────────
required_pkgs <- c("readr", "dplyr", "tidyr", "zoo", "tseries", "lmtest",
                   "vars", "urca", "ggplot2", "scales", "xtable",
                   "gridExtra", "sandwich", "strucchange")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

WD <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) getwd()
)
setwd(WD)

cat("================================================================\n")
cat("  STEP 5 — RISK APPETITE EXTENSION\n")
cat("  Working directory:", WD, "\n")
cat("================================================================\n\n")


###############################################################################
#  PART A — ETF-BASED CREDIT MEASURES
###############################################################################

cat("================================================================\n")
cat("  PART A: ETF-BASED CREDIT MEASURES\n")
cat("================================================================\n\n")

# ── Load master data ─────────────────────────────────────────────────────────
df <- read_csv(file.path(WD, "..", "data_master.csv"), show_col_types = FALSE)
df$Date <- as.Date(df$Date)

# ── Construct ETF-based variables (post-2018) ────────────────────────────────
df_etf <- df %>%
  filter(Date >= as.Date("2018-04-01")) %>%
  arrange(Date) %>%
  mutate(
    # Log returns (daily, percent)
    r_HYG = c(NA, 100 * diff(log(HYG))),
    r_LQD = c(NA, 100 * diff(log(LQD))),
    r_IEF = c(NA, 100 * diff(log(IEF))),
    r_TLT = c(NA, 100 * diff(log(TLT))),
    r_EMB = c(NA, 100 * diff(log(EMB))),

    # Duration-hedged credit returns
    #   HYG – β*IEF  isolates credit component (strip duration)
    #   LQD – β*IEF  same for IG

    # Log ratio (credit risk appetite proxies)
    log_HYG_IEF = log(HYG / IEF),
    log_LQD_IEF = log(LQD / IEF),
    log_HYG_LQD = log(HYG / LQD),

    # First differences of log ratios (stationary)
    d_log_HYG_IEF = c(NA, diff(log_HYG_IEF)),
    d_log_LQD_IEF = c(NA, diff(log_LQD_IEF)),
    d_log_HYG_LQD = c(NA, diff(log_HYG_LQD))
  )

cat(sprintf("  ETF sample: %d rows, %s to %s\n",
            nrow(df_etf), min(df_etf$Date), max(df_etf$Date)))

# ── ADF tests on ETF variables ──────────────────────────────────────────────
etf_vars <- c("r_HYG", "r_LQD", "r_IEF",
              "d_log_HYG_IEF", "d_log_LQD_IEF", "d_log_HYG_LQD")

cat("\n  ADF tests (ETF variables):\n")
adf_results <- data.frame()
for (v in etf_vars) {
  x <- df_etf[[v]][!is.na(df_etf[[v]])]
  if (length(x) < 50) next
  adf <- adf.test(x, alternative = "stationary")
  cat(sprintf("    %s: ADF = %.3f, p = %.4f %s\n",
              v, adf$statistic, adf$p.value,
              ifelse(adf$p.value < 0.05, "(stationary)", "(NON-stationary)")))
  adf_results <- rbind(adf_results, data.frame(
    Variable = v, ADF = round(adf$statistic, 3),
    p_value = round(adf$p.value, 4),
    Stationary = adf$p.value < 0.05,
    stringsAsFactors = FALSE
  ))
}

# ── Correlation matrix: funding vs ETF credit ────────────────────────────────
cat("\n  Correlations: funding variables vs ETF credit measures\n")

corr_vars <- c("SOFR_EFFR", "d_Reserves", "d_TGA",
               "d_Aaa", "d_Baa", "d_Baa_Aaa",
               "r_HYG", "r_LQD", "d_log_HYG_IEF",
               "d_log_LQD_IEF", "d_log_HYG_LQD")
corr_data <- df_etf %>% dplyr::select(all_of(corr_vars)) %>% drop_na()
corr_mat <- cor(corr_data)

cat("\n  Key correlations:\n")
pairs <- list(
  c("SOFR_EFFR", "r_HYG"), c("SOFR_EFFR", "d_log_HYG_IEF"),
  c("d_Reserves", "r_HYG"), c("d_Reserves", "d_log_HYG_IEF"),
  c("r_HYG", "d_Baa"), c("d_log_HYG_IEF", "d_Baa_Aaa"),
  c("d_log_HYG_LQD", "d_Baa_Aaa")
)
for (pr in pairs) {
  r <- corr_mat[pr[1], pr[2]]
  cat(sprintf("    %s vs %s: r = %.4f\n", pr[[1]], pr[[2]], r))
}


###############################################################################
#  A.1 — GRANGER CAUSALITY: FUNDING → ETF CREDIT
###############################################################################

cat("\n================================================================\n")
cat("  A.1: GRANGER CAUSALITY — FUNDING → ETF CREDIT\n")
cat("================================================================\n\n")

run_granger <- function(data, x_name, y_name, lags = c(1, 5, 10), min_obs = 100) {
  results <- data.frame()
  for (lag in lags) {
    cols <- c(y_name, x_name)
    sub <- data %>% dplyr::select(all_of(cols)) %>% drop_na()
    if (nrow(sub) < min_obs) next
    tryCatch({
      gt <- grangertest(as.formula(paste(y_name, "~", x_name)), order = lag, data = sub)
      f_val <- gt$F[2]
      p_val <- gt$`Pr(>F)`[2]
      sig <- ifelse(p_val < 0.01, "***",
             ifelse(p_val < 0.05, " **",
             ifelse(p_val < 0.10, "  *", "   ")))
      results <- rbind(results, data.frame(
        X = x_name, Y = y_name, Lag = lag, N = nrow(sub),
        F_stat = round(f_val, 3), p_value = round(p_val, 4),
        Sig = sig, stringsAsFactors = FALSE
      ))
    }, error = function(e) NULL)
  }
  return(results)
}

# Funding → ETF credit tests
funding_vars <- c("SOFR_EFFR", "d_Reserves", "d_TGA")
etf_credit_vars <- c("r_HYG", "r_LQD", "d_log_HYG_IEF",
                      "d_log_LQD_IEF", "d_log_HYG_LQD")

granger_etf <- data.frame()
for (xv in funding_vars) {
  for (yv in etf_credit_vars) {
    granger_etf <- rbind(granger_etf,
                         run_granger(df_etf, xv, yv, lags = c(1, 5, 10)))
  }
}

cat("  Results:\n")
print(granger_etf %>% filter(Sig != "   "), row.names = FALSE)
cat(sprintf("\n  Total tests: %d, Significant at 10%%: %d\n",
            nrow(granger_etf),
            sum(granger_etf$Sig != "   ")))

# Reverse causality: ETF credit → Funding
cat("\n  Reverse causality tests:\n")
granger_rev <- data.frame()
for (xv in etf_credit_vars) {
  for (yv in funding_vars) {
    granger_rev <- rbind(granger_rev,
                         run_granger(df_etf, xv, yv, lags = c(1, 5)))
  }
}
rev_sig <- granger_rev %>% filter(Sig != "   ")
if (nrow(rev_sig) > 0) {
  print(rev_sig, row.names = FALSE)
} else {
  cat("    No significant reverse causality\n")
}


###############################################################################
#  A.2 — VAR: FUNDING → ETF CREDIT (HYG/IEF ratio)
###############################################################################

cat("\n================================================================\n")
cat("  A.2: VAR — FUNDING → ETF CREDIT\n")
cat("================================================================\n\n")

# VAR E: d_TGA, d_Reserves, SOFR_EFFR, d_log_LQD_IEF, d_log_HYG_IEF
# Ordering: funding (most exogenous) → IG credit → HY credit
var_cols_e <- c("d_TGA", "d_Reserves", "SOFR_EFFR",
                "d_log_LQD_IEF", "d_log_HYG_IEF")
var_data_e <- df_etf %>% dplyr::select(all_of(var_cols_e)) %>% drop_na()
cat(sprintf("  VAR E observations: %d\n", nrow(var_data_e)))

lag_sel_e <- VARselect(var_data_e, lag.max = 20, type = "const")
cat("  Lag selection:\n")
print(lag_sel_e$selection)
opt_lag_e <- lag_sel_e$selection["AIC(n)"]
cat(sprintf("  Selected lag: %d (AIC)\n\n", opt_lag_e))

var_e <- VAR(var_data_e, p = opt_lag_e, type = "const")

roots_e <- roots(var_e)
cat(sprintf("  Max eigenvalue modulus: %.4f (Stable: %s)\n",
            max(roots_e), ifelse(max(roots_e) < 1, "YES", "NO")))

cat("\n  R-squared (adj):\n")
for (eq_name in names(var_e$varresult)) {
  r2 <- summary(var_e$varresult[[eq_name]])$adj.r.squared
  cat(sprintf("    %s: %.4f\n", eq_name, r2))
}

# Granger causality within VAR
cat("\n  Granger causality in VAR framework:\n")
for (cause_var in c("d_TGA", "d_Reserves", "SOFR_EFFR")) {
  gc <- causality(var_e, cause = cause_var)
  cat(sprintf("    %s → others: F = %.3f, p = %.4f %s\n",
              cause_var, gc$Granger$statistic, gc$Granger$p.value,
              ifelse(gc$Granger$p.value < 0.05, "**", "")))
}

# Portmanteau
pt_e <- serial.test(var_e, lags.pt = 20, type = "PT.asymptotic")
cat(sprintf("\n  Portmanteau (20 lags): Chi-sq = %.1f, p = %.4f\n",
            pt_e$serial$statistic, pt_e$serial$p.value))

# IRFs
n_ahead <- 30
n_boot <- 500
irf_e <- irf(var_e, impulse = NULL, response = NULL,
              n.ahead = n_ahead, boot = TRUE, runs = n_boot,
              ci = 0.95, ortho = TRUE)

make_irf_panel_df <- function(irf_obj, impulse, response, label) {
  idx <- which(names(irf_obj$irf) == impulse)
  irf_mat <- irf_obj$irf[[impulse]]
  low_mat <- irf_obj$Lower[[impulse]]
  up_mat  <- irf_obj$Upper[[impulse]]
  col_idx <- which(colnames(irf_mat) == response)
  data.frame(
    Horizon = 0:n_ahead,
    IRF     = irf_mat[, col_idx],
    Lower   = low_mat[, col_idx],
    Upper   = up_mat[, col_idx],
    Panel   = label,
    stringsAsFactors = FALSE
  )
}

# Panel IRF plot: funding → ETF credit
panel_e_df <- rbind(
  make_irf_panel_df(irf_e, "SOFR_EFFR",  "d_log_HYG_IEF", "SOFR-EFFR → Δlog(HYG/IEF)"),
  make_irf_panel_df(irf_e, "SOFR_EFFR",  "d_log_LQD_IEF", "SOFR-EFFR → Δlog(LQD/IEF)"),
  make_irf_panel_df(irf_e, "d_Reserves", "d_log_HYG_IEF", "ΔReserves → Δlog(HYG/IEF)"),
  make_irf_panel_df(irf_e, "d_Reserves", "d_log_LQD_IEF", "ΔReserves → Δlog(LQD/IEF)"),
  make_irf_panel_df(irf_e, "d_TGA",      "d_log_HYG_IEF", "ΔTGA → Δlog(HYG/IEF)"),
  make_irf_panel_df(irf_e, "d_TGA",      "d_log_LQD_IEF", "ΔTGA → Δlog(LQD/IEF)")
)

panel_e_df$Panel <- factor(panel_e_df$Panel,
  levels = c("SOFR-EFFR → Δlog(HYG/IEF)", "SOFR-EFFR → Δlog(LQD/IEF)",
             "ΔReserves → Δlog(HYG/IEF)", "ΔReserves → Δlog(LQD/IEF)",
             "ΔTGA → Δlog(HYG/IEF)",      "ΔTGA → Δlog(LQD/IEF)"))

p_irf_etf <- ggplot(panel_e_df, aes(x = Horizon)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange3", alpha = 0.15) +
  geom_line(aes(y = IRF), color = "darkorange3", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  facet_wrap(~ Panel, scales = "free_y", ncol = 2) +
  labs(title = "Impulse Responses — Funding Shocks to ETF Credit Measures (Post-2018)",
       subtitle = "HYG/IEF = high-yield risk appetite; LQD/IEF = investment-grade risk appetite",
       x = "Horizon (business days)", y = "Response") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        strip.text = element_text(face = "bold", size = 7.5),
        panel.grid.minor = element_blank())

ggsave("fig_irf_panel_etf.pdf", p_irf_etf, width = 8, height = 7, device = "pdf")
cat("  Saved: fig_irf_panel_etf.pdf\n")

# FEVD
fevd_e <- fevd(var_e, n.ahead = 30)
cat("\n  FEVD for ETF credit measures:\n")
for (resp in c("d_log_HYG_IEF", "d_log_LQD_IEF")) {
  cat(sprintf("  %s:\n", resp))
  fm <- fevd_e[[resp]]
  for (h in c(1, 5, 10, 20, 30)) {
    if (h <= nrow(fm)) {
      cat(sprintf("    h=%2d: TGA=%.1f%% Res=%.1f%% SOFR=%.1f%% LQD/IEF=%.1f%% HYG/IEF=%.1f%%\n",
                  h,
                  100*fm[h, "d_TGA"], 100*fm[h, "d_Reserves"],
                  100*fm[h, "SOFR_EFFR"],
                  100*fm[h, "d_log_LQD_IEF"], 100*fm[h, "d_log_HYG_IEF"]))
    }
  }
}

# FEVD plot
fevd_e_plot <- data.frame()
for (resp in c("d_log_HYG_IEF", "d_log_LQD_IEF")) {
  fm <- fevd_e[[resp]]
  for (h in 1:nrow(fm)) {
    for (src in var_cols_e) {
      fevd_e_plot <- rbind(fevd_e_plot, data.frame(
        Response = resp, Horizon = h, Source = src,
        Share = 100 * fm[h, src], stringsAsFactors = FALSE
      ))
    }
  }
}
fevd_e_plot$Source <- factor(fevd_e_plot$Source, levels = var_cols_e,
  labels = c("TGA", "Reserves", "SOFR-EFFR", "LQD/IEF", "HYG/IEF"))
fevd_e_plot$Response <- factor(fevd_e_plot$Response,
  levels = c("d_log_HYG_IEF", "d_log_LQD_IEF"),
  labels = c("Var. of Δlog(HYG/IEF)", "Var. of Δlog(LQD/IEF)"))

p_fevd_etf <- ggplot(fevd_e_plot, aes(x = Horizon, y = Share, fill = Source)) +
  geom_area(alpha = 0.8) +
  facet_wrap(~ Response, ncol = 2) +
  scale_fill_brewer(palette = "Set2", name = "Shock source") +
  labs(x = "Horizon (business days)", y = "Share of forecast error variance (%)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave("fig_fevd_etf.pdf", p_fevd_etf, width = 8, height = 4, device = "pdf")
cat("  Saved: fig_fevd_etf.pdf\n")


###############################################################################
#  A.3 — COMPARISON: ETF vs MOODY'S RESULTS
###############################################################################

cat("\n================================================================\n")
cat("  A.3: ETF vs MOODY'S COMPARISON\n")
cat("================================================================\n\n")

# VAR F: HYG/LQD ratio as market-based quality differential
# Compare with Baa-Aaa (Moody's quality diff)
var_cols_f <- c("d_TGA", "d_Reserves", "SOFR_EFFR", "d_log_HYG_LQD")
var_data_f <- df_etf %>% dplyr::select(all_of(var_cols_f)) %>% drop_na()
cat(sprintf("  VAR F (HYG/LQD quality diff) observations: %d\n", nrow(var_data_f)))

lag_sel_f <- VARselect(var_data_f, lag.max = 20, type = "const")
opt_lag_f <- lag_sel_f$selection["AIC(n)"]
cat(sprintf("  Selected lag: %d (AIC)\n", opt_lag_f))

var_f <- VAR(var_data_f, p = opt_lag_f, type = "const")
cat(sprintf("  Max eigenvalue: %.4f (Stable: %s)\n",
            max(roots(var_f)), ifelse(max(roots(var_f)) < 1, "YES", "NO")))

cat("\n  R-squared:\n")
for (eq_name in names(var_f$varresult)) {
  cat(sprintf("    %s: %.4f\n", eq_name,
              summary(var_f$varresult[[eq_name]])$adj.r.squared))
}

# Granger within VAR
cat("\n  Granger in VAR:\n")
for (cause_var in c("d_TGA", "d_Reserves", "SOFR_EFFR")) {
  gc <- causality(var_f, cause = cause_var)
  cat(sprintf("    %s → others: F = %.3f, p = %.4f %s\n",
              cause_var, gc$Granger$statistic, gc$Granger$p.value,
              ifelse(gc$Granger$p.value < 0.05, "**", "")))
}

# IRF: funding → HYG/LQD
irf_f <- irf(var_f, impulse = NULL, response = NULL,
              n.ahead = n_ahead, boot = TRUE, runs = n_boot,
              ci = 0.95, ortho = TRUE)

panel_f_df <- rbind(
  make_irf_panel_df(irf_f, "SOFR_EFFR",  "d_log_HYG_LQD", "SOFR-EFFR → Δlog(HYG/LQD)"),
  make_irf_panel_df(irf_f, "d_Reserves", "d_log_HYG_LQD", "ΔReserves → Δlog(HYG/LQD)"),
  make_irf_panel_df(irf_f, "d_TGA",      "d_log_HYG_LQD", "ΔTGA → Δlog(HYG/LQD)")
)

p_irf_quality <- ggplot(panel_f_df, aes(x = Horizon)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "purple4", alpha = 0.15) +
  geom_line(aes(y = IRF), color = "purple4", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  facet_wrap(~ Panel, ncol = 3) +
  labs(title = "Funding Shocks → HYG/LQD Ratio (Market-Based Quality Differential)",
       x = "Horizon (business days)", y = "Response") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        strip.text = element_text(face = "bold", size = 8),
        panel.grid.minor = element_blank())

ggsave("fig_irf_quality_etf.pdf", p_irf_quality, width = 9, height = 3.5, device = "pdf")
cat("  Saved: fig_irf_quality_etf.pdf\n")

# FEVD comparison
fevd_f <- fevd(var_f, n.ahead = 30)
cat("\n  FEVD for Δlog(HYG/LQD):\n")
fm <- fevd_f[["d_log_HYG_LQD"]]
for (h in c(1, 5, 10, 20, 30)) {
  cat(sprintf("    h=%2d: TGA=%.1f%% Res=%.1f%% SOFR=%.1f%% Own=%.1f%%\n",
              h, 100*fm[h,"d_TGA"], 100*fm[h,"d_Reserves"],
              100*fm[h,"SOFR_EFFR"], 100*fm[h,"d_log_HYG_LQD"]))
}


###############################################################################
#  PART B — CROSS-CURRENCY BASIS (JPY/USD)
###############################################################################

cat("\n\n================================================================\n")
cat("  PART B: CROSS-CURRENCY BASIS (JPY/USD)\n")
cat("================================================================\n\n")

# ── Load and merge ───────────────────────────────────────────────────────────
ccb <- read_csv(file.path(WD, "..", "cross_currency_basis.csv"), show_col_types = FALSE)
ccb$Date <- as.Date(ccb$Date)

# Merge with master data on overlapping dates
df_ccb <- df %>%
  inner_join(ccb %>% dplyr::select(Date, JPY), by = "Date") %>%
  arrange(Date) %>%
  mutate(
    d_JPY_basis = c(NA, diff(JPY)),
    # Negative basis = dollar funding premium
    # More negative = more expensive to borrow dollars via FX swaps
    JPY_basis_abs = abs(JPY)
  )

cat(sprintf("  Merged sample: %d rows, %s to %s\n",
            nrow(df_ccb), min(df_ccb$Date), max(df_ccb$Date)))
cat(sprintf("  JPY basis: mean = %.2f bp, sd = %.2f bp, min = %.2f, max = %.2f\n",
            mean(df_ccb$JPY, na.rm=TRUE), sd(df_ccb$JPY, na.rm=TRUE),
            min(df_ccb$JPY, na.rm=TRUE), max(df_ccb$JPY, na.rm=TRUE)))

# ── ADF tests ────────────────────────────────────────────────────────────────
cat("\n  ADF tests:\n")
for (v in c("JPY", "d_JPY_basis")) {
  x <- df_ccb[[v]][!is.na(df_ccb[[v]])]
  if (length(x) < 30) next
  adf <- adf.test(x, alternative = "stationary")
  cat(sprintf("    %s: ADF = %.3f, p = %.4f %s\n",
              v, adf$statistic, adf$p.value,
              ifelse(adf$p.value < 0.05, "(stationary)", "(non-stationary)")))
}

# ── Correlations ─────────────────────────────────────────────────────────────
cat("\n  Correlations (JPY basis vs funding/credit):\n")
corr_pairs_ccb <- list(
  c("JPY", "SOFR_EFFR"), c("JPY", "Baa_spread"), c("JPY", "Aaa_spread"),
  c("d_JPY_basis", "SOFR_EFFR"), c("d_JPY_basis", "d_Baa"),
  c("d_JPY_basis", "d_Aaa"), c("d_JPY_basis", "d_Reserves"),
  c("d_JPY_basis", "d_Baa_Aaa")
)
for (pr in corr_pairs_ccb) {
  sub <- df_ccb %>% dplyr::select(all_of(pr)) %>% drop_na()
  if (nrow(sub) < 30) next
  ct <- cor.test(sub[[1]], sub[[2]])
  cat(sprintf("    %s vs %s: r = %.4f (p = %.4f, n = %d) %s\n",
              pr[1], pr[2], ct$estimate, ct$p.value, nrow(sub),
              ifelse(ct$p.value < 0.05, "**", "")))
}


###############################################################################
#  B.1 — GRANGER CAUSALITY: CCB
###############################################################################

cat("\n================================================================\n")
cat("  B.1: GRANGER CAUSALITY — CROSS-CURRENCY BASIS\n")
cat("================================================================\n\n")

# With only ~520 obs, use lags 1 and 5 only
ccb_granger <- data.frame()

# (i) Funding → d_JPY_basis (does domestic funding stress widen the basis?)
for (xv in c("SOFR_EFFR", "d_Reserves", "d_TGA")) {
  ccb_granger <- rbind(ccb_granger,
                       run_granger(df_ccb, xv, "d_JPY_basis", lags = c(1, 5), min_obs = 50))
}

# (ii) d_JPY_basis → credit (does dollar funding premium predict credit?)
for (yv in c("d_Baa", "d_Aaa", "d_Baa_Aaa")) {
  ccb_granger <- rbind(ccb_granger,
                       run_granger(df_ccb, "d_JPY_basis", yv, lags = c(1, 5), min_obs = 50))
}

# (iii) JPY level → credit (does the level of dollar premium predict?)
for (yv in c("d_Baa", "d_Aaa", "d_Baa_Aaa")) {
  ccb_granger <- rbind(ccb_granger,
                       run_granger(df_ccb, "JPY", yv, lags = c(1, 5), min_obs = 50))
}

# (iv) Reverse: credit → d_JPY_basis
for (xv in c("d_Baa", "d_Aaa", "d_Baa_Aaa")) {
  ccb_granger <- rbind(ccb_granger,
                       run_granger(df_ccb, xv, "d_JPY_basis", lags = c(1, 5), min_obs = 50))
}

# (v) Funding → credit within this 2-year window (benchmark)
for (xv in c("SOFR_EFFR", "d_Reserves", "d_TGA")) {
  for (yv in c("d_Baa", "d_Aaa", "d_Baa_Aaa")) {
    ccb_granger <- rbind(ccb_granger,
                         run_granger(df_ccb, xv, yv, lags = c(1, 5), min_obs = 50))
  }
}

cat("  All Granger results (cross-currency basis sample):\n")
print(ccb_granger, row.names = FALSE)

cat(sprintf("\n  Significant at 10%%: %d / %d tests\n",
            sum(ccb_granger$Sig != "   "), nrow(ccb_granger)))

sig_ccb <- ccb_granger %>% filter(Sig != "   ")
if (nrow(sig_ccb) > 0) {
  cat("\n  Significant results:\n")
  print(sig_ccb, row.names = FALSE)
}


###############################################################################
#  B.2 — BIVARIATE VAR: SOFR-EFFR & JPY BASIS
###############################################################################

cat("\n================================================================\n")
cat("  B.2: VAR — SOFR-EFFR & JPY BASIS\n")
cat("================================================================\n\n")

# Small sample → keep it parsimonious. Try bivariate and trivariate.

# Bivariate: SOFR_EFFR ↔ d_JPY_basis
var_cols_g1 <- c("SOFR_EFFR", "d_JPY_basis")
var_data_g1 <- df_ccb %>% dplyr::select(all_of(var_cols_g1)) %>% drop_na()
cat(sprintf("  Bivariate VAR (SOFR-EFFR, d_JPY_basis): %d obs\n", nrow(var_data_g1)))

lag_sel_g1 <- VARselect(var_data_g1, lag.max = 15, type = "const")
# Use BIC for small sample (more conservative)
opt_lag_g1 <- lag_sel_g1$selection["SC(n)"]
cat(sprintf("  BIC lag: %d, AIC lag: %d — using BIC for small sample\n",
            opt_lag_g1, lag_sel_g1$selection["AIC(n)"]))

var_g1 <- VAR(var_data_g1, p = max(1, opt_lag_g1), type = "const")
cat(sprintf("  Max eigenvalue: %.4f\n", max(roots(var_g1))))

cat("\n  R-squared:\n")
for (eq in names(var_g1$varresult)) {
  cat(sprintf("    %s: %.4f\n", eq,
              summary(var_g1$varresult[[eq]])$adj.r.squared))
}

# IRF
irf_g1 <- irf(var_g1, impulse = NULL, response = NULL,
               n.ahead = 20, boot = TRUE, runs = n_boot,
               ci = 0.95, ortho = TRUE)

# Both directions
irf_g1_df <- rbind(
  data.frame(
    Horizon = 0:20,
    IRF = irf_g1$irf$SOFR_EFFR[, "d_JPY_basis"],
    Lower = irf_g1$Lower$SOFR_EFFR[, "d_JPY_basis"],
    Upper = irf_g1$Upper$SOFR_EFFR[, "d_JPY_basis"],
    Panel = "SOFR-EFFR → Δ(JPY basis)", stringsAsFactors = FALSE
  ),
  data.frame(
    Horizon = 0:20,
    IRF = irf_g1$irf$d_JPY_basis[, "SOFR_EFFR"],
    Lower = irf_g1$Lower$d_JPY_basis[, "SOFR_EFFR"],
    Upper = irf_g1$Upper$d_JPY_basis[, "SOFR_EFFR"],
    Panel = "Δ(JPY basis) → SOFR-EFFR", stringsAsFactors = FALSE
  )
)

p_irf_ccb_bivar <- ggplot(irf_g1_df, aes(x = Horizon)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkred", alpha = 0.15) +
  geom_line(aes(y = IRF), color = "darkred", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ Panel, scales = "free_y", ncol = 2) +
  labs(title = "Bivariate IRFs: SOFR-EFFR and JPY/USD Cross-Currency Basis",
       subtitle = paste0("Sample: ", min(df_ccb$Date), " to ", max(df_ccb$Date),
                         " (", nrow(var_data_g1), " obs, BIC lag = ", max(1, opt_lag_g1), ")"),
       x = "Horizon (business days)", y = "Response") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave("fig_irf_ccb_bivariate.pdf", p_irf_ccb_bivar, width = 8, height = 3.5, device = "pdf")
cat("  Saved: fig_irf_ccb_bivariate.pdf\n")


###############################################################################
#  B.3 — TRIVARIATE VAR: JPY BASIS, FUNDING, CREDIT
###############################################################################

cat("\n================================================================\n")
cat("  B.3: TRIVARIATE VAR — BASIS, FUNDING, CREDIT\n")
cat("================================================================\n\n")

# d_JPY_basis → SOFR_EFFR → d_Baa
# Ordering: JPY basis (global, most exogenous) → SOFR-EFFR → d_Baa
var_cols_g2 <- c("d_JPY_basis", "SOFR_EFFR", "d_Baa")
var_data_g2 <- df_ccb %>% dplyr::select(all_of(var_cols_g2)) %>% drop_na()
cat(sprintf("  Trivariate VAR observations: %d\n", nrow(var_data_g2)))

lag_sel_g2 <- VARselect(var_data_g2, lag.max = 15, type = "const")
opt_lag_g2 <- lag_sel_g2$selection["SC(n)"]
cat(sprintf("  BIC lag: %d\n", opt_lag_g2))

var_g2 <- VAR(var_data_g2, p = max(1, opt_lag_g2), type = "const")
cat(sprintf("  Max eigenvalue: %.4f (Stable: %s)\n",
            max(roots(var_g2)), ifelse(max(roots(var_g2)) < 1, "YES", "NO")))

cat("\n  R-squared:\n")
for (eq in names(var_g2$varresult)) {
  cat(sprintf("    %s: %.4f\n", eq,
              summary(var_g2$varresult[[eq]])$adj.r.squared))
}

# Granger
cat("\n  Granger in trivariate VAR:\n")
for (cv in var_cols_g2) {
  gc <- causality(var_g2, cause = cv)
  cat(sprintf("    %s → others: F = %.3f, p = %.4f %s\n",
              cv, gc$Granger$statistic, gc$Granger$p.value,
              ifelse(gc$Granger$p.value < 0.05, "**", "")))
}

# IRF
irf_g2 <- irf(var_g2, impulse = NULL, response = NULL,
               n.ahead = 20, boot = TRUE, runs = n_boot,
               ci = 0.95, ortho = TRUE)

panel_g2_df <- rbind(
  data.frame(Horizon = 0:20,
    IRF = irf_g2$irf$d_JPY_basis[, "SOFR_EFFR"],
    Lower = irf_g2$Lower$d_JPY_basis[, "SOFR_EFFR"],
    Upper = irf_g2$Upper$d_JPY_basis[, "SOFR_EFFR"],
    Panel = "Δ(JPY basis) → SOFR-EFFR", stringsAsFactors = FALSE),
  data.frame(Horizon = 0:20,
    IRF = irf_g2$irf$d_JPY_basis[, "d_Baa"],
    Lower = irf_g2$Lower$d_JPY_basis[, "d_Baa"],
    Upper = irf_g2$Upper$d_JPY_basis[, "d_Baa"],
    Panel = "Δ(JPY basis) → ΔBaa", stringsAsFactors = FALSE),
  data.frame(Horizon = 0:20,
    IRF = irf_g2$irf$SOFR_EFFR[, "d_Baa"],
    Lower = irf_g2$Lower$SOFR_EFFR[, "d_Baa"],
    Upper = irf_g2$Upper$SOFR_EFFR[, "d_Baa"],
    Panel = "SOFR-EFFR → ΔBaa", stringsAsFactors = FALSE)
)

p_irf_ccb_tri <- ggplot(panel_g2_df, aes(x = Horizon)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkred", alpha = 0.15) +
  geom_line(aes(y = IRF), color = "darkred", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ Panel, scales = "free_y", ncol = 3) +
  labs(title = "Trivariate IRFs: JPY Basis → SOFR-EFFR → Baa Spread",
       subtitle = paste0("Cholesky ordering: Δ(JPY basis) → SOFR-EFFR → ΔBaa | ",
                         nrow(var_data_g2), " obs, BIC lag = ", max(1, opt_lag_g2)),
       x = "Horizon (business days)", y = "Response") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        strip.text = element_text(face = "bold", size = 8),
        panel.grid.minor = element_blank())

ggsave("fig_irf_ccb_trivariate.pdf", p_irf_ccb_tri, width = 9, height = 3.5, device = "pdf")
cat("  Saved: fig_irf_ccb_trivariate.pdf\n")

# FEVD
fevd_g2 <- fevd(var_g2, n.ahead = 20)
cat("\n  FEVD for d_Baa (trivariate CCB VAR):\n")
fm <- fevd_g2[["d_Baa"]]
for (h in c(1, 5, 10, 20)) {
  cat(sprintf("    h=%2d: JPY_basis=%.1f%% SOFR=%.1f%% Own=%.1f%%\n",
              h, 100*fm[h,"d_JPY_basis"], 100*fm[h,"SOFR_EFFR"], 100*fm[h,"d_Baa"]))
}


###############################################################################
#  B.4 — ROLLING CORRELATION: JPY BASIS vs SOFR-EFFR
###############################################################################

cat("\n")

df_ccb_roll <- df_ccb %>%
  filter(!is.na(JPY), !is.na(SOFR_EFFR)) %>%
  arrange(Date) %>%
  mutate(
    roll_corr_jpy_sofr = zoo::rollapply(
      cbind(JPY, SOFR_EFFR), width = 60,
      FUN = function(x) cor(x[,1], x[,2], use = "complete.obs"),
      by.column = FALSE, fill = NA, align = "right"
    ),
    roll_corr_djpy_dbaa = zoo::rollapply(
      cbind(d_JPY_basis, d_Baa), width = 60,
      FUN = function(x) cor(x[,1], x[,2], use = "complete.obs"),
      by.column = FALSE, fill = NA, align = "right"
    )
  )

p_roll_ccb <- ggplot(df_ccb_roll %>% filter(!is.na(roll_corr_jpy_sofr)),
       aes(x = Date)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  geom_line(aes(y = roll_corr_jpy_sofr, color = "JPY basis vs SOFR-EFFR (levels)"),
            linewidth = 0.6) +
  geom_line(aes(y = roll_corr_djpy_dbaa, color = "Δ(JPY basis) vs ΔBaa"),
            linewidth = 0.6) +
  scale_color_manual(values = c("JPY basis vs SOFR-EFFR (levels)" = "darkred",
                                 "Δ(JPY basis) vs ΔBaa" = "steelblue")) +
  labs(title = "60-Day Rolling Correlations: JPY/USD Basis vs Domestic Funding/Credit",
       x = NULL, y = "Correlation", color = NULL) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("fig_rolling_corr_ccb.pdf", p_roll_ccb, width = 8, height = 4, device = "pdf")
cat("  Saved: fig_rolling_corr_ccb.pdf\n")


###############################################################################
#  B.5 — TIME SERIES: JPY BASIS vs SOFR-EFFR
###############################################################################

p_ts_ccb <- ggplot(df_ccb, aes(x = Date)) +
  geom_line(aes(y = JPY, color = "JPY/USD 3M basis (bp)"), linewidth = 0.5) +
  geom_line(aes(y = SOFR_EFFR * 100, color = "SOFR-EFFR (bp)"), linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  scale_color_manual(values = c("JPY/USD 3M basis (bp)" = "darkred",
                                 "SOFR-EFFR (bp)" = "steelblue")) +
  labs(title = "JPY/USD Cross-Currency Basis and SOFR-EFFR Spread",
       subtitle = "Negative basis = dollar funding premium in FX swap market",
       x = NULL, y = "Basis points", color = NULL) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("fig_ts_ccb_vs_sofr.pdf", p_ts_ccb, width = 8, height = 4, device = "pdf")
cat("  Saved: fig_ts_ccb_vs_sofr.pdf\n")


###############################################################################
#  LATEX TABLES
###############################################################################

cat("\n================================================================\n")
cat("  GENERATING LATEX TABLES\n")
cat("================================================================\n\n")

# ── Table: Granger ETF results (significant only) ──
sig_etf_all <- rbind(
  granger_etf %>% mutate(Battery = "Funding → ETF credit"),
  granger_rev %>% mutate(Battery = "ETF credit → Funding (reverse)")
)

etf_tbl <- sig_etf_all %>%
  filter(Sig != "   " | Battery == "Funding → ETF credit") %>%
  dplyr::select(Battery, X, Y, Lag, N, F_stat, p_value, Sig)

etf_xt <- xtable(etf_tbl,
  caption = "Granger Causality: Funding Variables and ETF-Based Credit Measures (Post-2018)",
  label = "tab:granger_etf",
  digits = c(0, 0, 0, 0, 0, 0, 3, 4, 0))
print(etf_xt, file = "table_granger_etf.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp")
cat("  Saved: table_granger_etf.tex\n")

# ── Table: CCB Granger results ──
ccb_xt <- xtable(ccb_granger,
  caption = "Granger Causality: Cross-Currency Basis and Funding/Credit Variables (May 2024--May 2026)",
  label = "tab:granger_ccb",
  digits = c(0, 0, 0, 0, 0, 3, 4, 0))
print(ccb_xt, file = "table_granger_ccb.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp")
cat("  Saved: table_granger_ccb.tex\n")

# ── Table: VAR comparison ──
comp_df <- data.frame(
  Specification = c(
    "VAR A: Moody's (Post-2018)",
    "VAR E: ETF HYG/IEF, LQD/IEF (Post-2018)",
    "VAR F: ETF HYG/LQD quality diff (Post-2018)",
    "Trivariate: JPY basis, SOFR, Baa (2024-26)"
  ),
  Variables = c(
    "d_TGA, d_Res, SOFR-EFFR, d_Aaa, d_Baa",
    "d_TGA, d_Res, SOFR-EFFR, dlog(LQD/IEF), dlog(HYG/IEF)",
    "d_TGA, d_Res, SOFR-EFFR, dlog(HYG/LQD)",
    "d_JPY_basis, SOFR-EFFR, d_Baa"
  ),
  N = c(
    nrow(var_data_e),  # reusing same post-2018 sample size
    nrow(var_data_e),
    nrow(var_data_f),
    nrow(var_data_g2)
  ),
  Lags = c(opt_lag_e, opt_lag_e, opt_lag_f, max(1, opt_lag_g2)),
  stringsAsFactors = FALSE
)

comp_xt <- xtable(comp_df,
  caption = "VAR Specifications --- Risk Appetite Extension",
  label = "tab:var_extension",
  digits = c(0, 0, 0, 0, 0))
print(comp_xt, file = "table_var_extension.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp")
cat("  Saved: table_var_extension.tex\n")

# ── Save CSV ──
write_csv(rbind(
  granger_etf %>% mutate(Battery = "ETF_funding_to_credit"),
  granger_rev %>% mutate(Battery = "ETF_reverse"),
  ccb_granger %>% mutate(Battery = "CCB")
), "granger_extension_results.csv")
cat("  Saved: granger_extension_results.csv\n")


cat("\n================================================================\n")
cat("  STEP 5 COMPLETE\n")
cat("================================================================\n")
