###############################################################################
#  STEP 5 — THE DOLLAR AND RISK APPETITE (Section 5.6)
#  Thesis: "The Plumbing of Liquidity"
#  Author: Jordi Camps Triay
#
#  This script tests whether the US dollar (DXY) is the missing link in the
#  funding-to-credit transmission chain documented in Sections 5.3–5.5.
#
#  Structure:
#    Part A — The Dollar Channel: DXY → credit spreads, DXY ↔ funding
#    Part B — ETF Robustness: traded credit measures vs Moody's OAS
#    Part C — Cross-Currency Basis: global dollar funding (JPY/USD, 2 years)
#
#  All variables in this script are daily — no forward-filled weekly data.
#  This addresses the frequency mismatch limitation of the reserve/TGA
#  channels analysed in Sections 5.3–5.4 (see discussion in Section 5.7).
#
#  Input : ../data_master.csv, ../cross_currency_basis.csv
#  Output: figures, LaTeX tables, CSV — all written to this folder
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

WD <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
setwd(WD)

cat("================================================================\n")
cat("  STEP 5 — THE DOLLAR AND RISK APPETITE\n")
cat("  Working directory:", WD, "\n")
cat("================================================================\n\n")

# ── SHARED DATA & HELPERS ────────────────────────────────────────────────────
df <- read_csv(file.path(WD, "..", "data_master.csv"), show_col_types = FALSE)
df$Date <- as.Date(df$Date)

df_post <- df %>%
  filter(Date >= as.Date("2018-04-01")) %>%
  arrange(Date) %>%
  mutate(
    r_DXY       = c(NA, 100 * diff(log(DXY))),
    r_HYG       = c(NA, 100 * diff(log(HYG))),
    r_LQD       = c(NA, 100 * diff(log(LQD))),
    r_IEF       = c(NA, 100 * diff(log(IEF))),
    r_HYG_SHV   = c(NA, 100 * diff(log(HYG / SHV))),
    r_HYG_LQD   = c(NA, 100 * diff(log(HYG / LQD))),
    d_log_HYG_IEF = c(NA, diff(log(HYG / IEF))),
    d_log_LQD_IEF = c(NA, diff(log(LQD / IEF)))
  )

cat(sprintf("  Post-2018 sample: %d rows, %s to %s\n",
            nrow(df_post), min(df_post$Date), max(df_post$Date)))

n_ahead <- 30
n_boot  <- 500

run_granger <- function(data, x_name, y_name, lags = c(1, 5, 10), min_obs = 100) {
  results <- data.frame()
  for (lag in lags) {
    sub <- data %>% dplyr::select(all_of(c(y_name, x_name))) %>% drop_na()
    if (nrow(sub) < min_obs) next
    tryCatch({
      gt <- grangertest(as.formula(paste(y_name, "~", x_name)),
                        order = lag, data = sub)
      f_val <- gt$F[2]; p_val <- gt$`Pr(>F)`[2]
      sig <- ifelse(p_val < 0.01, "***",
             ifelse(p_val < 0.05, " **",
             ifelse(p_val < 0.10, "  *", "   ")))
      results <- rbind(results, data.frame(
        X = x_name, Y = y_name, Lag = lag, N = nrow(sub),
        F_stat = round(f_val, 3), p_value = round(p_val, 4),
        Sig = sig, stringsAsFactors = FALSE))
    }, error = function(e) NULL)
  }
  return(results)
}

make_irf_df <- function(irf_obj, impulse, response, label) {
  data.frame(
    Horizon = 0:n_ahead,
    IRF   = irf_obj$irf[[impulse]][, response],
    Lower = irf_obj$Lower[[impulse]][, response],
    Upper = irf_obj$Upper[[impulse]][, response],
    Panel = label, stringsAsFactors = FALSE)
}


###############################################################################
#  PART A — THE DOLLAR CHANNEL
###############################################################################

cat("================================================================\n")
cat("  PART A: THE DOLLAR CHANNEL\n")
cat("================================================================\n\n")

# ── A.1  Correlations ────────────────────────────────────────────────────────

cat("  A.1 — Contemporaneous correlations (r_DXY vs):\n")
corr_targets <- c("d_Baa", "d_Aaa", "d_Baa_Aaa",
                   "r_HYG_SHV", "r_HYG_LQD",
                   "SOFR_EFFR", "d_Reserves")
for (yv in corr_targets) {
  sub <- df_post %>% dplyr::select(r_DXY, all_of(yv)) %>% drop_na()
  ct <- cor.test(sub$r_DXY, sub[[yv]])
  sig <- ifelse(ct$p.value < 0.01, "***",
         ifelse(ct$p.value < 0.05, " **", ""))
  cat(sprintf("    r_DXY vs %-14s r = %7.4f  (p = %.4f) %s\n",
              yv, ct$estimate, ct$p.value, sig))
}

# ── A.2  Granger causality ───────────────────────────────────────────────────

cat("\n  A.2 — Granger causality:\n\n")

granger_dxy <- data.frame()

# DXY → credit
cat("  (i) DXY → credit spreads:\n")
for (yv in c("d_Baa", "d_Aaa", "d_Baa_Aaa")) {
  res <- run_granger(df_post, "r_DXY", yv)
  granger_dxy <- rbind(granger_dxy, res)
}
print(granger_dxy %>% filter(Sig != "   "), row.names = FALSE)

# DXY → ETF ratios
cat("\n  (ii) DXY → ETF credit ratios:\n")
for (yv in c("r_HYG_SHV", "r_HYG_LQD")) {
  res <- run_granger(df_post, "r_DXY", yv)
  granger_dxy <- rbind(granger_dxy, res)
}
print(granger_dxy %>% filter(Y %in% c("r_HYG_SHV","r_HYG_LQD"), Sig != "   "),
      row.names = FALSE)

# Reverse: credit → DXY
cat("\n  (iii) Credit → DXY (reverse):\n")
gc_rev <- data.frame()
for (xv in c("d_Baa", "d_Aaa", "r_HYG_SHV", "r_HYG_LQD")) {
  gc_rev <- rbind(gc_rev, run_granger(df_post, xv, "r_DXY"))
}
print(gc_rev %>% filter(Sig != "   "), row.names = FALSE)

# Reserves → DXY
cat("\n  (iv) Funding → DXY:\n")
gc_fund_dxy <- data.frame()
for (xv in c("d_Reserves", "d_TGA", "SOFR_EFFR")) {
  gc_fund_dxy <- rbind(gc_fund_dxy, run_granger(df_post, xv, "r_DXY"))
}
print(gc_fund_dxy, row.names = FALSE)

# DXY → Funding
cat("\n  (v) DXY → Funding:\n")
gc_dxy_fund <- data.frame()
for (yv in c("d_Reserves", "SOFR_EFFR", "d_TGA")) {
  gc_dxy_fund <- rbind(gc_dxy_fund, run_granger(df_post, "r_DXY", yv))
}
print(gc_dxy_fund %>% filter(Sig != "   "), row.names = FALSE)


# ── A.3  VAR: DXY in the funding-credit system ──────────────────────────────

cat("\n================================================================\n")
cat("  A.3 — VAR: DXY IN THE FUNDING-CREDIT SYSTEM\n")
cat("================================================================\n\n")

# Ordering: r_DXY (global, most exogenous) → d_TGA → d_Reserves → SOFR_EFFR → d_Baa
var_cols_g <- c("r_DXY", "d_TGA", "d_Reserves", "SOFR_EFFR", "d_Baa")
var_data_g <- df_post %>% dplyr::select(all_of(var_cols_g)) %>% drop_na()
cat(sprintf("  VAR G observations: %d\n", nrow(var_data_g)))

lag_sel_g <- VARselect(var_data_g, lag.max = 20, type = "const")
cat("  Lag selection:\n"); print(lag_sel_g$selection)
opt_lag_g <- lag_sel_g$selection["AIC(n)"]
cat(sprintf("  Selected lag: %d (AIC)\n\n", opt_lag_g))

var_g <- VAR(var_data_g, p = opt_lag_g, type = "const")
roots_g <- roots(var_g)
cat(sprintf("  Max eigenvalue: %.4f (Stable: %s)\n",
            max(roots_g), ifelse(max(roots_g) < 1, "YES", "NO")))

cat("\n  R-squared (adj):\n")
for (eq in names(var_g$varresult))
  cat(sprintf("    %s: %.4f\n", eq, summary(var_g$varresult[[eq]])$adj.r.squared))

cat("\n  Granger causality in VAR:\n")
for (cv in c("r_DXY", "d_TGA", "d_Reserves", "SOFR_EFFR")) {
  gc <- causality(var_g, cause = cv)
  cat(sprintf("    %s → others: F = %.3f, p = %.4f %s\n",
              cv, gc$Granger$statistic, gc$Granger$p.value,
              ifelse(gc$Granger$p.value < 0.05, "**", "")))
}

# Portmanteau
pt_g <- serial.test(var_g, lags.pt = 20, type = "PT.asymptotic")
cat(sprintf("\n  Portmanteau (20 lags): Chi-sq = %.1f, p = %.4f\n",
            pt_g$serial$statistic, pt_g$serial$p.value))

# IRFs
irf_g <- irf(var_g, impulse = NULL, response = NULL,
              n.ahead = n_ahead, boot = TRUE, runs = n_boot,
              ci = 0.95, ortho = TRUE)

panel_g_df <- rbind(
  make_irf_df(irf_g, "r_DXY",      "d_Baa",      "DXY → ΔBaa"),
  make_irf_df(irf_g, "r_DXY",      "SOFR_EFFR",  "DXY → SOFR-EFFR"),
  make_irf_df(irf_g, "r_DXY",      "d_Reserves",  "DXY → ΔReserves"),
  make_irf_df(irf_g, "d_Reserves", "d_Baa",       "ΔReserves → ΔBaa"),
  make_irf_df(irf_g, "d_Reserves", "r_DXY",       "ΔReserves → DXY"),
  make_irf_df(irf_g, "d_Baa",      "r_DXY",       "ΔBaa → DXY")
)

p_irf_dxy <- ggplot(panel_g_df, aes(x = Horizon)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.15) +
  geom_line(aes(y = IRF), color = "steelblue", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  facet_wrap(~ Panel, scales = "free_y", ncol = 3) +
  labs(title = "Impulse Responses — Dollar Index in the Funding-Credit System (Post-2018)",
       subtitle = "Cholesky: DXY → ΔTGA → ΔReserves → SOFR-EFFR → ΔBaa",
       x = "Horizon (business days)", y = "Response") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        strip.text = element_text(face = "bold", size = 8),
        panel.grid.minor = element_blank())

ggsave("fig_irf_dxy_system.pdf", p_irf_dxy, width = 9, height = 6, device = "pdf")
cat("  Saved: fig_irf_dxy_system.pdf\n")

# FEVD
fevd_g <- fevd(var_g, n.ahead = 30)

cat("\n  FEVD for d_Baa (VAR G):\n")
fm_baa <- fevd_g[["d_Baa"]]
for (h in c(1, 5, 10, 20, 30))
  cat(sprintf("    h=%2d: DXY=%.1f%%  TGA=%.1f%%  Res=%.1f%%  SOFR=%.1f%%  Own=%.1f%%\n",
              h, 100*fm_baa[h,"r_DXY"], 100*fm_baa[h,"d_TGA"],
              100*fm_baa[h,"d_Reserves"], 100*fm_baa[h,"SOFR_EFFR"],
              100*fm_baa[h,"d_Baa"]))

cat("\n  FEVD for r_DXY (VAR G):\n")
fm_dxy <- fevd_g[["r_DXY"]]
for (h in c(1, 5, 10, 30))
  cat(sprintf("    h=%2d: Own=%.1f%%  TGA=%.1f%%  Res=%.1f%%  SOFR=%.1f%%  Baa=%.1f%%\n",
              h, 100*fm_dxy[h,"r_DXY"], 100*fm_dxy[h,"d_TGA"],
              100*fm_dxy[h,"d_Reserves"], 100*fm_dxy[h,"SOFR_EFFR"],
              100*fm_dxy[h,"d_Baa"]))

# FEVD stacked area plot for d_Baa
fevd_g_plot <- data.frame()
for (h in 1:30) {
  for (src in var_cols_g)
    fevd_g_plot <- rbind(fevd_g_plot, data.frame(
      Horizon = h, Source = src,
      Share = 100 * fm_baa[h, src], stringsAsFactors = FALSE))
}
fevd_g_plot$Source <- factor(fevd_g_plot$Source, levels = var_cols_g,
  labels = c("DXY", "TGA", "Reserves", "SOFR-EFFR", "Own (Baa)"))

p_fevd_dxy <- ggplot(fevd_g_plot, aes(x = Horizon, y = Share, fill = Source)) +
  geom_area(alpha = 0.8) +
  scale_fill_brewer(palette = "Set2", name = "Shock source") +
  labs(title = "FEVD of ΔBaa: Dollar Shocks Dominate All Funding Variables",
       x = "Horizon (business days)", y = "Share of forecast error variance (%)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))

ggsave("fig_fevd_dxy_baa.pdf", p_fevd_dxy, width = 8, height = 4, device = "pdf")
cat("  Saved: fig_fevd_dxy_baa.pdf\n")


# ── A.4  VAR: DXY → ETF ratios ──────────────────────────────────────────────

cat("\n================================================================\n")
cat("  A.4 — VAR: DXY → ETF CREDIT RATIOS\n")
cat("================================================================\n\n")

# DXY, Reserves, HYG/SHV, HYG/LQD
var_cols_h <- c("r_DXY", "d_Reserves", "r_HYG_SHV", "r_HYG_LQD")
var_data_h <- df_post %>% dplyr::select(all_of(var_cols_h)) %>% drop_na()
cat(sprintf("  VAR H observations: %d\n", nrow(var_data_h)))

lag_sel_h <- VARselect(var_data_h, lag.max = 20, type = "const")
opt_lag_h <- lag_sel_h$selection["AIC(n)"]
cat(sprintf("  AIC lag: %d\n", opt_lag_h))

var_h <- VAR(var_data_h, p = opt_lag_h, type = "const")
cat(sprintf("  Max eigenvalue: %.4f (Stable: %s)\n",
            max(roots(var_h)), ifelse(max(roots(var_h)) < 1, "YES", "NO")))

cat("\n  R-squared:\n")
for (eq in names(var_h$varresult))
  cat(sprintf("    %s: %.4f\n", eq, summary(var_h$varresult[[eq]])$adj.r.squared))

cat("\n  Granger in VAR:\n")
for (cv in c("r_DXY", "d_Reserves")) {
  gc <- causality(var_h, cause = cv)
  cat(sprintf("    %s → others: F = %.3f, p = %.4f %s\n",
              cv, gc$Granger$statistic, gc$Granger$p.value,
              ifelse(gc$Granger$p.value < 0.05, "**", "")))
}

# FEVD
fevd_h <- fevd(var_h, n.ahead = 30)
cat("\n  FEVD for r_HYG_SHV:\n")
fh <- fevd_h[["r_HYG_SHV"]]
for (h in c(1, 5, 10, 30))
  cat(sprintf("    h=%2d: DXY=%.1f%%  Res=%.1f%%  Own=%.1f%%  HYG/LQD=%.1f%%\n",
              h, 100*fh[h,"r_DXY"], 100*fh[h,"d_Reserves"],
              100*fh[h,"r_HYG_SHV"], 100*fh[h,"r_HYG_LQD"]))

cat("\n  FEVD for r_HYG_LQD:\n")
fh2 <- fevd_h[["r_HYG_LQD"]]
for (h in c(1, 5, 10, 30))
  cat(sprintf("    h=%2d: DXY=%.1f%%  Res=%.1f%%  HYG/SHV=%.1f%%  Own=%.1f%%\n",
              h, 100*fh2[h,"r_DXY"], 100*fh2[h,"d_Reserves"],
              100*fh2[h,"r_HYG_SHV"], 100*fh2[h,"r_HYG_LQD"]))

# IRFs
irf_h <- irf(var_h, impulse = NULL, response = NULL,
              n.ahead = n_ahead, boot = TRUE, runs = n_boot,
              ci = 0.95, ortho = TRUE)

panel_h_df <- rbind(
  make_irf_df(irf_h, "r_DXY",      "r_HYG_SHV", "DXY → HYG/SHV"),
  make_irf_df(irf_h, "r_DXY",      "r_HYG_LQD", "DXY → HYG/LQD"),
  make_irf_df(irf_h, "d_Reserves", "r_HYG_SHV", "ΔReserves → HYG/SHV"),
  make_irf_df(irf_h, "d_Reserves", "r_HYG_LQD", "ΔReserves → HYG/LQD")
)

p_irf_etf <- ggplot(panel_h_df, aes(x = Horizon)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange3", alpha = 0.15) +
  geom_line(aes(y = IRF), color = "darkorange3", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  facet_wrap(~ Panel, scales = "free_y", ncol = 2) +
  labs(title = "Impulse Responses: DXY and Reserves → ETF Credit Ratios (Post-2018)",
       x = "Horizon (business days)", y = "Response") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        strip.text = element_text(face = "bold", size = 8),
        panel.grid.minor = element_blank())

ggsave("fig_irf_dxy_etf_ratios.pdf", p_irf_etf, width = 8, height = 5, device = "pdf")
cat("  Saved: fig_irf_dxy_etf_ratios.pdf\n")


###############################################################################
#  PART B — ETF ROBUSTNESS (brief)
###############################################################################

cat("\n\n================================================================\n")
cat("  PART B: ETF ROBUSTNESS — MOODY'S vs TRADED CREDIT\n")
cat("================================================================\n\n")

# Quick Granger: funding → ETF credit (confirming Moody's results)
granger_etf <- data.frame()
for (xv in c("SOFR_EFFR", "d_Reserves", "d_TGA")) {
  for (yv in c("r_HYG", "r_LQD", "d_log_HYG_IEF", "d_log_LQD_IEF")) {
    granger_etf <- rbind(granger_etf,
                         run_granger(df_post, xv, yv, lags = c(1, 5)))
  }
}

cat("  Significant results (funding → ETF credit):\n")
print(granger_etf %>% filter(Sig != "   "), row.names = FALSE)
cat(sprintf("\n  Significant at 10%%: %d / %d tests\n",
            sum(granger_etf$Sig != "   "), nrow(granger_etf)))


###############################################################################
#  PART C — CROSS-CURRENCY BASIS (JPY/USD)
###############################################################################

cat("\n\n================================================================\n")
cat("  PART C: CROSS-CURRENCY BASIS (JPY/USD)\n")
cat("================================================================\n\n")

ccb <- read_csv(file.path(WD, "..", "cross_currency_basis.csv"), show_col_types = FALSE)
ccb$Date <- as.Date(ccb$Date)

df_ccb <- df %>%
  inner_join(ccb %>% dplyr::select(Date, JPY), by = "Date") %>%
  arrange(Date) %>%
  mutate(
    d_JPY_basis = c(NA, diff(JPY)),
    r_DXY = c(NA, 100 * diff(log(DXY)))
  )

cat(sprintf("  Merged sample: %d rows, %s to %s\n",
            nrow(df_ccb), min(df_ccb$Date), max(df_ccb$Date)))
cat(sprintf("  JPY basis: mean = %.1f bp, sd = %.1f bp\n",
            mean(df_ccb$JPY, na.rm = TRUE), sd(df_ccb$JPY, na.rm = TRUE)))

# Key correlations
cat("\n  Correlations:\n")
for (pr in list(
  c("JPY", "SOFR_EFFR"), c("JPY", "Aaa_spread"), c("JPY", "DXY"),
  c("d_JPY_basis", "d_Baa"), c("d_JPY_basis", "d_Aaa"),
  c("d_JPY_basis", "r_DXY")
)) {
  sub <- df_ccb %>% dplyr::select(all_of(pr)) %>% drop_na()
  if (nrow(sub) < 30) next
  ct <- cor.test(sub[[1]], sub[[2]])
  cat(sprintf("    %-14s vs %-14s r = %7.4f (p = %.4f) %s\n",
              pr[1], pr[2], ct$estimate, ct$p.value,
              ifelse(ct$p.value < 0.05, "**", "")))
}

# Granger tests (both directions)
cat("\n  Granger causality:\n")
ccb_granger <- data.frame()

for (xv in c("SOFR_EFFR", "d_Reserves", "r_DXY")) {
  ccb_granger <- rbind(ccb_granger,
    run_granger(df_ccb, xv, "d_JPY_basis", lags = c(1, 5), min_obs = 50))
}
for (yv in c("d_Baa", "d_Aaa", "d_Baa_Aaa")) {
  ccb_granger <- rbind(ccb_granger,
    run_granger(df_ccb, "d_JPY_basis", yv, lags = c(1, 5), min_obs = 50))
}
for (xv in c("d_Baa", "d_Aaa")) {
  ccb_granger <- rbind(ccb_granger,
    run_granger(df_ccb, xv, "d_JPY_basis", lags = c(1, 5), min_obs = 50))
}
# DXY ↔ basis
ccb_granger <- rbind(ccb_granger,
  run_granger(df_ccb, "r_DXY", "d_JPY_basis", lags = c(1, 5), min_obs = 50))
ccb_granger <- rbind(ccb_granger,
  run_granger(df_ccb, "d_JPY_basis", "r_DXY", lags = c(1, 5), min_obs = 50))

print(ccb_granger, row.names = FALSE)

cat("\n  Significant results:\n")
sig_ccb <- ccb_granger %>% filter(Sig != "   ")
if (nrow(sig_ccb) > 0) print(sig_ccb, row.names = FALSE) else cat("    None\n")

# Time series plot
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

# Rolling correlation: DXY vs JPY basis
df_ccb_roll <- df_ccb %>%
  filter(!is.na(r_DXY), !is.na(d_JPY_basis)) %>%
  arrange(Date) %>%
  mutate(
    roll_dxy_basis = zoo::rollapply(
      cbind(r_DXY, d_JPY_basis), width = 60,
      FUN = function(x) cor(x[,1], x[,2], use = "complete.obs"),
      by.column = FALSE, fill = NA, align = "right"))

p_roll <- ggplot(df_ccb_roll %>% filter(!is.na(roll_dxy_basis)),
                 aes(x = Date, y = roll_dxy_basis)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  geom_line(color = "darkred", linewidth = 0.6) +
  labs(title = "60-Day Rolling Correlation: DXY Returns vs Δ(JPY Basis)",
       x = NULL, y = "Correlation") +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_blank())

ggsave("fig_rolling_corr_dxy_basis.pdf", p_roll, width = 8, height = 3.5, device = "pdf")
cat("  Saved: fig_rolling_corr_dxy_basis.pdf\n")


###############################################################################
#  LATEX TABLES
###############################################################################

cat("\n================================================================\n")
cat("  GENERATING LATEX TABLES\n")
cat("================================================================\n\n")

# Table: DXY Granger results
dxy_all <- rbind(
  granger_dxy %>% mutate(Battery = "DXY_to_credit"),
  gc_rev %>% mutate(Battery = "credit_to_DXY"),
  gc_fund_dxy %>% mutate(Battery = "funding_to_DXY"),
  gc_dxy_fund %>% mutate(Battery = "DXY_to_funding")
)

dxy_xt <- xtable(dxy_all %>% dplyr::select(-Battery),
  caption = "Granger Causality: Dollar Index and Credit/Funding Variables (Post-2018)",
  label = "tab:granger_dxy",
  digits = c(0, 0, 0, 0, 0, 3, 4, 0))
print(dxy_xt, file = "table_granger_dxy.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp")
cat("  Saved: table_granger_dxy.tex\n")

# Table: FEVD comparison
fevd_comp <- data.frame(
  Shock = c("DXY", "Reserves", "TGA", "SOFR-EFFR", "Own (Baa)"),
  h1  = sprintf("%.1f", 100*fm_baa[1,]),
  h5  = sprintf("%.1f", 100*fm_baa[5,]),
  h10 = sprintf("%.1f", 100*fm_baa[10,]),
  h30 = sprintf("%.1f", 100*fm_baa[30,]),
  stringsAsFactors = FALSE
)
fevd_xt <- xtable(fevd_comp,
  caption = "FEVD of $\\Delta$Baa: Dollar Index vs Funding Variables (VAR G, Post-2018)",
  label = "tab:fevd_dxy",
  digits = 0)
print(fevd_xt, file = "table_fevd_dxy.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.colnames.function = function(x) gsub("h", "$h=$", x))
cat("  Saved: table_fevd_dxy.tex\n")

# Table: CCB Granger
ccb_xt <- xtable(ccb_granger,
  caption = "Granger Causality: Cross-Currency Basis Sample (May 2024 -- May 2026)",
  label = "tab:granger_ccb",
  digits = c(0, 0, 0, 0, 0, 3, 4, 0))
print(ccb_xt, file = "table_granger_ccb.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp")
cat("  Saved: table_granger_ccb.tex\n")

# Save all results
write_csv(rbind(
  dxy_all,
  granger_etf %>% mutate(Battery = "ETF_robustness"),
  ccb_granger %>% mutate(Battery = "CCB")
), "granger_extension_results.csv")
cat("  Saved: granger_extension_results.csv\n")


cat("\n================================================================\n")
cat("  STEP 5 COMPLETE\n")
cat("================================================================\n")
