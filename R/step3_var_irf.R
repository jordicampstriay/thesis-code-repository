###############################################################################
#  STEP 3 — VAR ESTIMATION AND IMPULSE RESPONSE FUNCTIONS  (Section 5.4)
#  Thesis: "The Plumbing of Liquidity"
#  Author: Jordi Camps Triay
#  Last updated: 2026-05-12
#
#  Prerequisites: data_master.csv from Step 1, granger_all_results.csv from Step 2
#
#  DESIGN:
#    The Granger tests (Step 2) established that:
#      (a) d_Reserves → d_Aaa at all lags (collateral channel)
#      (b) d_TGA → d_Baa_Aaa at lags 5,10 (fiscal channel)
#      (c) SOFR-EFFR → nothing post-2018 (nonlinearity, motivates Step 4)
#      (d) No reverse causality credit → SOFR-EFFR (exogeneity confirmed)
#      (e) Bidirectional causality credit ↔ d_Reserves at longer lags
#      (f) ON_RRP insignificant everywhere (excluded)
#
#    We estimate a VAR to:
#      1. Formally select the lag order via information criteria (AIC, BIC, HQ)
#      2. Trace impulse response functions (IRFs) showing the dynamic path
#         of credit spread responses to funding shocks
#      3. Compute forecast error variance decomposition (FEVD) to quantify
#         how much credit spread variation is attributable to each channel
#      4. Provide structural identification via Cholesky decomposition
#
#    Cholesky ordering (most exogenous → most endogenous):
#      d_TGA → d_Reserves → SOFR_EFFR → d_Aaa → d_Baa
#
#    Rationale: TGA is driven by fiscal cash flows (exogenous). Reserves
#    respond to TGA and monetary policy. SOFR-EFFR is market-determined
#    but exogenous to credit (no reverse causality). Aaa is ordered before
#    Baa because it responds first to collateral shocks (closer Treasury
#    substitute). Baa is last (most endogenous, responds to all channels).
#
#  OUTPUT: LaTeX tables, PDF figures, var_results.RData
###############################################################################

# ── 0. PACKAGES ──────────────────────────────────────────────────────────────
required_pkgs <- c("readr", "dplyr", "tidyr", "vars", "xtable", "ggplot2",
                   "scales", "gridExtra", "zoo")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

source(file.path(dirname(if (interactive()) rstudioapi::getSourceEditorContext()$path else sys.frame(1)$ofile), "config.R"))

cat("================================================================\n")
cat("  STEP 3 — VAR ESTIMATION AND IMPULSE RESPONSE FUNCTIONS\n")
cat("================================================================\n\n")


# ── 1. LOAD AND PREPARE DATA ────────────────────────────────────────────────
df <- read_csv(file.path(DATA_DIR, "data_master.csv"), show_col_types = FALSE)
df$Date <- as.Date(df$Date)

# Post-2018 sample (primary)
df_post <- df %>% filter(Date >= as.Date("2018-04-02"))

# Extended sample (robustness)
df_extended <- df %>% filter(Date >= as.Date("2003-05-01"))

cat(sprintf("  Post-2018: %d rows | Extended: %d rows\n\n",
            nrow(df_post), nrow(df_extended)))


###############################################################################
#  VAR SPECIFICATION A — 5-VARIABLE MODEL (Post-2018)
#
#  Ordering: d_TGA, d_Reserves, SOFR_EFFR, d_Aaa, d_Baa
#
#  All variables are stationary (confirmed in Step 1 ADF tests):
#    - d_TGA, d_Reserves, d_Baa, d_Aaa: first differences of I(1) series
#    - SOFR_EFFR: stationary in levels (ADF p = 0.01)
###############################################################################

cat("================================================================\n")
cat("  VAR A: 5-Variable Model (Post-2018)\n")
cat("================================================================\n\n")

# Prepare the data matrix
var_cols_a <- c("d_TGA", "d_Reserves", "SOFR_EFFR", "d_Aaa", "d_Baa")
var_data_a <- df_post %>%
  dplyr::select(all_of(var_cols_a)) %>%
  drop_na()

cat(sprintf("  VAR A observations: %d (after dropping NAs)\n", nrow(var_data_a)))
cat(sprintf("  Variables: %s\n\n", paste(var_cols_a, collapse = ", ")))


# ── 1a. LAG SELECTION ────────────────────────────────────────────────────────
cat("--- Lag Selection ---\n")
max_lag <- 15  # Test up to 15 lags (3 trading weeks)
lag_sel_a <- VARselect(var_data_a, lag.max = max_lag, type = "const")
cat("  Information criteria:\n")
print(lag_sel_a$selection)
cat("\n")

# Extract optimal lags from each criterion
aic_lag <- lag_sel_a$selection["AIC(n)"]
bic_lag <- lag_sel_a$selection["SC(n)"]
hq_lag  <- lag_sel_a$selection["HQ(n)"]

# Use AIC lag for primary spec (tends to be more generous, captures dynamics)
# BIC often selects lag 1 which may miss the weekly transmission documented in Step 2
opt_lag_a <- aic_lag
cat(sprintf("  Selected lag: %d (AIC). BIC suggests: %d. HQ suggests: %d\n\n",
            opt_lag_a, bic_lag, hq_lag))


# ── 1b. VAR ESTIMATION ──────────────────────────────────────────────────────
cat("--- VAR Estimation ---\n")
var_a <- VAR(var_data_a, p = opt_lag_a, type = "const")

# Summary statistics
cat(sprintf("  Estimated VAR(%d) with %d variables\n", opt_lag_a, length(var_cols_a)))
cat(sprintf("  Total parameters per equation: %d\n", opt_lag_a * length(var_cols_a) + 1))
cat(sprintf("  Effective observations: %d\n\n", nrow(var_data_a) - opt_lag_a))

# Print adjusted R-squared for each equation
for (eq_name in names(var_a$varresult)) {
  r2 <- summary(var_a$varresult[[eq_name]])$adj.r.squared
  cat(sprintf("  R-squared (adj) for %s equation: %.4f\n", eq_name, r2))
}
cat("\n")


# ── 1c. STABILITY CHECK ─────────────────────────────────────────────────────
cat("--- Stability Check ---\n")
roots_a <- roots(var_a)
cat(sprintf("  Max eigenvalue modulus: %.4f\n", max(roots_a)))
if (max(roots_a) < 1) {
  cat("  VAR is STABLE (all roots inside unit circle)\n\n")
} else {
  cat("  WARNING: VAR is NOT stable\n\n")
}


# ── 1d. RESIDUAL DIAGNOSTICS ────────────────────────────────────────────────
cat("--- Residual Diagnostics ---\n")
# Portmanteau test for serial correlation
pt_test <- serial.test(var_a, lags.pt = 20, type = "PT.asymptotic")
cat(sprintf("  Portmanteau test (20 lags): Chi-sq = %.2f, df = %d, p = %.4f\n",
            pt_test$serial$statistic, pt_test$serial$parameter, pt_test$serial$p.value))
if (pt_test$serial$p.value > 0.05) {
  cat("  No serial correlation at 5% level (good)\n\n")
} else {
  cat("  Serial correlation detected at 5% level (residuals not white noise)\n\n")
}


# ── 1e. GRANGER CAUSALITY WITHIN VAR ────────────────────────────────────────
cat("--- Granger Causality in VAR Framework ---\n")
for (cause_var in c("d_TGA", "d_Reserves", "SOFR_EFFR")) {
  gc <- causality(var_a, cause = cause_var)
  cat(sprintf("  %s → others: F = %.3f, p = %.4f %s\n",
              cause_var,
              gc$Granger$statistic,
              gc$Granger$p.value,
              ifelse(gc$Granger$p.value < 0.05, "**", "")))
}
cat("\n")


###############################################################################
#  IMPULSE RESPONSE FUNCTIONS — VAR A (Post-2018)
###############################################################################

cat("================================================================\n")
cat("  IMPULSE RESPONSE FUNCTIONS — VAR A\n")
cat("================================================================\n\n")

n_ahead <- 30  # 30 business days (6 weeks) horizon
n_boot  <- 500 # Bootstrap replications for confidence intervals

# Compute structural (Cholesky) IRFs
irf_a <- irf(var_a, impulse = NULL, response = NULL,
              n.ahead = n_ahead, boot = TRUE, runs = n_boot,
              ci = 0.95, ortho = TRUE)

# ── Plot IRFs: Funding shocks → Credit responses ────────────────────────────

# Custom plotting function for publication-quality IRFs
plot_irf <- function(irf_obj, impulse, response, title, filename,
                     y_label = "Response (bp)", scale_factor = 1) {
  idx <- which(names(irf_obj$irf) == impulse)
  resp_idx <- which(colnames(irf_obj$irf[[impulse]]) == response)

  irf_vals  <- irf_obj$irf[[impulse]][, resp_idx] * scale_factor
  lower     <- irf_obj$Lower[[impulse]][, resp_idx] * scale_factor
  upper     <- irf_obj$Upper[[impulse]][, resp_idx] * scale_factor
  horizon   <- 0:n_ahead

  plot_df <- data.frame(
    Horizon = horizon,
    IRF = irf_vals,
    Lower = lower,
    Upper = upper
  )

  p <- ggplot(plot_df, aes(x = Horizon)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.2) +
    geom_line(aes(y = IRF), color = "steelblue", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    labs(title = title,
         x = "Horizon (business days)",
         y = y_label) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(FIG_DIR, filename), p, width = 6, height = 3.5, device = "pdf")
  cat(sprintf("  Saved: %s\n", filename))
  return(p)
}

# Key IRFs: Funding → Credit
p1 <- plot_irf(irf_a, "d_Reserves", "d_Aaa",
               "Response of d_Aaa to d_Reserves shock (Cholesky, Post-2018)",
               "fig_irf_reserves_aaa.pdf")

p2 <- plot_irf(irf_a, "d_Reserves", "d_Baa",
               "Response of d_Baa to d_Reserves shock (Cholesky, Post-2018)",
               "fig_irf_reserves_baa.pdf")

p3 <- plot_irf(irf_a, "d_TGA", "d_Aaa",
               "Response of d_Aaa to d_TGA shock (Cholesky, Post-2018)",
               "fig_irf_tga_aaa.pdf")

p4 <- plot_irf(irf_a, "d_TGA", "d_Baa",
               "Response of d_Baa to d_TGA shock (Cholesky, Post-2018)",
               "fig_irf_tga_baa.pdf")

p5 <- plot_irf(irf_a, "SOFR_EFFR", "d_Aaa",
               "Response of d_Aaa to SOFR-EFFR shock (Cholesky, Post-2018)",
               "fig_irf_sofreffr_aaa.pdf")

p6 <- plot_irf(irf_a, "SOFR_EFFR", "d_Baa",
               "Response of d_Baa to SOFR-EFFR shock (Cholesky, Post-2018)",
               "fig_irf_sofreffr_baa.pdf")

# Cross-channel: Reserves → SOFR-EFFR (does quantity affect price?)
p7 <- plot_irf(irf_a, "d_Reserves", "SOFR_EFFR",
               "Response of SOFR-EFFR to d_Reserves shock (Cholesky, Post-2018)",
               "fig_irf_reserves_sofreffr.pdf")

# TGA → Reserves (fiscal plumbing)
p8 <- plot_irf(irf_a, "d_TGA", "d_Reserves",
               "Response of d_Reserves to d_TGA shock (Cholesky, Post-2018)",
               "fig_irf_tga_reserves.pdf")

cat("\n")


# ── Combined panel figure ────────────────────────────────────────────────────
# Create a 2x3 panel of key IRFs for the report

make_irf_panel_df <- function(irf_obj, impulse, response, label) {
  idx <- which(names(irf_obj$irf) == impulse)
  resp_idx <- which(colnames(irf_obj$irf[[impulse]]) == response)
  data.frame(
    Horizon = 0:n_ahead,
    IRF = irf_obj$irf[[impulse]][, resp_idx],
    Lower = irf_obj$Lower[[impulse]][, resp_idx],
    Upper = irf_obj$Upper[[impulse]][, resp_idx],
    Panel = label,
    stringsAsFactors = FALSE
  )
}

panel_df <- rbind(
  make_irf_panel_df(irf_a, "d_Reserves", "d_Aaa",   "d_Reserves -> d_Aaa"),
  make_irf_panel_df(irf_a, "d_Reserves", "d_Baa",   "d_Reserves -> d_Baa"),
  make_irf_panel_df(irf_a, "d_TGA",      "d_Aaa",   "d_TGA -> d_Aaa"),
  make_irf_panel_df(irf_a, "d_TGA",      "d_Baa",   "d_TGA -> d_Baa"),
  make_irf_panel_df(irf_a, "SOFR_EFFR",  "d_Aaa",   "SOFR_EFFR -> d_Aaa"),
  make_irf_panel_df(irf_a, "SOFR_EFFR",  "d_Baa",   "SOFR_EFFR -> d_Baa")
)

# Order panels logically
panel_df$Panel <- factor(panel_df$Panel,
  levels = c("d_Reserves -> d_Aaa", "d_Reserves -> d_Baa",
             "d_TGA -> d_Aaa",      "d_TGA -> d_Baa",
             "SOFR_EFFR -> d_Aaa",  "SOFR_EFFR -> d_Baa"))

p_panel <- ggplot(panel_df, aes(x = Horizon)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = IRF), color = "steelblue", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  facet_wrap(~ Panel, scales = "free_y", ncol = 2) +
  labs(x = "Horizon (business days)", y = "Response") +
  theme_minimal(base_size = 9) +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 10, 5, 5)
  )

ggsave(file.path(FIG_DIR, "fig_irf_panel_post2018.pdf"), p_panel, width = 8, height = 7, device = "pdf")
cat("  Saved: fig_irf_panel_post2018.pdf\n\n")


###############################################################################
#  FORECAST ERROR VARIANCE DECOMPOSITION — VAR A
###############################################################################

cat("================================================================\n")
cat("  FORECAST ERROR VARIANCE DECOMPOSITION — VAR A\n")
cat("================================================================\n\n")

fevd_a <- fevd(var_a, n.ahead = 30)

# Print FEVD at horizons 1, 5, 10, 20, 30
for (resp in c("d_Aaa", "d_Baa")) {
  cat(sprintf("  FEVD for %s:\n", resp))
  fevd_mat <- fevd_a[[resp]]
  for (h in c(1, 5, 10, 20, 30)) {
    if (h <= nrow(fevd_mat)) {
      cat(sprintf("    h=%2d: d_TGA=%.1f%% d_Res=%.1f%% SOFR=%.1f%% d_Aaa=%.1f%% d_Baa=%.1f%%\n",
                  h,
                  100*fevd_mat[h, "d_TGA"],
                  100*fevd_mat[h, "d_Reserves"],
                  100*fevd_mat[h, "SOFR_EFFR"],
                  100*fevd_mat[h, "d_Aaa"],
                  100*fevd_mat[h, "d_Baa"]))
    }
  }
  cat("\n")
}

# FEVD plot
fevd_plot_df <- data.frame()
for (resp in c("d_Aaa", "d_Baa")) {
  fevd_mat <- fevd_a[[resp]]
  for (h in 1:nrow(fevd_mat)) {
    for (src in var_cols_a) {
      fevd_plot_df <- rbind(fevd_plot_df, data.frame(
        Response = resp,
        Horizon = h,
        Source = src,
        Share = 100 * fevd_mat[h, src],
        stringsAsFactors = FALSE
      ))
    }
  }
}

fevd_plot_df$Source <- factor(fevd_plot_df$Source,
  levels = var_cols_a,
  labels = c("TGA", "Reserves", "SOFR-EFFR", "Aaa", "Baa"))

fevd_plot_df$Response <- factor(fevd_plot_df$Response,
  levels = c("d_Aaa", "d_Baa"),
  labels = c("Variance of d_Aaa", "Variance of d_Baa"))

p_fevd <- ggplot(fevd_plot_df, aes(x = Horizon, y = Share, fill = Source)) +
  geom_area(alpha = 0.8) +
  facet_wrap(~ Response, ncol = 2) +
  scale_fill_brewer(palette = "Set2", name = "Shock source") +
  labs(x = "Horizon (business days)", y = "Share of forecast error variance (%)") +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIG_DIR, "fig_fevd_post2018.pdf"), p_fevd, width = 8, height = 4, device = "pdf")
cat("  Saved: fig_fevd_post2018.pdf\n\n")


###############################################################################
#  VAR SPECIFICATION B — QUALITY DIFFERENTIAL (Post-2018)
#
#  Ordering: d_TGA, d_Reserves, SOFR_EFFR, d_Baa_Aaa
#  Focuses on the Baa-Aaa quality differential to isolate the collateral
#  channel (which affects Aaa more than Baa, compressing the differential)
###############################################################################

cat("================================================================\n")
cat("  VAR B: Quality Differential Model (Post-2018)\n")
cat("================================================================\n\n")

var_cols_b <- c("d_TGA", "d_Reserves", "SOFR_EFFR", "d_Baa_Aaa")
var_data_b <- df_post %>%
  dplyr::select(all_of(var_cols_b)) %>%
  drop_na()

cat(sprintf("  VAR B observations: %d\n", nrow(var_data_b)))

# Lag selection
lag_sel_b <- VARselect(var_data_b, lag.max = max_lag, type = "const")
cat("  Lag selection:\n")
print(lag_sel_b$selection)
opt_lag_b <- lag_sel_b$selection["AIC(n)"]
cat(sprintf("  Selected lag: %d (AIC)\n\n", opt_lag_b))

# Estimate
var_b <- VAR(var_data_b, p = opt_lag_b, type = "const")

# Stability
roots_b <- roots(var_b)
cat(sprintf("  Max eigenvalue modulus: %.4f (Stable: %s)\n\n",
            max(roots_b), ifelse(max(roots_b) < 1, "YES", "NO")))

# R-squared
for (eq_name in names(var_b$varresult)) {
  r2 <- summary(var_b$varresult[[eq_name]])$adj.r.squared
  cat(sprintf("  R-squared (adj) for %s: %.4f\n", eq_name, r2))
}
cat("\n")

# IRFs for quality differential
irf_b <- irf(var_b, impulse = NULL, response = NULL,
              n.ahead = n_ahead, boot = TRUE, runs = n_boot,
              ci = 0.95, ortho = TRUE)

# Panel: all shocks → d_Baa_Aaa
panel_b_df <- rbind(
  make_irf_panel_df(irf_b, "d_Reserves", "d_Baa_Aaa", "d_Reserves -> d_Baa_Aaa"),
  make_irf_panel_df(irf_b, "d_TGA",      "d_Baa_Aaa", "d_TGA -> d_Baa_Aaa"),
  make_irf_panel_df(irf_b, "SOFR_EFFR",  "d_Baa_Aaa", "SOFR_EFFR -> d_Baa_Aaa")
)

panel_b_df$Panel <- factor(panel_b_df$Panel,
  levels = c("d_Reserves -> d_Baa_Aaa", "d_TGA -> d_Baa_Aaa", "SOFR_EFFR -> d_Baa_Aaa"))

p_panel_b <- ggplot(panel_b_df, aes(x = Horizon)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange", alpha = 0.2) +
  geom_line(aes(y = IRF), color = "darkorange3", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  facet_wrap(~ Panel, scales = "free_y", ncol = 3) +
  labs(x = "Horizon (business days)", y = "Response of d(Baa-Aaa)") +
  theme_minimal(base_size = 9) +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIG_DIR, "fig_irf_panel_quality_diff.pdf"), p_panel_b, width = 9, height = 3.5, device = "pdf")
cat("  Saved: fig_irf_panel_quality_diff.pdf\n\n")

# FEVD for quality differential
fevd_b <- fevd(var_b, n.ahead = 30)
cat("  FEVD for d_Baa_Aaa:\n")
fevd_baa_aaa <- fevd_b[["d_Baa_Aaa"]]
for (h in c(1, 5, 10, 20, 30)) {
  if (h <= nrow(fevd_baa_aaa)) {
    cat(sprintf("    h=%2d: d_TGA=%.1f%% d_Res=%.1f%% SOFR=%.1f%% d_Baa_Aaa=%.1f%%\n",
                h,
                100*fevd_baa_aaa[h, "d_TGA"],
                100*fevd_baa_aaa[h, "d_Reserves"],
                100*fevd_baa_aaa[h, "SOFR_EFFR"],
                100*fevd_baa_aaa[h, "d_Baa_Aaa"]))
  }
}
cat("\n")


###############################################################################
#  VAR C — EXTENDED SAMPLE (2003-2026) — ROBUSTNESS
###############################################################################

cat("================================================================\n")
cat("  VAR C: Extended Sample (2003-2026)\n")
cat("================================================================\n\n")

var_cols_c <- c("d_TGA", "d_Reserves", "SOFR_EFFR", "d_Aaa", "d_Baa")
var_data_c <- df_extended %>%
  dplyr::select(all_of(var_cols_c)) %>%
  drop_na()

cat(sprintf("  VAR C observations: %d\n", nrow(var_data_c)))

# Lag selection
lag_sel_c <- VARselect(var_data_c, lag.max = max_lag, type = "const")
cat("  Lag selection:\n")
print(lag_sel_c$selection)
opt_lag_c <- lag_sel_c$selection["AIC(n)"]
cat(sprintf("  Selected lag: %d (AIC)\n\n", opt_lag_c))

# Estimate
var_c <- VAR(var_data_c, p = opt_lag_c, type = "const")

# Stability
roots_c <- roots(var_c)
cat(sprintf("  Max eigenvalue modulus: %.4f (Stable: %s)\n\n",
            max(roots_c), ifelse(max(roots_c) < 1, "YES", "NO")))

# R-squared
for (eq_name in names(var_c$varresult)) {
  r2 <- summary(var_c$varresult[[eq_name]])$adj.r.squared
  cat(sprintf("  R-squared (adj) for %s: %.4f\n", eq_name, r2))
}
cat("\n")

# IRFs extended
irf_c <- irf(var_c, impulse = NULL, response = NULL,
              n.ahead = n_ahead, boot = TRUE, runs = n_boot,
              ci = 0.95, ortho = TRUE)

# Panel for extended sample
panel_c_df <- rbind(
  make_irf_panel_df(irf_c, "d_Reserves", "d_Aaa",   "d_Reserves -> d_Aaa"),
  make_irf_panel_df(irf_c, "d_Reserves", "d_Baa",   "d_Reserves -> d_Baa"),
  make_irf_panel_df(irf_c, "d_TGA",      "d_Aaa",   "d_TGA -> d_Aaa"),
  make_irf_panel_df(irf_c, "d_TGA",      "d_Baa",   "d_TGA -> d_Baa"),
  make_irf_panel_df(irf_c, "SOFR_EFFR",  "d_Aaa",   "SOFR_EFFR -> d_Aaa"),
  make_irf_panel_df(irf_c, "SOFR_EFFR",  "d_Baa",   "SOFR_EFFR -> d_Baa")
)

panel_c_df$Panel <- factor(panel_c_df$Panel,
  levels = c("d_Reserves -> d_Aaa", "d_Reserves -> d_Baa",
             "d_TGA -> d_Aaa",      "d_TGA -> d_Baa",
             "SOFR_EFFR -> d_Aaa",  "SOFR_EFFR -> d_Baa"))

p_panel_c <- ggplot(panel_c_df, aes(x = Horizon)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkred", alpha = 0.15) +
  geom_line(aes(y = IRF), color = "darkred", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  facet_wrap(~ Panel, scales = "free_y", ncol = 2) +
  labs(title = "Impulse Responses — Extended Sample (2003-2026, SOFR proxy)",
       x = "Horizon (business days)", y = "Response") +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIG_DIR, "fig_irf_panel_extended.pdf"), p_panel_c, width = 8, height = 7, device = "pdf")
cat("  Saved: fig_irf_panel_extended.pdf\n\n")

# FEVD extended
fevd_c <- fevd(var_c, n.ahead = 30)
for (resp in c("d_Aaa", "d_Baa")) {
  cat(sprintf("  FEVD for %s (extended):\n", resp))
  fevd_mat <- fevd_c[[resp]]
  for (h in c(1, 5, 10, 20, 30)) {
    if (h <= nrow(fevd_mat)) {
      cat(sprintf("    h=%2d: d_TGA=%.1f%% d_Res=%.1f%% SOFR=%.1f%% d_Aaa=%.1f%% d_Baa=%.1f%%\n",
                  h,
                  100*fevd_mat[h, "d_TGA"],
                  100*fevd_mat[h, "d_Reserves"],
                  100*fevd_mat[h, "SOFR_EFFR"],
                  100*fevd_mat[h, "d_Aaa"],
                  100*fevd_mat[h, "d_Baa"]))
    }
  }
  cat("\n")
}


###############################################################################
#  CHOLESKY ORDERING ROBUSTNESS
#
#  Test alternative ordering: d_Reserves, d_TGA, SOFR_EFFR, d_Aaa, d_Baa
#  (Reserves before TGA — tests whether ordering drives the results)
###############################################################################

cat("================================================================\n")
cat("  CHOLESKY ORDERING ROBUSTNESS\n")
cat("================================================================\n\n")

var_cols_alt <- c("d_Reserves", "d_TGA", "SOFR_EFFR", "d_Aaa", "d_Baa")
var_data_alt <- df_post %>%
  dplyr::select(all_of(var_cols_alt)) %>%
  drop_na()

var_alt <- VAR(var_data_alt, p = opt_lag_a, type = "const")
irf_alt <- irf(var_alt, impulse = NULL, response = NULL,
                n.ahead = n_ahead, boot = TRUE, runs = n_boot,
                ci = 0.95, ortho = TRUE)

# Compare key IRFs: Reserves → d_Aaa under both orderings
compare_df <- rbind(
  cbind(make_irf_panel_df(irf_a, "d_Reserves", "d_Aaa", "Primary ordering"),
        Ordering = "TGA-first"),
  cbind(make_irf_panel_df(irf_alt, "d_Reserves", "d_Aaa", "Alt ordering"),
        Ordering = "Reserves-first")
)

p_robust_order <- ggplot(compare_df, aes(x = Horizon, color = Ordering, fill = Ordering)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.1, color = NA) +
  geom_line(aes(y = IRF), linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("steelblue", "darkorange3")) +
  scale_fill_manual(values = c("steelblue", "darkorange3")) +
  labs(title = "Robustness: Response of d_Aaa to d_Reserves shock under alternative Cholesky orderings",
       x = "Horizon (business days)", y = "Response") +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 9, face = "bold")
  )

ggsave(file.path(FIG_DIR, "fig_irf_ordering_robustness.pdf"), p_robust_order, width = 7, height = 4, device = "pdf")
cat("  Saved: fig_irf_ordering_robustness.pdf\n\n")


###############################################################################
#  LATEX TABLES
###############################################################################

cat("================================================================\n")
cat("  GENERATING LATEX TABLES\n")
cat("================================================================\n\n")

# ── Table: Lag Selection ──
lag_tbl <- data.frame(
  Criterion = c("AIC", "SC (BIC)", "HQ", "FPE"),
  `VAR_A` = as.integer(lag_sel_a$selection),
  `VAR_B` = as.integer(lag_sel_b$selection),
  `VAR_C` = as.integer(lag_sel_c$selection),
  stringsAsFactors = FALSE
)
colnames(lag_tbl) <- c("Criterion", "VAR A (Post-2018)", "VAR B (Diff, Post-2018)",
                        "VAR C (Extended)")

lag_xt <- xtable(lag_tbl,
  caption = "VAR Lag Order Selection by Information Criteria",
  label = "tab:var_lag_selection",
  digits = 0)
print(lag_xt, file = file.path(TBL_DIR, "table_var_lag_selection.tex"),
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp")
cat("  Saved: table_var_lag_selection.tex\n")


# ── Table: VAR Diagnostics ──
diag_tbl <- data.frame(
  Model = c("VAR A: 5-var (Post-2018)", "VAR B: 4-var Diff (Post-2018)",
            "VAR C: 5-var (Extended)"),
  Lags = c(opt_lag_a, opt_lag_b, opt_lag_c),
  Obs = c(nrow(var_data_a) - opt_lag_a,
          nrow(var_data_b) - opt_lag_b,
          nrow(var_data_c) - opt_lag_c),
  MaxRoot = c(max(roots_a), max(roots_b), max(roots_c)),
  Stable = c(max(roots_a) < 1, max(roots_b) < 1, max(roots_c) < 1),
  stringsAsFactors = FALSE
)
diag_tbl$Stable <- ifelse(diag_tbl$Stable, "Yes", "No")
diag_tbl$MaxRoot <- round(diag_tbl$MaxRoot, 4)

diag_xt <- xtable(diag_tbl,
  caption = "VAR Model Diagnostics",
  label = "tab:var_diagnostics",
  digits = c(0, 0, 0, 0, 4, 0))
print(diag_xt, file = file.path(TBL_DIR, "table_var_diagnostics.tex"),
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp")
cat("  Saved: table_var_diagnostics.tex\n")


# ── Table: FEVD at selected horizons ──
fevd_summary <- data.frame()
for (h in c(1, 5, 10, 20, 30)) {
  for (resp in c("d_Aaa", "d_Baa")) {
    row <- data.frame(
      Response = resp,
      Horizon = h,
      d_TGA = round(100 * fevd_a[[resp]][h, "d_TGA"], 1),
      d_Reserves = round(100 * fevd_a[[resp]][h, "d_Reserves"], 1),
      SOFR_EFFR = round(100 * fevd_a[[resp]][h, "SOFR_EFFR"], 1),
      Own = round(100 * fevd_a[[resp]][h, resp], 1),
      stringsAsFactors = FALSE
    )
    fevd_summary <- rbind(fevd_summary, row)
  }
}
colnames(fevd_summary) <- c("Response", "$h$", "$\\Delta$TGA (\\%)",
                             "$\\Delta$Reserves (\\%)", "SOFR--EFFR (\\%)",
                             "Own (\\%)")

fevd_xt <- xtable(fevd_summary,
  caption = "Forecast Error Variance Decomposition (Post-2018, VAR A)",
  label = "tab:var_fevd",
  digits = c(0, 0, 0, 1, 1, 1, 1))
print(fevd_xt, file = file.path(TBL_DIR, "table_var_fevd.tex"),
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.colnames.function = identity)
cat("  Saved: table_var_fevd.tex\n")


# ── Table: R-squared for each VAR equation ──
rsq_tbl <- data.frame()
for (eq_name in names(var_a$varresult)) {
  rsq_tbl <- rbind(rsq_tbl, data.frame(
    Equation = eq_name,
    R2_adj_A = round(summary(var_a$varresult[[eq_name]])$adj.r.squared, 4),
    stringsAsFactors = FALSE
  ))
}
# Add VAR C
rsq_c <- data.frame()
for (eq_name in names(var_c$varresult)) {
  rsq_c <- rbind(rsq_c, data.frame(
    Equation = eq_name,
    R2_adj_C = round(summary(var_c$varresult[[eq_name]])$adj.r.squared, 4),
    stringsAsFactors = FALSE
  ))
}
rsq_tbl <- merge(rsq_tbl, rsq_c, by = "Equation", all = TRUE)
colnames(rsq_tbl) <- c("Equation", "Adj. $R^2$ (Post-2018)", "Adj. $R^2$ (Extended)")

rsq_xt <- xtable(rsq_tbl,
  caption = "Adjusted $R^2$ by VAR Equation",
  label = "tab:var_rsquared",
  digits = c(0, 0, 4, 4))
print(rsq_xt, file = file.path(TBL_DIR, "table_var_rsquared.tex"),
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity)
cat("  Saved: table_var_rsquared.tex\n\n")


###############################################################################
#  VAR D — TED SPREAD (2003-2022)
#
#  RATIONALE: LIBOR was the dominant benchmark interest rate before its
#  cessation in January 2022. The TED spread (3M LIBOR - 3M T-bill)
#  provides continuous variation in the funding price variable, unlike
#  SOFR-EFFR which is near-zero post-2018. The TED spread produced the
#  strongest Granger results in Step 2 (F=102.9 for Baa-Aaa at lag 1).
#
#  CAVEAT: TED contains bank credit risk via LIBOR, so its predictive
#  power for credit spreads is partly mechanical (shared credit exposure).
#  This VAR provides a complementary view to the SOFR-EFFR specification.
#
#  Ordering: d_TGA, d_Reserves, TED_spread, d_Aaa, d_Baa
###############################################################################

cat("================================================================\n")
cat("  VAR D: TED Spread (2003-2022)\n")
cat("================================================================\n\n")

df_ted <- df %>%
  filter(Date >= as.Date("2003-05-01"), Date <= as.Date("2022-01-21"))

var_cols_d <- c("d_TGA", "d_Reserves", "TED_spread", "d_Aaa", "d_Baa")
var_data_d <- df_ted %>%
  dplyr::select(all_of(var_cols_d)) %>%
  drop_na()

cat(sprintf("  VAR D observations: %d\n", nrow(var_data_d)))

# Lag selection
lag_sel_d <- VARselect(var_data_d, lag.max = max_lag, type = "const")
cat("  Lag selection:\n")
print(lag_sel_d$selection)
opt_lag_d <- lag_sel_d$selection["AIC(n)"]
cat(sprintf("  Selected lag: %d (AIC)\n\n", opt_lag_d))

# Estimate
var_d <- VAR(var_data_d, p = opt_lag_d, type = "const")

# Stability
roots_d <- roots(var_d)
cat(sprintf("  Max eigenvalue modulus: %.4f (Stable: %s)\n\n",
            max(roots_d), ifelse(max(roots_d) < 1, "YES", "NO")))

# R-squared
for (eq_name in names(var_d$varresult)) {
  r2 <- summary(var_d$varresult[[eq_name]])$adj.r.squared
  cat(sprintf("  R-squared (adj) for %s: %.4f\n", eq_name, r2))
}
cat("\n")

# Residual diagnostics
pt_test_d <- serial.test(var_d, lags.pt = 20, type = "PT.asymptotic")
cat(sprintf("  Portmanteau test (20 lags): Chi-sq = %.2f, p = %.4f\n\n",
            pt_test_d$serial$statistic, pt_test_d$serial$p.value))

# Granger causality within VAR
cat("  Granger causality in VAR framework:\n")
for (cause_var in c("d_TGA", "d_Reserves", "TED_spread")) {
  gc <- causality(var_d, cause = cause_var)
  cat(sprintf("    %s → others: F = %.3f, p = %.4f %s\n",
              cause_var,
              gc$Granger$statistic,
              gc$Granger$p.value,
              ifelse(gc$Granger$p.value < 0.05, "**", "")))
}
cat("\n")

# IRFs
irf_d <- irf(var_d, impulse = NULL, response = NULL,
              n.ahead = n_ahead, boot = TRUE, runs = n_boot,
              ci = 0.95, ortho = TRUE)

# Panel: TED spread IRFs
panel_d_df <- rbind(
  make_irf_panel_df(irf_d, "TED_spread", "d_Aaa",  "TED_spread -> d_Aaa"),
  make_irf_panel_df(irf_d, "TED_spread", "d_Baa",  "TED_spread -> d_Baa"),
  make_irf_panel_df(irf_d, "d_Reserves", "d_Aaa",  "d_Reserves -> d_Aaa"),
  make_irf_panel_df(irf_d, "d_Reserves", "d_Baa",  "d_Reserves -> d_Baa"),
  make_irf_panel_df(irf_d, "d_TGA",      "d_Aaa",  "d_TGA -> d_Aaa"),
  make_irf_panel_df(irf_d, "d_TGA",      "d_Baa",  "d_TGA -> d_Baa")
)

panel_d_df$Panel <- factor(panel_d_df$Panel,
  levels = c("TED_spread -> d_Aaa", "TED_spread -> d_Baa",
             "d_Reserves -> d_Aaa", "d_Reserves -> d_Baa",
             "d_TGA -> d_Aaa",      "d_TGA -> d_Baa"))

p_panel_d <- ggplot(panel_d_df, aes(x = Horizon)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "forestgreen", alpha = 0.15) +
  geom_line(aes(y = IRF), color = "forestgreen", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  facet_wrap(~ Panel, scales = "free_y", ncol = 2) +
  labs(title = "Impulse Responses — TED Spread VAR (2003-2022)",
       x = "Horizon (business days)", y = "Response") +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIG_DIR, "fig_irf_panel_ted.pdf"), p_panel_d, width = 8, height = 7, device = "pdf")
cat("  Saved: fig_irf_panel_ted.pdf\n\n")

# FEVD for TED VAR
fevd_d <- fevd(var_d, n.ahead = 30)
for (resp in c("d_Aaa", "d_Baa")) {
  cat(sprintf("  FEVD for %s (TED VAR):\n", resp))
  fevd_mat <- fevd_d[[resp]]
  for (h in c(1, 5, 10, 20, 30)) {
    if (h <= nrow(fevd_mat)) {
      cat(sprintf("    h=%2d: d_TGA=%.1f%% d_Res=%.1f%% TED=%.1f%% d_Aaa=%.1f%% d_Baa=%.1f%%\n",
                  h,
                  100*fevd_mat[h, "d_TGA"],
                  100*fevd_mat[h, "d_Reserves"],
                  100*fevd_mat[h, "TED_spread"],
                  100*fevd_mat[h, "d_Aaa"],
                  100*fevd_mat[h, "d_Baa"]))
    }
  }
  cat("\n")
}

# FEVD plot for TED VAR
fevd_d_plot <- data.frame()
for (resp in c("d_Aaa", "d_Baa")) {
  fevd_mat <- fevd_d[[resp]]
  for (h in 1:nrow(fevd_mat)) {
    for (src in var_cols_d) {
      fevd_d_plot <- rbind(fevd_d_plot, data.frame(
        Response = resp, Horizon = h, Source = src,
        Share = 100 * fevd_mat[h, src], stringsAsFactors = FALSE
      ))
    }
  }
}
fevd_d_plot$Source <- factor(fevd_d_plot$Source,
  levels = var_cols_d,
  labels = c("TGA", "Reserves", "TED spread", "Aaa", "Baa"))
fevd_d_plot$Response <- factor(fevd_d_plot$Response,
  levels = c("d_Aaa", "d_Baa"),
  labels = c("Variance of d_Aaa", "Variance of d_Baa"))

p_fevd_d <- ggplot(fevd_d_plot, aes(x = Horizon, y = Share, fill = Source)) +
  geom_area(alpha = 0.8) +
  facet_wrap(~ Response, ncol = 2) +
  scale_fill_brewer(palette = "Set2", name = "Shock source") +
  labs(x = "Horizon (business days)", y = "Share of forecast error variance (%)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "fig_fevd_ted.pdf"), p_fevd_d, width = 8, height = 4, device = "pdf")
cat("  Saved: fig_fevd_ted.pdf\n\n")


###############################################################################
#  SAVE RESULTS
###############################################################################

save(var_a, var_b, var_c, var_d, var_alt,
     irf_a, irf_b, irf_c, irf_d, irf_alt,
     fevd_a, fevd_b, fevd_c, fevd_d,
     lag_sel_a, lag_sel_b, lag_sel_c, lag_sel_d,
     file = file.path(OUT_DIR, "var_results.RData"))
cat("  Saved: var_results.RData\n\n")


###############################################################################
#  INTERPRETATION SUMMARY
###############################################################################

cat("================================================================\n")
cat("  INTERPRETATION SUMMARY\n")
cat("================================================================\n\n")

cat("  1. LAG ORDER:\n")
cat(sprintf("     VAR A (post-2018): AIC selects %d lags\n", opt_lag_a))
cat(sprintf("     VAR C (extended):  AIC selects %d lags\n", opt_lag_c))
cat(sprintf("     VAR D (TED 2003-2022): AIC selects %d lags\n\n", opt_lag_d))

cat("  2. TED SPREAD VAR (Key comparison):\n")
cat(sprintf("     TED_spread R2 = %.4f (vs SOFR_EFFR post-2018 R2 = %.4f)\n",
            summary(var_d$varresult[["TED_spread"]])$adj.r.squared,
            summary(var_a$varresult[["SOFR_EFFR"]])$adj.r.squared))
cat(sprintf("     d_Baa R2 in TED VAR = %.4f (vs %.4f in SOFR VAR)\n\n",
            summary(var_d$varresult[["d_Baa"]])$adj.r.squared,
            summary(var_a$varresult[["d_Baa"]])$adj.r.squared))

cat("  3. VARIANCE DECOMPOSITION COMPARISON:\n")
cat(sprintf("     TED → d_Baa at h=10: %.1f%% (vs SOFR → d_Baa: %.1f%%)\n",
            100*fevd_d[["d_Baa"]][10, "TED_spread"],
            100*fevd_a[["d_Baa"]][10, "SOFR_EFFR"]))
cat(sprintf("     TED → d_Aaa at h=10: %.1f%% (vs SOFR → d_Aaa: %.1f%%)\n\n",
            100*fevd_d[["d_Aaa"]][10, "TED_spread"],
            100*fevd_a[["d_Aaa"]][10, "SOFR_EFFR"]))

cat("================================================================\n")
cat("  STEP 3 COMPLETE\n")
cat("================================================================\n")
