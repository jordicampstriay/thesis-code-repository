###############################################################################
#  STEP 1 — DATA PREPARATION (Sections 5.1 & 5.2)
#  Thesis: "The Plumbing of Liquidity"
#  Author: Jordi Camps Triay
#  Last updated: 2026-05-11
#
#  This script is fully self-contained and reproducible.
#  Input : raw_data_clean.csv  (exported from "RAW DATA.xlsx", first sheet)
#  Output: data_master.csv          — cleaned panel, ready for analysis
#          table_data_sources.tex    — variable/source catalogue
#          table_summary_full.tex    — summary stats, full sample
#          table_summary_post2018.tex— summary stats, post-2018 sample
#          table_data_availability.tex — data coverage per variable
#          table_adf_results.tex     — ADF unit root tests
#          table_correlation_fd.tex  — correlation matrix in first differences
#          fig_ts_*.pdf              — time-series diagnostic plots
#          fig_corr_*.pdf            — correlation heat maps
###############################################################################

# ── 0.  PACKAGES ─────────────────────────────────────────────────────────────
required_pkgs <- c("readr", "dplyr", "tidyr", "zoo", "tseries", "ggplot2",
                   "scales", "corrplot", "xtable", "lubridate", "gridExtra",
                   "grDevices")

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

setwd("/Users/jordi/Downloads/University/TFG/Data/R")

cat("================================================================\n")
cat("  STEP 1 — DATA PREPARATION\n")
cat("================================================================\n\n")

# ── 1.  LOAD RAW DATA ────────────────────────────────────────────────────────
raw <- read_csv("raw_data_clean.csv",
                na = c("", "NA", "-"),
                show_col_types = FALSE)

# Standardise column names
colnames(raw) <- c("Date", "SOFR", "EFFR", "OBFR", "TED_spread", "DXY",
                    "HYG", "LQD", "TGCR", "Liq_Swaps", "M2",
                    "BGCR", "ON_RRP", "TGA", "Fed_Assets",
                    "Reserves", "Baa_spread", "Aaa_spread",
                    "TLT", "IEF", "MOVE", "SOFR_EFFR",
                    "SOFR_OBFR", "Aaa_Baa", "BGCR_TGCR")

raw$Date <- as.Date(raw$Date)

# Force all non-Date columns to numeric
for (col in colnames(raw)[-1]) {
  raw[[col]] <- suppressWarnings(as.numeric(raw[[col]]))
}

# Sort chronologically (raw file is reverse-chronological)
raw <- raw %>% arrange(Date)

cat(sprintf("Raw data loaded: %d rows, %d columns\n", nrow(raw), ncol(raw)))
cat(sprintf("Date range: %s to %s\n\n",
            min(raw$Date, na.rm = TRUE), max(raw$Date, na.rm = TRUE)))


# ── 2.  DATA AVAILABILITY AUDIT ─────────────────────────────────────────────
# For each variable, find the first and last non-NA date and count valid obs.
# This is essential for defining nested sample periods.

cat("================================================================\n")
cat("  DATA AVAILABILITY AUDIT\n")
cat("================================================================\n\n")

vars_of_interest <- c("SOFR_EFFR", "TED_spread", "Baa_spread", "Aaa_spread",
                       "Aaa_Baa", "DXY", "HYG", "LQD", "TLT", "IEF",
                       "Reserves", "TGA", "Fed_Assets", "ON_RRP",
                       "Liq_Swaps", "MOVE", "SOFR", "EFFR", "OBFR",
                       "TGCR", "BGCR", "M2", "SOFR_OBFR", "BGCR_TGCR")

availability <- data.frame(
  Variable   = character(),
  First_date = as.Date(character()),
  Last_date  = as.Date(character()),
  N_valid    = integer(),
  N_total    = integer(),
  Pct_valid  = numeric(),
  stringsAsFactors = FALSE
)

for (v in vars_of_interest) {
  if (v %in% colnames(raw)) {
    valid_idx <- which(!is.na(raw[[v]]))
    if (length(valid_idx) > 0) {
      availability <- rbind(availability, data.frame(
        Variable   = v,
        First_date = raw$Date[min(valid_idx)],
        Last_date  = raw$Date[max(valid_idx)],
        N_valid    = length(valid_idx),
        N_total    = nrow(raw),
        Pct_valid  = round(length(valid_idx) / nrow(raw) * 100, 1)
      ))
    }
  }
}

print(availability, row.names = FALSE)
cat("\n")

# ── 3.  DATA CLEANING ────────────────────────────────────────────────────────

cat("================================================================\n")
cat("  DATA CLEANING\n")
cat("================================================================\n\n")

df <- raw  # work on a copy

# 3a. Holiday zeros — days when both Aaa and Baa spreads are exactly zero.
#     These are bond-market holidays where FRED carries forward a zero.
#     Setting to NA prevents spurious zero-returns.
holiday_mask <- !is.na(df$Baa_spread) & df$Baa_spread == 0 &
                !is.na(df$Aaa_spread) & df$Aaa_spread == 0
n_holidays <- sum(holiday_mask)

market_cols <- c("Baa_spread", "Aaa_spread", "Aaa_Baa",
                 "DXY", "HYG", "LQD", "TLT", "IEF",
                 "TED_spread", "MOVE")
for (col in market_cols) {
  df[[col]][holiday_mask] <- NA
}
cat(sprintf("  Holiday zeros set to NA: %d rows\n", n_holidays))

# 3b. DXY = 0 on non-trading days (separate from bond holidays)
dxy_zero <- !is.na(df$DXY) & df$DXY == 0
n_dxy_zero <- sum(dxy_zero)
df$DXY[dxy_zero] <- NA
cat(sprintf("  Additional DXY zeros set to NA: %d rows\n", n_dxy_zero))

# 3c. Drop SOFR_OBFR — correlation 0.997 with SOFR_EFFR (near-perfect
#     multicollinearity; keeping both inflates SE and makes individual
#     coefficients uninterpretable).
df$SOFR_OBFR <- NULL
cat("  Dropped: SOFR_OBFR (correlation 0.997 with SOFR_EFFR)\n")

# 3d. Drop BGCR_TGCR — 97% of observations are zero.
#     No variation to exploit; contributes only noise.
df$BGCR_TGCR <- NULL
cat("  Dropped: BGCR_TGCR (97% zeros, no exploitable variation)\n")

# 3e. Drop M2 — monthly frequency, cannot be meaningfully interpolated
#     to daily without introducing measurement error.
df$M2 <- NULL
cat("  Dropped: M2 (monthly frequency, unusable at daily)\n")

# 3f. Drop individual rate levels (SOFR, EFFR, OBFR, TGCR, BGCR) —
#     we use the constructed spreads, not the levels.
df$SOFR <- NULL
df$EFFR <- NULL
df$OBFR <- NULL
df$TGCR <- NULL
df$BGCR <- NULL
cat("  Dropped: individual rate levels (SOFR, EFFR, OBFR, TGCR, BGCR)\n")
cat("           — the analysis uses spread variables, not raw rates\n")

cat(sprintf("\n  Remaining columns: %d\n", ncol(df)))
cat(sprintf("  Columns: %s\n\n", paste(colnames(df), collapse = ", ")))

# 3g. Forward-fill weekly series.
#     Reserves, TGA, Fed_Assets, Liq_Swaps are released on Wednesdays.
#     Each release is carried through the following business days.
#     Precedent: Adrian, Boyarchenko and Shachar (2017).
weekly_cols <- c("Reserves", "TGA", "Fed_Assets", "Liq_Swaps")
for (col in weekly_cols) {
  n_before <- sum(is.na(df[[col]]))
  df[[col]] <- na.locf(df[[col]], na.rm = FALSE)
  n_after <- sum(is.na(df[[col]]))
  cat(sprintf("  Forward-filled %s: %d NA -> %d NA\n",
              col, n_before, n_after))
}
cat("\n")


# ── 4.  CONSTRUCT DERIVED VARIABLES ──────────────────────────────────────────

cat("================================================================\n")
cat("  VARIABLE CONSTRUCTION\n")
cat("================================================================\n\n")

# 4a. Baa-Aaa quality differential
#     NOTE: the raw data has Aaa_Baa = Aaa - Baa (negative when Baa > Aaa).
#     For interpretability, we define Baa_Aaa = Baa_spread - Aaa_spread > 0.
df$Baa_Aaa <- df$Baa_spread - df$Aaa_spread
cat("  Created: Baa_Aaa = Baa_spread - Aaa_spread (quality premium, positive)\n")

# Drop the original Aaa_Baa column to avoid confusion
df$Aaa_Baa <- NULL

# 4b. First differences for non-stationary quantity variables
df <- df %>%
  mutate(
    d_Reserves   = c(NA, diff(Reserves)),
    d_TGA        = c(NA, diff(TGA)),
    d_ON_RRP     = c(NA, diff(ON_RRP)),
    d_Fed_Assets = c(NA, diff(Fed_Assets))
  )
cat("  Created: d_Reserves, d_TGA, d_ON_RRP, d_Fed_Assets (first differences)\n")

# 4c. First differences for credit spreads
df <- df %>%
  mutate(
    d_Baa     = c(NA, diff(Baa_spread)),
    d_Aaa     = c(NA, diff(Aaa_spread)),
    d_Baa_Aaa = c(NA, diff(Baa_Aaa))
  )
cat("  Created: d_Baa, d_Aaa, d_Baa_Aaa (first differences of credit spreads)\n")

# 4d. Log-returns for price variables
#     Log-returns are approximately normal, additive over time, and more
#     suitable for regression analysis than simple price differences.
df <- df %>%
  mutate(
    lr_HYG = c(NA, diff(log(HYG))) * 100,
    lr_LQD = c(NA, diff(log(LQD))) * 100,
    lr_DXY = c(NA, diff(log(DXY))) * 100,
    lr_TLT = c(NA, diff(log(TLT))) * 100,
    lr_IEF = c(NA, diff(log(IEF))) * 100
  )
cat("  Created: lr_HYG, lr_LQD, lr_DXY, lr_TLT, lr_IEF (log-returns, x100)\n")

# 4e. HYG/LQD risk appetite ratio (log)
df <- df %>%
  mutate(
    log_HYG_LQD   = log(HYG / LQD),
    d_log_HYG_LQD = c(NA, diff(log_HYG_LQD))
  )
cat("  Created: log_HYG_LQD, d_log_HYG_LQD (risk appetite proxy)\n\n")


# ── 5.  DEFINE SAMPLE PERIODS ────────────────────────────────────────────────

cat("================================================================\n")
cat("  SAMPLE DEFINITIONS\n")
cat("================================================================\n\n")

# Full sample: starts when TED spread and credit spreads are jointly available
df_full <- df %>% filter(Date >= as.Date("2003-05-01"))

# Post-2018 sample: starts when SOFR publication began (April 3, 2018)
df_post <- df %>% filter(Date >= as.Date("2018-04-03"))

# Count business days with valid SOFR_EFFR in post-2018 sample
n_sofr_valid <- sum(!is.na(df_post$SOFR_EFFR))

cat(sprintf("  Full sample:      %s to %s  (%d rows)\n",
            min(df_full$Date), max(df_full$Date), nrow(df_full)))
cat(sprintf("  Post-2018 sample: %s to %s  (%d rows)\n",
            min(df_post$Date), max(df_post$Date), nrow(df_post)))
cat(sprintf("  Valid SOFR-EFFR obs in post-2018: %d\n\n", n_sofr_valid))


# ── 6.  SAVE MASTER DATA ─────────────────────────────────────────────────────

write.csv(df, "data_master.csv", row.names = FALSE)
cat("  Saved: data_master.csv\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#                         DIAGNOSTIC TABLES & FIGURES
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════


# ── 7.  DATA SOURCE TABLE (for Section 5.1) ──────────────────────────────────

cat("================================================================\n")
cat("  TABLE: DATA SOURCES\n")
cat("================================================================\n\n")

source_table <- data.frame(
  Variable = c(
    "SOFR--EFFR spread",
    "TED spread",
    "Moody's Baa corporate spread",
    "Moody's Aaa corporate spread",
    "Baa--Aaa quality differential",
    "Reserve balances at Fed",
    "TGA balance",
    "Total assets (Fed)",
    "ON Reverse Repo operations",
    "Fed dollar liquidity swaps",
    "DXY index",
    "HYG (iShares High Yield ETF)",
    "LQD (iShares IG Corporate ETF)",
    "MOVE index"
  ),
  Source = c(
    "FRBNY; FRED (SOFR, EFFR)",
    "FRED (TEDRATE)",
    "FRED (BAA10Y)",
    "FRED (AAA10Y)",
    "Constructed (Baa -- Aaa)",
    "FRED (WRESBAL)",
    "FRED (WTREGEN)",
    "FRED (WALCL)",
    "FRED (RRPONTSYD)",
    "FRED (SWPT)",
    "FRED (DTWEXBGS)",
    "Bloomberg (HYG US Equity)",
    "Bloomberg (LQD US Equity)",
    "Bloomberg (MOVE Index)"
  ),
  Frequency = c(
    "Daily", "Daily", "Daily", "Daily", "Daily",
    "Weekly (Wed)", "Weekly (Wed)", "Weekly (Wed)",
    "Daily", "Weekly (Wed)",
    "Daily", "Daily", "Daily", "Daily"
  ),
  Available = c(
    "Apr 2018--Apr 2026",
    "May 2003--Dec 2024",
    "May 2003--Apr 2026",
    "May 2003--Apr 2026",
    "May 2003--Apr 2026",
    "May 2003--Apr 2026",
    "May 2003--Apr 2026",
    "May 2003--Apr 2026",
    "May 2003--Apr 2026",
    "May 2003--Apr 2026",
    "Nov 2005--Apr 2026",
    "Apr 2007--Apr 2026",
    "Jul 2005--Apr 2026",
    "Oct 2020--Apr 2026"
  ),
  Role = c(
    "Funding input (price)",
    "Funding input (extended)",
    "Credit output",
    "Credit output",
    "Credit output (collateral)",
    "Funding input (quantity)",
    "Funding input (exogenous)",
    "Funding input (quantity)",
    "Funding input (quantity)",
    "Descriptive only",
    "Risk appetite extension",
    "Risk appetite extension",
    "Risk appetite extension",
    "Validation"
  ),
  stringsAsFactors = FALSE
)

# Print to console
print(source_table, row.names = FALSE, right = FALSE)

# Export LaTeX
source_xt <- xtable(
  source_table,
  caption = "Data Sources and Variable Descriptions",
  label   = "tab:data_sources",
  align   = c("l", "l", "l", "c", "c", "l")
)
print(source_xt,
      file = "table_data_sources.tex",
      include.rownames = FALSE,
      booktabs = TRUE,
      caption.placement = "top",
      sanitize.text.function = function(x) gsub("--", "\\\\textendash{}", x),
      tabular.environment = "tabularx",
      width = "\\textwidth",
      floating = TRUE,
      table.placement = "htbp",
      scalebox = 0.78)
cat("\n  Saved: table_data_sources.tex\n\n")


# ── 8.  DATA AVAILABILITY TABLE ──────────────────────────────────────────────

cat("================================================================\n")
cat("  TABLE: DATA AVAILABILITY\n")
cat("================================================================\n\n")

# Re-compute on cleaned data (before dropping columns)
avail_vars <- c("SOFR_EFFR", "TED_spread", "Baa_spread", "Aaa_spread",
                "Baa_Aaa", "DXY", "HYG", "LQD",
                "Reserves", "TGA", "Fed_Assets", "ON_RRP",
                "Liq_Swaps", "MOVE")

avail_df <- data.frame(
  Variable  = character(),
  First     = character(),
  Last      = character(),
  N_obs     = integer(),
  stringsAsFactors = FALSE
)

for (v in avail_vars) {
  if (v %in% colnames(df)) {
    valid <- which(!is.na(df[[v]]))
    if (length(valid) > 0) {
      avail_df <- rbind(avail_df, data.frame(
        Variable = v,
        First    = as.character(df$Date[min(valid)]),
        Last     = as.character(df$Date[max(valid)]),
        N_obs    = length(valid)
      ))
    }
  }
}

print(avail_df, row.names = FALSE)

avail_xt <- xtable(
  avail_df,
  caption = "Data Availability After Cleaning (Business Days)",
  label   = "tab:data_availability",
  digits  = c(0, 0, 0, 0, 0)
)
print(avail_xt,
      file = "table_data_availability.tex",
      include.rownames = FALSE,
      booktabs = TRUE,
      caption.placement = "top",
      sanitize.text.function = identity,
      table.placement = "htbp")
cat("\n  Saved: table_data_availability.tex\n\n")


# ── 9.  SUMMARY STATISTICS ───────────────────────────────────────────────────

cat("================================================================\n")
cat("  TABLE: SUMMARY STATISTICS\n")
cat("================================================================\n\n")

compute_summary <- function(data, vars, label) {
  out <- data.frame(
    Variable = character(), N = integer(),
    Mean = numeric(), SD = numeric(),
    Min = numeric(), Q25 = numeric(),
    Median = numeric(), Q75 = numeric(),
    Max = numeric(),
    stringsAsFactors = FALSE
  )
  for (v in vars) {
    if (v %in% colnames(data)) {
      x <- data[[v]]
      x <- x[!is.na(x)]
      if (length(x) > 0) {
        q <- quantile(x, probs = c(0.25, 0.50, 0.75))
        out <- rbind(out, data.frame(
          Variable = v, N = length(x),
          Mean = round(mean(x), 4), SD = round(sd(x), 4),
          Min = round(min(x), 4), Q25 = round(q[1], 4),
          Median = round(q[2], 4), Q75 = round(q[3], 4),
          Max = round(max(x), 4)
        ))
      }
    }
  }
  out
}

# Variables for the summary: levels of key series
summary_vars_levels <- c("SOFR_EFFR", "TED_spread",
                          "Baa_spread", "Aaa_spread", "Baa_Aaa",
                          "DXY", "HYG", "LQD", "MOVE")

summary_vars_quantities <- c("Reserves", "TGA", "Fed_Assets", "ON_RRP")

# Post-2018 summary
summ_post_levels <- compute_summary(df_post, summary_vars_levels, "Post-2018")
summ_post_quant  <- compute_summary(df_post, summary_vars_quantities, "Post-2018")

cat("--- Summary Statistics: Post-2018 Sample (Levels) ---\n")
print(summ_post_levels, row.names = FALSE)
cat("\n--- Summary Statistics: Post-2018 Sample (Quantities, millions $) ---\n")
print(summ_post_quant, row.names = FALSE)

# Full-sample summary
summ_full_levels <- compute_summary(df_full, summary_vars_levels, "Full")
summ_full_quant  <- compute_summary(df_full, summary_vars_quantities, "Full")

cat("\n--- Summary Statistics: Full Sample (Levels) ---\n")
print(summ_full_levels, row.names = FALSE)

# Combine for LaTeX export
summ_post_all <- rbind(summ_post_levels, summ_post_quant)

summ_xt <- xtable(
  summ_post_all,
  caption = "Summary Statistics --- Post-2018 Sample (April 2018 -- April 2026)",
  label   = "tab:summary_post2018",
  digits  = c(0, 0, 0, 4, 4, 4, 4, 4, 4, 4)
)
print(summ_xt,
      file = "table_summary_post2018.tex",
      include.rownames = FALSE,
      booktabs = TRUE,
      caption.placement = "top",
      sanitize.text.function = function(x) gsub("_", "\\\\_", x),
      table.placement = "htbp",
      scalebox = 0.82)

summ_full_all <- rbind(summ_full_levels, summ_full_quant)
summ_full_xt <- xtable(
  summ_full_all,
  caption = "Summary Statistics --- Full Sample (May 2003 -- April 2026)",
  label   = "tab:summary_full",
  digits  = c(0, 0, 0, 4, 4, 4, 4, 4, 4, 4)
)
print(summ_full_xt,
      file = "table_summary_full.tex",
      include.rownames = FALSE,
      booktabs = TRUE,
      caption.placement = "top",
      sanitize.text.function = function(x) gsub("_", "\\\\_", x),
      table.placement = "htbp",
      scalebox = 0.82)

cat("\n  Saved: table_summary_post2018.tex\n")
cat("  Saved: table_summary_full.tex\n\n")


# ── 10. AUGMENTED DICKEY-FULLER TESTS ────────────────────────────────────────

cat("================================================================\n")
cat("  TABLE: ADF UNIT ROOT TESTS\n")
cat("  H0: series has a unit root (non-stationary)\n")
cat("  Reject at p < 0.05 => stationary\n")
cat("================================================================\n\n")

adf_vars <- c("SOFR_EFFR", "TED_spread",
               "Baa_spread", "Aaa_spread", "Baa_Aaa",
               "DXY", "HYG", "LQD",
               "Reserves", "TGA", "Fed_Assets", "ON_RRP",
               "MOVE")

adf_results <- data.frame(
  Variable        = character(),
  N               = integer(),
  ADF_stat_level  = numeric(),
  p_value_level   = numeric(),
  Decision_level  = character(),
  ADF_stat_diff   = numeric(),
  p_value_diff    = numeric(),
  Decision_diff   = character(),
  stringsAsFactors = FALSE
)

for (v in adf_vars) {
  # Use post-2018 sample for SOFR_EFFR, MOVE; full sample for the rest
  if (v %in% c("SOFR_EFFR", "MOVE")) {
    x <- df_post[[v]]
  } else {
    x <- df_full[[v]]
  }
  x <- x[!is.na(x)]
  if (length(x) < 50) next

  tryCatch({
    # Levels
    adf_lev <- adf.test(x, alternative = "stationary")
    # First differences
    dx <- diff(x)
    dx <- dx[!is.na(dx)]
    adf_dif <- adf.test(dx, alternative = "stationary")

    adf_results <- rbind(adf_results, data.frame(
      Variable        = v,
      N               = length(x),
      ADF_stat_level  = round(adf_lev$statistic, 3),
      p_value_level   = round(adf_lev$p.value, 4),
      Decision_level  = ifelse(adf_lev$p.value < 0.05, "Stationary", "Unit root"),
      ADF_stat_diff   = round(adf_dif$statistic, 3),
      p_value_diff    = round(adf_dif$p.value, 4),
      Decision_diff   = ifelse(adf_dif$p.value < 0.05, "Stationary", "Unit root")
    ))
  }, error = function(e) {
    cat(sprintf("  ADF error for %s: %s\n", v, e$message))
  })
}

print(adf_results, row.names = FALSE)

# Interpretation
cat("\n--- Stationary in LEVELS (enter VAR as-is): ---\n")
stat_lev <- adf_results$Variable[adf_results$Decision_level == "Stationary"]
cat(sprintf("  %s\n", paste(stat_lev, collapse = ", ")))

cat("\n--- NON-STATIONARY in levels (use first differences): ---\n")
nonstat <- adf_results$Variable[adf_results$Decision_level == "Unit root"]
cat(sprintf("  %s\n", paste(nonstat, collapse = ", ")))

all_diff_stat <- all(adf_results$Decision_diff[adf_results$Decision_level == "Unit root"] == "Stationary")
cat(sprintf("  All become stationary after differencing: %s\n\n", all_diff_stat))

# LaTeX export
adf_latex <- adf_results
colnames(adf_latex) <- c("Variable", "$N$", "ADF (Level)", "$p$-value",
                          "Decision", "ADF ($\\Delta$)", "$p$-value ",
                          "Decision ")
adf_xt <- xtable(
  adf_latex,
  caption = "Augmented Dickey--Fuller Unit Root Tests",
  label   = "tab:adf_tests",
  digits  = c(0, 0, 0, 3, 4, 0, 3, 4, 0)
)
print(adf_xt,
      file = "table_adf_results.tex",
      include.rownames = FALSE,
      booktabs = TRUE,
      caption.placement = "top",
      sanitize.colnames.function = identity,
      sanitize.text.function = function(x) gsub("_", "\\\\_", x),
      table.placement = "htbp",
      scalebox = 0.82)
cat("  Saved: table_adf_results.tex\n\n")


# ── 11. CORRELATION DIAGNOSTICS ──────────────────────────────────────────────

cat("================================================================\n")
cat("  CORRELATION ANALYSIS\n")
cat("================================================================\n\n")

# 11a. SOFR-EFFR vs SOFR-OBFR redundancy check (on raw data, before dropping)
sofr_effr_raw <- raw$SOFR_EFFR[!is.na(raw$SOFR_EFFR) & !is.na(raw$SOFR_OBFR)]
sofr_obfr_raw <- raw$SOFR_OBFR[!is.na(raw$SOFR_EFFR) & !is.na(raw$SOFR_OBFR)]
rho_redundancy <- cor(sofr_effr_raw, sofr_obfr_raw)
cat(sprintf("  SOFR-EFFR vs SOFR-OBFR correlation: %.4f\n", rho_redundancy))
cat("  >> Confirms near-perfect multicollinearity. SOFR-OBFR correctly excluded.\n\n")

# 11b. Correlation matrix of LEVELS (funding + credit spreads), post-2018
corr_vars_levels <- c("SOFR_EFFR", "TED_spread",
                       "Baa_spread", "Aaa_spread", "Baa_Aaa")

corr_level_data <- df_post %>%
  select(all_of(corr_vars_levels)) %>%
  drop_na()

corr_level_mat <- cor(corr_level_data, use = "pairwise.complete.obs")
cat("--- Correlation: Spread Variables in Levels (Post-2018) ---\n")
print(round(corr_level_mat, 3))
cat("\n")

# 11c. Correlation matrix of FIRST DIFFERENCES — the key diagnostic.
#      This is what matters for regression / VAR analysis.
corr_vars_fd <- c("SOFR_EFFR", "d_Reserves", "d_TGA", "d_ON_RRP",
                   "d_Baa", "d_Aaa", "d_Baa_Aaa")

corr_fd_data <- df_post %>%
  select(all_of(corr_vars_fd)) %>%
  drop_na()

corr_fd_mat <- cor(corr_fd_data, use = "pairwise.complete.obs")
cat("--- Correlation: Key Variables in Stationary Transformations (Post-2018) ---\n")
print(round(corr_fd_mat, 3))
cat("\n")

# Key observation: near-zero contemporaneous correlation between funding
# and credit in first differences => transmission operates with a LAG
cat("  Key finding: SOFR_EFFR vs d_Baa correlation =",
    round(corr_fd_mat["SOFR_EFFR", "d_Baa"], 4), "\n")
cat("  Key finding: SOFR_EFFR vs d_Aaa correlation =",
    round(corr_fd_mat["SOFR_EFFR", "d_Aaa"], 4), "\n")
cat("  Key finding: d_Reserves vs d_Aaa correlation =",
    round(corr_fd_mat["d_Reserves", "d_Aaa"], 4), "\n")
cat("  >> Near-zero contemporaneous correlations suggest the transmission\n")
cat("     operates with a lag => Granger causality and VAR are the correct tools.\n\n")

# LaTeX export — correlation in first differences
corr_fd_rounded <- round(corr_fd_mat, 3)
# Nice variable labels
rownames(corr_fd_rounded) <- c("SOFR-EFFR", "$\\Delta$Reserves", "$\\Delta$TGA",
                                 "$\\Delta$ON RRP",
                                 "$\\Delta$Baa", "$\\Delta$Aaa", "$\\Delta$(Baa-Aaa)")
colnames(corr_fd_rounded) <- rownames(corr_fd_rounded)

corr_xt <- xtable(
  corr_fd_rounded,
  caption = "Correlation Matrix --- Key Variables in Stationary Transformations (Post-2018 Sample)",
  label   = "tab:correlation_fd",
  digits  = 3
)
print(corr_xt,
      file = "table_correlation_fd.tex",
      booktabs = TRUE,
      caption.placement = "top",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      sanitize.rownames.function = identity,
      table.placement = "htbp",
      scalebox = 0.78)
cat("  Saved: table_correlation_fd.tex\n\n")

# 11d. Full-sample correlation (TED spread as funding proxy)
corr_vars_full <- c("TED_spread", "d_Baa", "d_Aaa", "d_Baa_Aaa",
                     "d_Reserves", "d_TGA")
corr_full_data <- df_full %>%
  select(all_of(corr_vars_full)) %>%
  drop_na()

corr_full_mat <- cor(corr_full_data, use = "pairwise.complete.obs")
cat("--- Correlation: Full Sample with TED Spread ---\n")
print(round(corr_full_mat, 3))
cat("\n")


# ── 12. CORRELATION HEAT MAP PLOTS ──────────────────────────────────────────

cat("================================================================\n")
cat("  CORRELATION HEAT MAP PLOTS\n")
cat("================================================================\n\n")

# Post-2018 first differences
pdf("fig_corr_fd_post2018.pdf", width = 8, height = 7)
par(mar = c(1, 1, 3, 1))
corrplot(corr_fd_mat, method = "color", type = "upper",
         addCoef.col = "black", number.cex = 0.75,
         tl.col = "black", tl.srt = 45, tl.cex = 0.85,
         title = "",
         col = colorRampPalette(c("#B71C1C", "#FFCDD2", "white",
                                   "#BBDEFB", "#0D47A1"))(200))
mtext("Correlation Matrix: Stationary Transformations (Post-2018)", side = 3,
      line = 1.5, cex = 1.1, font = 2)
dev.off()
cat("  Saved: fig_corr_fd_post2018.pdf\n")

# Full-sample with TED
pdf("fig_corr_fd_fullsample.pdf", width = 7, height = 6)
par(mar = c(1, 1, 3, 1))
corrplot(corr_full_mat, method = "color", type = "upper",
         addCoef.col = "black", number.cex = 0.8,
         tl.col = "black", tl.srt = 45, tl.cex = 0.9,
         title = "",
         col = colorRampPalette(c("#B71C1C", "#FFCDD2", "white",
                                   "#BBDEFB", "#0D47A1"))(200))
mtext("Correlation Matrix: Stationary Transformations (Full Sample)", side = 3,
      line = 1.5, cex = 1.1, font = 2)
dev.off()
cat("  Saved: fig_corr_fd_fullsample.pdf\n\n")


# ── 13. TIME SERIES DIAGNOSTIC PLOTS ────────────────────────────────────────

cat("================================================================\n")
cat("  TIME SERIES PLOTS\n")
cat("================================================================\n\n")

# Common thesis-quality theme
theme_thesis <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(face = "bold", size = 13,
                                       margin = margin(b = 4)),
      plot.subtitle    = element_text(color = "grey40", size = 9,
                                       margin = margin(b = 8)),
      plot.caption     = element_text(size = 7.5, color = "grey50",
                                       hjust = 0, margin = margin(t = 8)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      axis.title       = element_text(size = 10, color = "grey30"),
      axis.text        = element_text(size = 9),
      legend.position  = "bottom",
      legend.text      = element_text(size = 9),
      plot.margin      = margin(12, 14, 8, 10)
    )
}

# Colour palette
COL_NAVY  <- "#1B4F72"
COL_RED   <- "#922B21"
COL_GREEN <- "#1E8449"
COL_ORANGE <- "#D4AC0D"
COL_GREY   <- "grey60"

# NBER recession bands
rec_bands <- data.frame(
  start = as.Date(c("2007-12-01", "2020-02-01")),
  end   = as.Date(c("2009-06-01", "2020-04-01"))
)

add_recessions <- function(p, from_date = as.Date("2003-01-01")) {
  bands <- rec_bands[rec_bands$end >= from_date, ]
  if (nrow(bands) > 0) {
    p <- p + geom_rect(data = bands,
                        aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                        fill = "grey80", alpha = 0.4, inherit.aes = FALSE)
  }
  p
}

# Stress episode bands (for post-2018 plots)
stress_bands <- data.frame(
  start = as.Date(c("2018-12-01", "2019-09-16", "2020-03-09", "2022-09-01")),
  end   = as.Date(c("2018-12-31", "2019-10-15", "2020-04-15", "2022-12-31")),
  label = c("Q4 2018\nselloff", "Sept 2019\nrepo", "COVID-19\ndash for cash",
            "QT\nstress")
)

add_stress <- function(p) {
  p + geom_rect(data = stress_bands,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                fill = "#FFCCCC", alpha = 0.35, inherit.aes = FALSE)
}

# ─── PLOT 1: SOFR-EFFR spread ───
p1 <- ggplot(df_post %>% filter(!is.na(SOFR_EFFR)),
             aes(x = Date, y = SOFR_EFFR * 100)) +
  geom_hline(yintercept = 0, color = COL_GREY, linewidth = 0.4, linetype = "dashed")

p1 <- add_stress(p1)
p1 <- p1 +
  geom_line(color = COL_NAVY, linewidth = 0.35, alpha = 0.85) +
  labs(title = "SOFR-EFFR Spread",
       subtitle = "Daily, April 2018 -- April 2026. Shaded areas: stress episodes.",
       x = NULL, y = "Basis points",
       caption = "Source: Federal Reserve Bank of New York, FRED.") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",
               expand = c(0.01, 0)) +
  theme_thesis()

# ─── PLOT 2: Credit spreads (Baa and Aaa) ───
p2_data <- df_full %>%
  select(Date, Baa_spread, Aaa_spread) %>%
  pivot_longer(-Date, names_to = "Series", values_to = "Value") %>%
  filter(!is.na(Value))

p2 <- ggplot(p2_data, aes(x = Date, y = Value, color = Series))
p2 <- add_recessions(p2)
p2 <- p2 +
  geom_line(linewidth = 0.35, alpha = 0.85) +
  scale_color_manual(
    values = c("Aaa_spread" = COL_NAVY, "Baa_spread" = COL_RED),
    labels = c("Aaa_spread" = "Moody's Aaa", "Baa_spread" = "Moody's Baa")
  ) +
  labs(title = "Corporate Credit Spreads",
       subtitle = "Moody's Baa and Aaa spreads over 10-year Treasury. Shaded: NBER recessions.",
       x = NULL, y = "Spread (percentage points)", color = NULL,
       caption = "Source: FRED (BAA10Y, AAA10Y).") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               expand = c(0.01, 0)) +
  theme_thesis()

# ─── PLOT 3: Baa-Aaa quality differential ───
p3 <- ggplot(df_full %>% filter(!is.na(Baa_Aaa)),
             aes(x = Date, y = Baa_Aaa))
p3 <- add_recessions(p3)
p3 <- p3 +
  geom_line(color = COL_GREEN, linewidth = 0.35, alpha = 0.85) +
  labs(title = "Baa -- Aaa Quality Differential",
       subtitle = "Positive = Baa wider than Aaa (normal). Higher = greater credit differentiation.",
       x = NULL, y = "Spread differential (pp)",
       caption = "Source: Constructed from FRED data.") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               expand = c(0.01, 0)) +
  theme_thesis()

# ─── PLOT 4: Reserves and TGA ───
p4_data <- df_full %>%
  select(Date, Reserves, TGA) %>%
  filter(!is.na(Reserves)) %>%
  mutate(Reserves = Reserves / 1e6, TGA = TGA / 1e6) %>%
  pivot_longer(-Date, names_to = "Series", values_to = "Value") %>%
  filter(!is.na(Value))

p4 <- ggplot(p4_data, aes(x = Date, y = Value, color = Series))
p4 <- add_recessions(p4)
p4 <- p4 +
  geom_line(linewidth = 0.45, alpha = 0.85) +
  scale_color_manual(values = c("Reserves" = COL_NAVY, "TGA" = COL_ORANGE)) +
  labs(title = "Reserve Balances and Treasury General Account",
       subtitle = "Weekly (forward-filled to daily). In trillions of USD.",
       x = NULL, y = "USD trillions", color = NULL,
       caption = "Source: FRED (WRESBAL, WTREGEN).") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               expand = c(0.01, 0)) +
  theme_thesis()

# ─── PLOT 5: TED spread (full sample) ───
p5 <- ggplot(df_full %>% filter(!is.na(TED_spread)),
             aes(x = Date, y = TED_spread * 100))
p5 <- add_recessions(p5)
p5 <- p5 +
  geom_line(color = COL_RED, linewidth = 0.35, alpha = 0.85) +
  labs(title = "TED Spread",
       subtitle = "3-month LIBOR minus 3-month Treasury bill rate. Shaded: NBER recessions.",
       x = NULL, y = "Basis points",
       caption = "Source: FRED (TEDRATE).") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               expand = c(0.01, 0)) +
  theme_thesis()

# ─── PLOT 6: HYG and LQD ETF prices ───
p6_data <- df_full %>%
  select(Date, HYG, LQD) %>%
  pivot_longer(-Date, names_to = "Series", values_to = "Value") %>%
  filter(!is.na(Value))

p6 <- ggplot(p6_data, aes(x = Date, y = Value, color = Series))
p6 <- add_recessions(p6)
p6 <- p6 +
  geom_line(linewidth = 0.35, alpha = 0.85) +
  scale_color_manual(values = c("HYG" = COL_RED, "LQD" = COL_NAVY)) +
  labs(title = "Bond ETF Prices: HYG and LQD",
       subtitle = "Daily closing prices. HYG = high yield, LQD = investment grade.",
       x = NULL, y = "Price (USD)", color = NULL,
       caption = "Source: Bloomberg.") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               expand = c(0.01, 0)) +
  theme_thesis()

# ─── PLOT 7: ON Reverse Repo ───
p7 <- ggplot(df_full %>% filter(!is.na(ON_RRP), Date >= as.Date("2013-01-01")),
             aes(x = Date, y = ON_RRP))
p7 <- p7 +
  geom_line(color = "#7B1FA2", linewidth = 0.4, alpha = 0.85) +
  labs(title = "Overnight Reverse Repo Facility",
       subtitle = "Daily take-up, in billions of USD.",
       x = NULL, y = "USD billions",
       caption = "Source: FRED (RRPONTSYD).") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               expand = c(0.01, 0)) +
  theme_thesis()

# ─── PLOT 8: DXY index ───
p8 <- ggplot(df_full %>% filter(!is.na(DXY)),
             aes(x = Date, y = DXY))
p8 <- add_recessions(p8)
p8 <- p8 +
  geom_line(color = COL_GREEN, linewidth = 0.35, alpha = 0.85) +
  labs(title = "US Dollar Index (DXY)",
       subtitle = "Daily. Higher = stronger dollar = tighter global USD funding.",
       x = NULL, y = "Index level",
       caption = "Source: Bloomberg.") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               expand = c(0.01, 0)) +
  theme_thesis()

# ─── PLOT 9: Funding vs Credit overlay (motivational) ───
p9_data <- df_post %>%
  filter(!is.na(SOFR_EFFR) & !is.na(Baa_spread))

p9 <- ggplot(p9_data)
p9 <- add_stress(p9)
p9 <- p9 +
  geom_line(aes(x = Date, y = SOFR_EFFR * 100, color = "SOFR-EFFR (bps)"),
            linewidth = 0.4, alpha = 0.85) +
  geom_line(aes(x = Date,
                y = (Baa_spread - mean(Baa_spread, na.rm = TRUE)) * 30,
                color = "Baa spread (demeaned, rescaled)"),
            linewidth = 0.4, alpha = 0.7) +
  scale_color_manual(values = c("SOFR-EFFR (bps)" = COL_NAVY,
                                 "Baa spread (demeaned, rescaled)" = COL_RED)) +
  labs(title = "Funding Spread vs Credit Spread",
       subtitle = "Post-2018 sample. Visual inspection for regime-dependent co-movement.",
       x = NULL, y = "SOFR-EFFR (bps)", color = NULL,
       caption = "Note: Baa spread demeaned and rescaled for visual comparison only.") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",
               expand = c(0.01, 0)) +
  theme_thesis()

# ─── Save all plots ───
plots <- list(
  list(p = p1, name = "fig_ts_sofr_effr.pdf",     w = 10, h = 4.5),
  list(p = p2, name = "fig_ts_credit_spreads.pdf", w = 10, h = 4.5),
  list(p = p3, name = "fig_ts_baa_aaa_diff.pdf",   w = 10, h = 4.5),
  list(p = p4, name = "fig_ts_reserves_tga.pdf",    w = 10, h = 4.5),
  list(p = p5, name = "fig_ts_ted_spread.pdf",      w = 10, h = 4.5),
  list(p = p6, name = "fig_ts_etf_prices.pdf",      w = 10, h = 4.5),
  list(p = p7, name = "fig_ts_on_rrp.pdf",          w = 10, h = 4.5),
  list(p = p8, name = "fig_ts_dxy.pdf",              w = 10, h = 4.5),
  list(p = p9, name = "fig_ts_funding_vs_credit.pdf", w = 10, h = 4.5)
)

for (item in plots) {
  ggsave(item$name, item$p, width = item$w, height = item$h, dpi = 300)
  cat(sprintf("  Saved: %s\n", item$name))
}
cat("\n")


# ── 14. ROLLING VOLATILITY DIAGNOSTIC ───────────────────────────────────────

cat("================================================================\n")
cat("  ROLLING VOLATILITY: Aaa vs Baa\n")
cat("================================================================\n\n")

df_vol <- df_full %>%
  select(Date, d_Aaa, d_Baa) %>%
  filter(!is.na(d_Aaa) & !is.na(d_Baa))

window_size <- 60  # 60 business days ~ 3 months

df_vol <- df_vol %>%
  mutate(
    vol_Aaa   = zoo::rollapply(d_Aaa, width = window_size, FUN = sd,
                                fill = NA, align = "right"),
    vol_Baa   = zoo::rollapply(d_Baa, width = window_size, FUN = sd,
                                fill = NA, align = "right"),
    vol_ratio = vol_Aaa / vol_Baa
  )

# Plot: volatility ratio
p_vol <- ggplot(df_vol %>% filter(!is.na(vol_ratio)),
                aes(x = Date, y = vol_ratio))
p_vol <- add_recessions(p_vol)
p_vol <- p_vol +
  geom_hline(yintercept = 1, color = COL_GREY, linewidth = 0.4, linetype = "dashed") +
  geom_line(color = COL_GREEN, linewidth = 0.4, alpha = 0.85) +
  annotate("rect",
           xmin = as.Date("2019-09-15"), xmax = as.Date("2019-10-15"),
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#E57373") +
  annotate("rect",
           xmin = as.Date("2020-03-01"), xmax = as.Date("2020-04-15"),
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#E57373") +
  labs(title = "Volatility Ratio: Aaa / Baa (60-Day Rolling Window)",
       subtitle = "Above 1 = Aaa more volatile than Baa (unusual, signals collateral-channel stress).",
       x = NULL, y = expression(sigma[Aaa] / sigma[Baa]),
       caption = "Source: Own calculations from FRED data. Red bands: Sept 2019, March 2020.") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               expand = c(0.01, 0)) +
  coord_cartesian(ylim = c(0, max(df_vol$vol_ratio, na.rm = TRUE) * 1.05)) +
  theme_thesis()

ggsave("fig_ts_vol_ratio.pdf", p_vol, width = 10, height = 4.5, dpi = 300)
cat("  Saved: fig_ts_vol_ratio.pdf\n\n")


# ── 15. MULTICOLLINEARITY CHECK (VIF-like) ──────────────────────────────────

cat("================================================================\n")
cat("  MULTICOLLINEARITY DIAGNOSTIC\n")
cat("================================================================\n\n")

# Check pairwise correlations among funding inputs that will enter the VAR
var_candidates <- c("SOFR_EFFR", "d_Reserves", "d_TGA")
mc_data <- df_post %>%
  select(all_of(var_candidates)) %>%
  drop_na()

mc_cor <- cor(mc_data)
cat("Pairwise correlations among VAR funding inputs (post-2018):\n")
print(round(mc_cor, 4))
cat("\n")

# All pairwise |r| < 0.30 => no multicollinearity concern
max_offdiag <- max(abs(mc_cor[upper.tri(mc_cor)]))
cat(sprintf("  Maximum off-diagonal |r|: %.4f\n", max_offdiag))
if (max_offdiag < 0.50) {
  cat("  >> No multicollinearity concern among VAR funding inputs.\n\n")
} else {
  cat("  >> WARNING: elevated pairwise correlation. Investigate before VAR.\n\n")
}


# ── 16. ZERO-VARIATION CHECK ON SOFR-EFFR ───────────────────────────────────

cat("================================================================\n")
cat("  SOFR-EFFR DISTRIBUTION DIAGNOSTIC\n")
cat("================================================================\n\n")

sofr_clean <- df_post$SOFR_EFFR[!is.na(df_post$SOFR_EFFR)]
n_near_zero <- sum(abs(sofr_clean) < 0.005)
pct_near_zero <- round(n_near_zero / length(sofr_clean) * 100, 1)

cat(sprintf("  Total SOFR-EFFR observations: %d\n", length(sofr_clean)))
cat(sprintf("  Near-zero (|spread| < 0.5 bp): %d (%.1f%%)\n",
            n_near_zero, pct_near_zero))
cat(sprintf("  Mean: %.4f  |  SD: %.4f  |  Skewness metric: %.4f\n",
            mean(sofr_clean), sd(sofr_clean),
            mean((sofr_clean - mean(sofr_clean))^3) / sd(sofr_clean)^3))
cat("  >> Concentrated near zero with fat tails.\n")
cat("     This motivates the regime-dependent analysis in Section 5.5:\n")
cat("     a linear Granger test averages rare spikes with the ~85% of days\n")
cat("     when the spread is effectively zero.\n\n")


# ── 17. SAMPLE PERIOD DOCUMENTATION ─────────────────────────────────────────

cat("================================================================\n")
cat("  SAMPLE PERIOD SUMMARY FOR THESIS\n")
cat("================================================================\n\n")

cat("  FULL SAMPLE (TED-based analysis):\n")
cat(sprintf("    %s to %s\n", min(df_full$Date), max(df_full$Date)))
cat(sprintf("    Total rows: %d\n",  nrow(df_full)))
ted_valid <- sum(!is.na(df_full$TED_spread))
cat(sprintf("    Valid TED spread obs: %d\n", ted_valid))
cat(sprintf("    Valid Baa spread obs: %d\n", sum(!is.na(df_full$Baa_spread))))
cat(sprintf("    Valid Aaa spread obs: %d\n\n", sum(!is.na(df_full$Aaa_spread))))

cat("  POST-2018 SAMPLE (SOFR-EFFR-based VAR):\n")
cat(sprintf("    %s to %s\n", min(df_post$Date), max(df_post$Date)))
cat(sprintf("    Total rows: %d\n", nrow(df_post)))
cat(sprintf("    Valid SOFR-EFFR obs: %d\n", sum(!is.na(df_post$SOFR_EFFR))))
cat(sprintf("    Valid d_Reserves obs: %d\n", sum(!is.na(df_post$d_Reserves))))
cat(sprintf("    Valid d_TGA obs: %d\n", sum(!is.na(df_post$d_TGA))))
cat(sprintf("    Valid d_Baa obs: %d\n", sum(!is.na(df_post$d_Baa))))
cat(sprintf("    Valid d_Aaa obs: %d\n\n", sum(!is.na(df_post$d_Aaa))))

cat("  RISK APPETITE EXTENSION:\n")
hyg_lqd_valid <- sum(!is.na(df_post$d_log_HYG_LQD))
dxy_valid_post <- sum(!is.na(df_post$lr_DXY))
cat(sprintf("    Valid d_log(HYG/LQD) obs (post-2018): %d\n", hyg_lqd_valid))
cat(sprintf("    Valid DXY log-return obs (post-2018): %d\n\n", dxy_valid_post))


# ── 18. FINAL DIAGNOSTIC REPORT ─────────────────────────────────────────────

cat("================================================================\n")
cat("  CLEANING REPORT SUMMARY\n")
cat("================================================================\n\n")

cat(sprintf("  Input file:  raw_data_clean.csv (%d rows)\n", nrow(raw)))
cat(sprintf("  Output file: data_master.csv    (%d rows, %d columns)\n",
            nrow(df), ncol(df)))
cat(sprintf("  Holiday zeros removed:  %d rows\n", n_holidays))
cat(sprintf("  DXY zeros removed:      %d rows\n", n_dxy_zero))
cat("  Variables dropped:      SOFR_OBFR, BGCR_TGCR, M2, SOFR, EFFR, OBFR, TGCR, BGCR\n")
cat("  Weekly series forward-filled: Reserves, TGA, Fed_Assets, Liq_Swaps\n")
cat("  Constructed variables:  Baa_Aaa, d_Reserves, d_TGA, d_ON_RRP, d_Fed_Assets,\n")
cat("                          d_Baa, d_Aaa, d_Baa_Aaa, lr_HYG, lr_LQD, lr_DXY,\n")
cat("                          lr_TLT, lr_IEF, log_HYG_LQD, d_log_HYG_LQD\n\n")

cat("  LaTeX tables generated:\n")
cat("    table_data_sources.tex       — variable catalogue (Section 5.1)\n")
cat("    table_data_availability.tex  — date coverage per variable\n")
cat("    table_summary_post2018.tex   — summary statistics (post-2018)\n")
cat("    table_summary_full.tex       — summary statistics (full sample)\n")
cat("    table_adf_results.tex        — ADF unit root tests\n")
cat("    table_correlation_fd.tex     — correlation in first differences\n\n")

cat("  PDF figures generated:\n")
cat("    fig_ts_sofr_effr.pdf         — SOFR-EFFR spread\n")
cat("    fig_ts_credit_spreads.pdf    — Baa & Aaa spreads\n")
cat("    fig_ts_baa_aaa_diff.pdf      — quality differential\n")
cat("    fig_ts_reserves_tga.pdf      — reserves & TGA\n")
cat("    fig_ts_ted_spread.pdf        — TED spread (full sample)\n")
cat("    fig_ts_etf_prices.pdf        — HYG & LQD\n")
cat("    fig_ts_on_rrp.pdf            — ON RRP facility\n")
cat("    fig_ts_dxy.pdf               — DXY index\n")
cat("    fig_ts_funding_vs_credit.pdf — overlay plot\n")
cat("    fig_ts_vol_ratio.pdf         — Aaa/Baa volatility ratio\n")
cat("    fig_corr_fd_post2018.pdf     — correlation heat map (post-2018)\n")
cat("    fig_corr_fd_fullsample.pdf   — correlation heat map (full sample)\n\n")

cat("================================================================\n")
cat("  STEP 1 COMPLETE — Data ready for Granger tests (Step 2)\n")
cat("================================================================\n")
