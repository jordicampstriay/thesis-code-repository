###############################################################################
#  STEP 2 — GRANGER CAUSALITY TESTS  (Section 5.3)
#  Thesis: "The Plumbing of Liquidity"
#  Author: Jordi Camps Triay
#  Last updated: 2026-05-12
#
#  Prerequisites: data_master.csv from Step 1
#  Output: LaTeX tables, PDF figures, report_granger.pdf
#
#  METHOD (Granger, 1969):
#    The Granger test asks whether lagged values of X help predict Y
#    beyond Y's own lagged values. Rejection of H0 at p < 0.05 means X
#    contains predictive information about Y that is not redundant with
#    Y's own history. We test at lags 1, 5, and 10 business days.
#
#  DESIGN:
#    Battery 1 — Funding price → credit spreads (post-2018, actual SOFR)
#    Battery 2 — Funding quantities → credit spreads (post-2018)
#    Battery 3 — Reverse causality checks (post-2018)
#    Battery 4 — Extended sample with SOFR proxy (2003–2026)
#    Battery 5 — TED spread robustness (2003–2022)
#    Battery 6 — Risk appetite extension (post-2018)
###############################################################################

# ── 0. PACKAGES ──────────────────────────────────────────────────────────────
required_pkgs <- c("readr", "dplyr", "tidyr", "lmtest", "xtable", "ggplot2",
                   "scales", "zoo")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

source(file.path(dirname(if (interactive()) rstudioapi::getSourceEditorContext()$path else sys.frame(1)$ofile), "config.R"))

cat("================================================================\n")
cat("  STEP 2 — GRANGER CAUSALITY TESTS\n")
cat("================================================================\n\n")

# ── 1. LOAD DATA ─────────────────────────────────────────────────────────────
df <- read_csv(file.path(DATA_DIR, "data_master.csv"), show_col_types = FALSE)
df$Date <- as.Date(df$Date)

df_extended <- df %>% filter(Date >= as.Date("2003-05-01"))
df_post     <- df %>% filter(Date >= as.Date("2018-04-02"))
df_ted      <- df %>% filter(Date >= as.Date("2003-05-01"),
                              Date <= as.Date("2022-01-21"))

cat(sprintf("  Extended: %d rows | Post-2018: %d rows | TED: %d rows\n\n",
            nrow(df_extended), nrow(df_post), nrow(df_ted)))


# ── 2. GRANGER TEST FUNCTION ─────────────────────────────────────────────────
#
#  Returns a data.frame with one row per lag tested.
#  Uses lmtest::grangertest, which estimates:
#    Restricted: Y_t = a0 + a1*Y_{t-1} + ... + ap*Y_{t-p} + e_t
#    Unrestricted: Y_t = a0 + a1*Y_{t-1} + ... + ap*Y_{t-p}
#                       + b1*X_{t-1} + ... + bp*X_{t-p} + u_t
#  H0: b1 = b2 = ... = bp = 0  (X does not Granger-cause Y)
#  Test statistic: F-test on the joint restriction.

run_granger <- function(data, x_name, y_name, lags = c(1, 5, 10),
                         min_obs = 100) {
  xy <- data %>%
    select(all_of(c(x_name, y_name))) %>%
    drop_na()

  if (nrow(xy) < min_obs) {
    return(data.frame(
      X = x_name, Y = y_name, Lag = NA, N = nrow(xy),
      F_stat = NA, p_value = NA, Sig = "insufficient data",
      stringsAsFactors = FALSE))
  }

  results <- data.frame()
  for (lag in lags) {
    # Need at least lag + 30 obs for meaningful inference
    if (nrow(xy) < lag + 30) next
    tryCatch({
      gt <- grangertest(as.formula(paste(y_name, "~", x_name)),
                         order = lag, data = xy)
      f_val <- gt$F[2]
      p_val <- gt$`Pr(>F)`[2]
      sig <- ifelse(p_val < 0.01, "***",
             ifelse(p_val < 0.05, "**",
             ifelse(p_val < 0.10, "*", "")))
      results <- rbind(results, data.frame(
        X = x_name, Y = y_name, Lag = lag,
        N = nrow(xy),
        F_stat = round(f_val, 3),
        p_value = round(p_val, 4),
        Sig = sig,
        stringsAsFactors = FALSE))
    }, error = function(e) {
      cat(sprintf("  Error: %s -> %s lag %d: %s\n", x_name, y_name, lag, e$message))
    })
  }
  return(results)
}


###############################################################################
#  BATTERY 1 — FUNDING PRICE → CREDIT SPREADS  (Post-2018)
###############################################################################

cat("================================================================\n")
cat("  BATTERY 1: SOFR-EFFR → Credit Spreads (Post-2018)\n")
cat("================================================================\n\n")

# SOFR-EFFR is stationary in levels (ADF p = 0.01)
# Credit spread changes (d_Baa, d_Aaa, d_Baa_Aaa) are stationary

b1 <- data.frame()
for (y in c("d_Baa", "d_Aaa", "d_Baa_Aaa")) {
  b1 <- rbind(b1, run_granger(df_post, "SOFR_EFFR", y))
}
b1$Battery <- "1: Price -> Credit (Post-2018)"
cat("SOFR_EFFR → Credit spread changes:\n")
print(b1 %>% select(X, Y, Lag, N, F_stat, p_value, Sig), row.names = FALSE)
cat("\n")


###############################################################################
#  BATTERY 2 — FUNDING QUANTITIES → CREDIT SPREADS  (Post-2018)
###############################################################################

cat("================================================================\n")
cat("  BATTERY 2: Quantities → Credit Spreads (Post-2018)\n")
cat("================================================================\n\n")

# d_Reserves, d_TGA are stationary (first differences of unit-root series)
# d_ON_RRP also stationary in first differences

b2 <- data.frame()
for (x in c("d_Reserves", "d_TGA", "d_ON_RRP")) {
  for (y in c("d_Baa", "d_Aaa", "d_Baa_Aaa")) {
    b2 <- rbind(b2, run_granger(df_post, x, y))
  }
}
b2$Battery <- "2: Quantity -> Credit (Post-2018)"
cat("Quantity variables → Credit spread changes:\n")
print(b2 %>% select(X, Y, Lag, N, F_stat, p_value, Sig), row.names = FALSE)
cat("\n")


###############################################################################
#  BATTERY 3 — REVERSE CAUSALITY (Post-2018)
###############################################################################

cat("================================================================\n")
cat("  BATTERY 3: Reverse Causality Checks (Post-2018)\n")
cat("================================================================\n\n")

# Do credit spread changes predict funding conditions?
# If yes: endogeneity concern (Brunnermeier-Pedersen spiral).
# If no: strengthens directional interpretation.

b3 <- data.frame()
for (x in c("d_Baa", "d_Aaa")) {
  b3 <- rbind(b3, run_granger(df_post, x, "SOFR_EFFR"))
  b3 <- rbind(b3, run_granger(df_post, x, "d_Reserves"))
  b3 <- rbind(b3, run_granger(df_post, x, "d_TGA"))
}
b3$Battery <- "3: Reverse Causality (Post-2018)"
cat("Credit → Funding (reverse causality checks):\n")
print(b3 %>% select(X, Y, Lag, N, F_stat, p_value, Sig), row.names = FALSE)
cat("\n")


###############################################################################
#  BATTERY 4 — EXTENDED SAMPLE WITH SOFR PROXY (2003–2026)
###############################################################################

cat("================================================================\n")
cat("  BATTERY 4: Extended Sample (SOFR proxy, 2003-2026)\n")
cat("================================================================\n\n")

# SOFR-EFFR (with GC Repo proxy) stationary (ADF p = 0.01 on extended sample)

b4 <- data.frame()
# Funding price → credit
for (y in c("d_Baa", "d_Aaa", "d_Baa_Aaa")) {
  b4 <- rbind(b4, run_granger(df_extended, "SOFR_EFFR", y))
}
# Quantities → credit
for (x in c("d_Reserves", "d_TGA")) {
  for (y in c("d_Baa", "d_Aaa", "d_Baa_Aaa")) {
    b4 <- rbind(b4, run_granger(df_extended, x, y))
  }
}
# Reverse
for (x in c("d_Baa", "d_Aaa")) {
  b4 <- rbind(b4, run_granger(df_extended, x, "SOFR_EFFR"))
}
b4$Battery <- "4: Extended Sample (2003-2026)"
cat("Extended sample results:\n")
print(b4 %>% select(X, Y, Lag, N, F_stat, p_value, Sig), row.names = FALSE)
cat("\n")


###############################################################################
#  BATTERY 5 — TED SPREAD ROBUSTNESS (2003–2022)
###############################################################################

cat("================================================================\n")
cat("  BATTERY 5: TED Spread → Credit (2003-2022)\n")
cat("================================================================\n\n")

# TED spread is stationary in levels (ADF p = 0.01)

b5 <- data.frame()
for (y in c("d_Baa", "d_Aaa", "d_Baa_Aaa")) {
  b5 <- rbind(b5, run_granger(df_ted, "TED_spread", y))
}
# Reverse
for (x in c("d_Baa", "d_Aaa")) {
  b5 <- rbind(b5, run_granger(df_ted, x, "TED_spread"))
}
b5$Battery <- "5: TED Spread (2003-2022)"
cat("TED spread results:\n")
print(b5 %>% select(X, Y, Lag, N, F_stat, p_value, Sig), row.names = FALSE)
cat("\n")


###############################################################################
#  BATTERY 6 — RISK APPETITE EXTENSION (Post-2018)
###############################################################################

cat("================================================================\n")
cat("  BATTERY 6: Risk Appetite Extension (Post-2018)\n")
cat("================================================================\n\n")

b6 <- data.frame()
# Funding → HYG/LQD ratio
for (x in c("SOFR_EFFR", "d_Reserves", "d_TGA")) {
  b6 <- rbind(b6, run_granger(df_post, x, "d_log_HYG_LQD"))
}
# Funding → DXY
for (x in c("SOFR_EFFR", "d_Reserves")) {
  b6 <- rbind(b6, run_granger(df_post, x, "lr_DXY"))
}
# Funding → EMB
b6 <- rbind(b6, run_granger(df_post, "d_Reserves", "lr_EMB"))
b6 <- rbind(b6, run_granger(df_post, "SOFR_EFFR", "lr_EMB"))

b6$Battery <- "6: Risk Appetite (Post-2018)"
cat("Risk appetite results:\n")
print(b6 %>% select(X, Y, Lag, N, F_stat, p_value, Sig), row.names = FALSE)
cat("\n")


###############################################################################
#  COMBINED RESULTS AND INTERPRETATION
###############################################################################

cat("================================================================\n")
cat("  COMBINED RESULTS\n")
cat("================================================================\n\n")

all_results <- bind_rows(b1, b2, b3, b4, b5, b6)

# Significant results
sig_results <- all_results %>% filter(p_value < 0.05)
cat("=== SIGNIFICANT AT 5% (p < 0.05) ===\n\n")
if (nrow(sig_results) > 0) {
  print(sig_results %>%
          select(Battery, X, Y, Lag, N, F_stat, p_value, Sig) %>%
          arrange(Battery, X, Y, Lag),
        row.names = FALSE)
} else {
  cat("  No significant results at 5% level.\n")
}
cat("\n")

# Marginally significant (5-10%)
marginal <- all_results %>% filter(p_value >= 0.05 & p_value < 0.10)
cat("=== MARGINALLY SIGNIFICANT (0.05 <= p < 0.10) ===\n\n")
if (nrow(marginal) > 0) {
  print(marginal %>%
          select(Battery, X, Y, Lag, N, F_stat, p_value, Sig) %>%
          arrange(Battery, X, Y, Lag),
        row.names = FALSE)
} else {
  cat("  None.\n")
}
cat("\n")


###############################################################################
#  LATEX TABLES
###############################################################################

cat("================================================================\n")
cat("  GENERATING LATEX TABLES\n")
cat("================================================================\n\n")

# ── Table A: Post-2018 primary results (Batteries 1 + 2) ──
primary <- bind_rows(b1, b2) %>%
  select(X, Y, Lag, N, F_stat, p_value, Sig) %>%
  mutate(Direction = paste0(X, " $\\rightarrow$ ", Y)) %>%
  select(Direction, Lag, N, F_stat, p_value, Sig)

colnames(primary) <- c("Direction", "Lag", "$N$", "$F$-stat", "$p$-value", "")
primary_xt <- xtable(primary,
  caption = "Granger Causality Tests: Funding Conditions $\\rightarrow$ Credit Spread Changes (Post-2018 Sample)",
  label = "tab:granger_post2018",
  digits = c(0, 0, 0, 0, 3, 4, 0))
print(primary_xt, file = file.path(TBL_DIR, "table_granger_post2018.tex"),
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78)
cat("  Saved: table_granger_post2018.tex\n")

# ── Table B: Reverse causality (Battery 3) ──
reverse <- b3 %>%
  select(X, Y, Lag, N, F_stat, p_value, Sig) %>%
  mutate(Direction = paste0(X, " $\\rightarrow$ ", Y)) %>%
  select(Direction, Lag, N, F_stat, p_value, Sig)

colnames(reverse) <- c("Direction", "Lag", "$N$", "$F$-stat", "$p$-value", "")
reverse_xt <- xtable(reverse,
  caption = "Reverse Causality Tests: Credit Spread Changes $\\rightarrow$ Funding Conditions (Post-2018)",
  label = "tab:granger_reverse",
  digits = c(0, 0, 0, 0, 3, 4, 0))
print(reverse_xt, file = file.path(TBL_DIR, "table_granger_reverse.tex"),
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78)
cat("  Saved: table_granger_reverse.tex\n")

# ── Table C: Extended sample (Battery 4) ──
ext_tbl <- b4 %>%
  select(X, Y, Lag, N, F_stat, p_value, Sig) %>%
  mutate(Direction = paste0(X, " $\\rightarrow$ ", Y)) %>%
  select(Direction, Lag, N, F_stat, p_value, Sig)

colnames(ext_tbl) <- c("Direction", "Lag", "$N$", "$F$-stat", "$p$-value", "")
ext_xt <- xtable(ext_tbl,
  caption = "Granger Causality Tests: Extended Sample with SOFR Proxy (2003--2026)",
  label = "tab:granger_extended",
  digits = c(0, 0, 0, 0, 3, 4, 0))
print(ext_xt, file = file.path(TBL_DIR, "table_granger_extended.tex"),
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78)
cat("  Saved: table_granger_extended.tex\n")

# ── Table D: TED spread robustness (Battery 5) ──
ted_tbl <- b5 %>%
  select(X, Y, Lag, N, F_stat, p_value, Sig) %>%
  mutate(Direction = paste0(X, " $\\rightarrow$ ", Y)) %>%
  select(Direction, Lag, N, F_stat, p_value, Sig)

colnames(ted_tbl) <- c("Direction", "Lag", "$N$", "$F$-stat", "$p$-value", "")
ted_xt <- xtable(ted_tbl,
  caption = "Granger Causality Tests: TED Spread (2003--2022, Robustness)",
  label = "tab:granger_ted",
  digits = c(0, 0, 0, 0, 3, 4, 0))
print(ted_xt, file = file.path(TBL_DIR, "table_granger_ted.tex"),
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78)
cat("  Saved: table_granger_ted.tex\n")

# ── Table E: Significant results summary ──
if (nrow(sig_results) > 0) {
  sig_tbl <- sig_results %>%
    mutate(Direction = paste0(X, " $\\rightarrow$ ", Y),
           Sample = gsub("^[0-9]+: ", "", Battery)) %>%
    select(Sample, Direction, Lag, N, F_stat, p_value, Sig) %>%
    arrange(Sample, Direction, Lag)

  colnames(sig_tbl) <- c("Sample", "Direction", "Lag", "$N$",
                          "$F$-stat", "$p$-value", "")
  sig_xt <- xtable(sig_tbl,
    caption = "Significant Granger Causality Results ($p < 0.05$)",
    label = "tab:granger_significant",
    digits = c(0, 0, 0, 0, 0, 3, 4, 0))
  print(sig_xt, file = file.path(TBL_DIR, "table_granger_significant.tex"),
        include.rownames = FALSE, booktabs = TRUE,
        caption.placement = "top", table.placement = "htbp",
        sanitize.text.function = identity,
        sanitize.colnames.function = identity,
        scalebox = 0.75)
  cat("  Saved: table_granger_significant.tex\n")
}

# ── Table F: Risk appetite (Battery 6) ──
ra_tbl <- b6 %>%
  select(X, Y, Lag, N, F_stat, p_value, Sig) %>%
  mutate(Direction = paste0(X, " $\\rightarrow$ ", Y)) %>%
  select(Direction, Lag, N, F_stat, p_value, Sig)

colnames(ra_tbl) <- c("Direction", "Lag", "$N$", "$F$-stat", "$p$-value", "")
ra_xt <- xtable(ra_tbl,
  caption = "Granger Causality Tests: Funding $\\rightarrow$ Risk Appetite (Post-2018)",
  label = "tab:granger_risk_appetite",
  digits = c(0, 0, 0, 0, 3, 4, 0))
print(ra_xt, file = file.path(TBL_DIR, "table_granger_risk_appetite.tex"),
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78)
cat("  Saved: table_granger_risk_appetite.tex\n")


###############################################################################
#  SAVE COMPLETE RESULTS FOR DOWNSTREAM USE
###############################################################################

write.csv(all_results, "granger_all_results.csv", row.names = FALSE)
cat("\n  Saved: granger_all_results.csv\n\n")


###############################################################################
#  INTERPRETATION SUMMARY
###############################################################################

cat("================================================================\n")
cat("  INTERPRETATION SUMMARY\n")
cat("================================================================\n\n")

# Count significant by battery
for (bat in unique(all_results$Battery)) {
  sub <- all_results %>% filter(Battery == bat)
  n_sig <- sum(sub$p_value < 0.05, na.rm = TRUE)
  n_tot <- sum(!is.na(sub$p_value))
  cat(sprintf("  %s: %d/%d significant at 5%%\n", bat, n_sig, n_tot))
}

cat("\n")
cat("  KEY INTERPRETATIONS:\n\n")

cat("  1. QUANTITY CHANNEL (Reserves → Credit):\n")
b2_res_aaa <- b2 %>% filter(X == "d_Reserves", Y == "d_Aaa")
b2_res_baa <- b2 %>% filter(X == "d_Reserves", Y == "d_Baa")
b2_res_diff <- b2 %>% filter(X == "d_Reserves", Y == "d_Baa_Aaa")
cat(sprintf("     d_Reserves → d_Aaa: p = %s at lag 1, %s at lag 5, %s at lag 10\n",
            paste(b2_res_aaa$p_value, collapse=", "),
            ifelse(any(b2_res_aaa$p_value < 0.05), "SIGNIFICANT", "not sig"),
            ""))
cat(sprintf("     d_Reserves → d_Baa: p = %s\n",
            paste(b2_res_baa$p_value, collapse=", ")))
cat(sprintf("     d_Reserves → d_Baa_Aaa: p = %s\n\n",
            paste(b2_res_diff$p_value, collapse=", ")))

cat("  2. PRICE CHANNEL (SOFR-EFFR → Credit, post-2018):\n")
b1_baa <- b1 %>% filter(Y == "d_Baa")
b1_aaa <- b1 %>% filter(Y == "d_Aaa")
cat(sprintf("     SOFR_EFFR → d_Baa: p = %s\n",
            paste(b1_baa$p_value, collapse=", ")))
cat(sprintf("     SOFR_EFFR → d_Aaa: p = %s\n\n",
            paste(b1_aaa$p_value, collapse=", ")))

cat("  3. TGA EXOGENEITY ADVANTAGE:\n")
b2_tga <- b2 %>% filter(X == "d_TGA")
cat(sprintf("     d_TGA → d_Baa_Aaa: p = %s\n",
            paste((b2_tga %>% filter(Y == "d_Baa_Aaa"))$p_value, collapse=", ")))
cat(sprintf("     d_TGA → d_Aaa: p = %s\n\n",
            paste((b2_tga %>% filter(Y == "d_Aaa"))$p_value, collapse=", ")))

cat("  4. TED SPREAD (full sample):\n")
b5_fwd <- b5 %>% filter(X == "TED_spread")
cat(sprintf("     TED → d_Baa: p = %s\n",
            paste((b5_fwd %>% filter(Y == "d_Baa"))$p_value, collapse=", ")))
cat(sprintf("     TED → d_Aaa: p = %s\n\n",
            paste((b5_fwd %>% filter(Y == "d_Aaa"))$p_value, collapse=", ")))

cat("  5. REVERSE CAUSALITY:\n")
b3_rev_sofr <- b3 %>% filter(Y == "SOFR_EFFR")
cat(sprintf("     d_Baa → SOFR_EFFR: p = %s\n",
            paste((b3_rev_sofr %>% filter(X == "d_Baa"))$p_value, collapse=", ")))
cat(sprintf("     d_Aaa → SOFR_EFFR: p = %s\n\n",
            paste((b3_rev_sofr %>% filter(X == "d_Aaa"))$p_value, collapse=", ")))

cat("  6. EXTENDED SAMPLE (with proxy):\n")
b4_sofr <- b4 %>% filter(X == "SOFR_EFFR")
cat(sprintf("     SOFR_EFFR(proxy) → d_Baa: p = %s\n",
            paste((b4_sofr %>% filter(Y == "d_Baa"))$p_value, collapse=", ")))
cat(sprintf("     SOFR_EFFR(proxy) → d_Aaa: p = %s\n",
            paste((b4_sofr %>% filter(Y == "d_Aaa"))$p_value, collapse=", ")))
cat(sprintf("     SOFR_EFFR(proxy) → d_Baa_Aaa: p = %s\n\n",
            paste((b4_sofr %>% filter(Y == "d_Baa_Aaa"))$p_value, collapse=", ")))

cat("================================================================\n")
cat("  STEP 2 COMPLETE\n")
cat("================================================================\n")
