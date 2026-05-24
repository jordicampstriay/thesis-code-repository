###############################################################################
#  STEP 2 (ENHANCED) — GRANGER CAUSALITY WITH STRUCTURAL BREAK ANALYSIS
#  Thesis: "Beyond Interest Rates: Liquidity Mechanics and Risk Appetite
#           in the Financial System"
#  Author: Jordi Camps Triay
#  Date: 2026-05-17
#
#  STRUCTURE:
#    Part A — Post-2018 Sample (SOFR-EFFR, Reserves, TGA, ON RRP)
#      A1. Summary statistics & ADF
#      A2. Granger causality (forward + reverse)
#      A3. Chow test at midpoint + SupF endogenous break
#      A4. Time-series and diagnostic plots
#
#    Part B — TED Spread Sample (2003-2022)
#      B1. Summary statistics & ADF
#      B2. Granger causality (forward + reverse)
#      B3. Structural break at Jan 2015 (Chow + SupF)
#      B4. Subsample Granger: pre-2015 vs post-2015
#      B5. Time-series and diagnostic plots
#
#    Part C — Comparison and synthesis tables
#
#  INPUT:  data_master.csv (from Step 1)
#  OUTPUT: LaTeX tables, PDF figures, CSV results
#
#  NOTES:
#  - The 2003-2026 extended sample using SOFR proxy is DROPPED from the
#    thesis. The proxy (PD Survey Rate) differs from actual SOFR in
#    construction (mean vs median), coverage (GC-only vs bilateral), and
#    exhibits a structural break at 2015 (FRBNY FEDS Note, Anbil et al. 2019).
#  - The TED spread is used as an independent robustness check, explicitly
#    acknowledging it measures unsecured (LIBOR-based) funding costs, not
#    secured (repo) funding.
#  - Zero-value observations in TED spread (89 days) and SOFR-EFFR (228 days)
#    are genuine data points (spread = 0 when rates coincide), not missing
#    values. They are retained in the analysis.
#  - Forward-filled weekly series (Reserves, TGA) introduce measurement error;
#    zeros in d_Reserves/d_TGA on non-release days are an artifact of this.
###############################################################################

# ── 0. PACKAGES ─────────────────────────────────────────────────────────────
required_pkgs <- c("readr", "dplyr", "tidyr", "lmtest", "tseries", "strucchange",
                   "moments", "xtable", "ggplot2", "scales", "zoo", "gridExtra")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

WD <- "/Users/jordi/Downloads/University/TFG/Data/Thesis/step2_granger_enhanced"
setwd(WD)

# Load from parent directory
df <- read_csv("../data_master.csv", show_col_types = FALSE)
df$Date <- as.Date(df$Date)

cat("================================================================\n")
cat("  STEP 2 (ENHANCED) — GRANGER CAUSALITY + STRUCTURAL BREAKS\n")
cat("================================================================\n\n")

# ── Plot theme ──
theme_thesis <- theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


###############################################################################
#  HELPER FUNCTIONS
###############################################################################

run_granger <- function(data, x_name, y_name, lags = c(1, 5, 10),
                         min_obs = 100) {
  xy <- data %>%
    select(all_of(c(x_name, y_name))) %>%
    tidyr::drop_na()

  if (nrow(xy) < min_obs) {
    return(data.frame(
      X = x_name, Y = y_name, Lag = NA, N = nrow(xy),
      F_stat = NA, p_value = NA, Sig = "insuff.",
      stringsAsFactors = FALSE))
  }

  results <- data.frame()
  for (lag in lags) {
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
        N = nrow(xy), F_stat = round(f_val, 3),
        p_value = round(p_val, 4), Sig = sig,
        stringsAsFactors = FALSE))
    }, error = function(e) {
      cat(sprintf("  Error: %s -> %s lag %d: %s\n", x_name, y_name, lag, e$message))
    })
  }
  return(results)
}

run_adf <- function(x, varname, period) {
  x <- x[!is.na(x)]
  if (length(x) < 50) return(data.frame(Period=period, Variable=varname,
                                          ADF_stat=NA, p_value=NA, N=length(x),
                                          Stationary="insuff."))
  test <- adf.test(x, alternative = "stationary")
  data.frame(
    Period   = period,
    Variable = varname,
    ADF_stat = round(test$statistic, 3),
    p_value  = round(test$p.value, 4),
    N        = length(x),
    Stationary = ifelse(test$p.value < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}

compute_stats <- function(data, vars, label) {
  stats <- data.frame()
  for (v in vars) {
    x <- data[[v]]
    x <- x[!is.na(x)]
    if (length(x) < 10) next
    stats <- rbind(stats, data.frame(
      Period   = label,
      Variable = v,
      N        = length(x),
      Mean     = round(mean(x), 4),
      SD       = round(sd(x), 4),
      Median   = round(median(x), 4),
      Min      = round(min(x), 4),
      Max      = round(max(x), 4),
      Skewness = round(skewness(x), 2),
      Kurtosis = round(kurtosis(x), 2),
      stringsAsFactors = FALSE
    ))
  }
  return(stats)
}

run_chow <- function(data, x_name, y_name, lag, break_point_idx) {
  xy <- data %>%
    select(Date, all_of(c(x_name, y_name))) %>%
    tidyr::drop_na() %>%
    arrange(Date)

  n <- nrow(xy)
  if (n < lag + 50) return(NULL)

  Y <- xy[[y_name]][(lag+1):n]
  X_mat <- data.frame(row.names = 1:length(Y))
  for (i in 1:lag) {
    X_mat[[paste0("Y_lag", i)]] <- xy[[y_name]][(lag+1-i):(n-i)]
    X_mat[[paste0("X_lag", i)]] <- xy[[x_name]][(lag+1-i):(n-i)]
  }
  X_mat$Y <- Y
  dates_used <- xy$Date[(lag+1):n]

  # Find index closest to break point
  if (is.numeric(break_point_idx)) {
    bp <- break_point_idx
  } else {
    bp <- which.min(abs(dates_used - as.Date(break_point_idx)))
  }

  tryCatch({
    sc <- sctest(Y ~ ., data = X_mat, type = "Chow", point = bp)
    data.frame(
      Direction = paste0(x_name, " -> ", y_name),
      Lag       = lag,
      Break_date = as.character(dates_used[bp]),
      Chow_F    = round(sc$statistic, 3),
      p_value   = round(sc$p.value, 4),
      Sig       = ifelse(sc$p.value < 0.01, "***",
                  ifelse(sc$p.value < 0.05, "**",
                  ifelse(sc$p.value < 0.10, "*", ""))),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat(sprintf("  Chow error: %s -> %s lag %d: %s\n", x_name, y_name, lag, e$message))
    NULL
  })
}

run_supf <- function(data, x_name, y_name, lag) {
  xy <- data %>%
    select(Date, all_of(c(x_name, y_name))) %>%
    tidyr::drop_na() %>%
    arrange(Date)

  n <- nrow(xy)
  if (n < lag + 50) return(NULL)

  Y <- xy[[y_name]][(lag+1):n]
  X_mat <- data.frame(row.names = 1:length(Y))
  for (i in 1:lag) {
    X_mat[[paste0("Y_lag", i)]] <- xy[[y_name]][(lag+1-i):(n-i)]
    X_mat[[paste0("X_lag", i)]] <- xy[[x_name]][(lag+1-i):(n-i)]
  }
  X_mat$Y <- Y

  tryCatch({
    sc <- sctest(Y ~ ., data = X_mat, type = "supF")
    data.frame(
      Direction = paste0(x_name, " -> ", y_name),
      Lag       = lag,
      SupF      = round(sc$statistic, 3),
      p_value   = round(sc$p.value, 4),
      Sig       = ifelse(sc$p.value < 0.01, "***",
                  ifelse(sc$p.value < 0.05, "**",
                  ifelse(sc$p.value < 0.10, "*", ""))),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat(sprintf("  SupF error: %s -> %s lag %d: %s\n", x_name, y_name, lag, e$message))
    NULL
  })
}

# Rolling Granger p-value
roll_granger <- function(data, x_name, y_name, lag, window = 500, step = 20) {
  xy <- data %>%
    select(Date, all_of(c(x_name, y_name))) %>%
    tidyr::drop_na() %>%
    arrange(Date)

  n <- nrow(xy)
  if (n < window + lag) return(data.frame(Date=as.Date(character()), p_value=numeric()))

  out <- data.frame()
  for (i in seq(window, n, by = step)) {
    sub <- xy[(i - window + 1):i, ]
    tryCatch({
      gt <- grangertest(as.formula(paste(y_name, "~", x_name)),
                         order = lag, data = sub)
      out <- rbind(out, data.frame(
        Date = sub$Date[nrow(sub)],
        p_value = gt$`Pr(>F)`[2]
      ))
    }, error = function(e) {})
  }
  return(out)
}


###############################################################################
###############################################################################
#                     PART A — POST-2018 SAMPLE (SOFR-EFFR)
###############################################################################
###############################################################################

cat("\n================================================================\n")
cat("  PART A: POST-2018 SAMPLE (SOFR-EFFR)\n")
cat("================================================================\n\n")

df_post <- df %>% filter(Date >= as.Date("2018-04-02"))

cat(sprintf("  Sample: %s to %s (%d business days)\n",
            min(df_post$Date), max(df_post$Date), nrow(df_post)))
cat(sprintf("  Valid SOFR-EFFR: %d | NA: %d | Zeros: %d\n\n",
            sum(!is.na(df_post$SOFR_EFFR)),
            sum(is.na(df_post$SOFR_EFFR)),
            sum(df_post$SOFR_EFFR == 0, na.rm = TRUE)))


# ── A1. Summary Statistics ──

a_vars <- c("SOFR_EFFR", "d_Reserves", "d_TGA", "d_ON_RRP",
            "d_Baa", "d_Aaa", "d_Baa_Aaa")

stats_a <- compute_stats(df_post, a_vars, "Post-2018")

cat("Summary statistics:\n")
print(stats_a %>% select(Variable, N, Mean, SD, Skewness, Kurtosis), row.names = FALSE)
cat("\n")

# ADF tests
adf_a <- data.frame()
for (v in a_vars) {
  adf_a <- rbind(adf_a, run_adf(df_post[[v]], v, "Post-2018"))
}
cat("ADF unit root tests:\n")
print(adf_a, row.names = FALSE)
cat("\n")


# ── A2. Granger Causality ──

cat("================================================================\n")
cat("  A2: GRANGER CAUSALITY (POST-2018)\n")
cat("================================================================\n\n")

# Forward: funding → credit
fwd_pairs_a <- list(
  c("SOFR_EFFR", "d_Baa"), c("SOFR_EFFR", "d_Aaa"), c("SOFR_EFFR", "d_Baa_Aaa"),
  c("d_Reserves", "d_Baa"), c("d_Reserves", "d_Aaa"), c("d_Reserves", "d_Baa_Aaa"),
  c("d_TGA", "d_Baa"), c("d_TGA", "d_Aaa"), c("d_TGA", "d_Baa_Aaa"),
  c("d_ON_RRP", "d_Baa"), c("d_ON_RRP", "d_Aaa"), c("d_ON_RRP", "d_Baa_Aaa")
)

# Reverse: credit → funding
rev_pairs_a <- list(
  c("d_Baa", "SOFR_EFFR"), c("d_Aaa", "SOFR_EFFR"),
  c("d_Baa", "d_Reserves"), c("d_Aaa", "d_Reserves"),
  c("d_Baa", "d_TGA"), c("d_Aaa", "d_TGA")
)

gc_fwd_a <- data.frame()
for (pair in fwd_pairs_a) {
  gc_fwd_a <- rbind(gc_fwd_a, run_granger(df_post, pair[1], pair[2]))
}
gc_fwd_a$Type <- "Forward"

gc_rev_a <- data.frame()
for (pair in rev_pairs_a) {
  gc_rev_a <- rbind(gc_rev_a, run_granger(df_post, pair[1], pair[2]))
}
gc_rev_a$Type <- "Reverse"

gc_a <- bind_rows(gc_fwd_a, gc_rev_a)

cat("FORWARD: Funding → Credit Spreads\n")
print(gc_fwd_a %>% select(X, Y, Lag, N, F_stat, p_value, Sig), row.names = FALSE)
cat("\nREVERSE: Credit Spreads → Funding\n")
print(gc_rev_a %>% select(X, Y, Lag, N, F_stat, p_value, Sig), row.names = FALSE)
cat("\n")


# ── A3. Structural Break Tests (Post-2018) ──

cat("================================================================\n")
cat("  A3: STRUCTURAL BREAK TESTS (POST-2018)\n")
cat("================================================================\n\n")

# Chow at midpoint
chow_a <- data.frame()
all_pairs_a <- c(fwd_pairs_a, rev_pairs_a)
for (pair in all_pairs_a) {
  for (lag in c(1, 5, 10)) {
    xy <- df_post %>% select(Date, all_of(c(pair[1], pair[2]))) %>%
      tidyr::drop_na() %>% arrange(Date)
    mid <- floor(nrow(xy) / 2)
    res <- run_chow(df_post, pair[1], pair[2], lag, mid)
    if (!is.null(res)) chow_a <- rbind(chow_a, res)
  }
}

cat("Chow test at sample midpoint (~Apr 2022):\n\n")
print(chow_a, row.names = FALSE)
cat(sprintf("\nSignificant at 5%%: %d / %d\n\n", sum(chow_a$p_value < 0.05), nrow(chow_a)))

# Chow at March 2020 (COVID)
cat("Chow test at March 2020 (COVID crisis):\n\n")
chow_a_covid <- data.frame()
for (pair in all_pairs_a) {
  for (lag in c(1, 5, 10)) {
    res <- run_chow(df_post, pair[1], pair[2], lag, "2020-03-15")
    if (!is.null(res)) chow_a_covid <- rbind(chow_a_covid, res)
  }
}
print(chow_a_covid, row.names = FALSE)
cat(sprintf("\nSignificant at 5%%: %d / %d\n\n", sum(chow_a_covid$p_value < 0.05), nrow(chow_a_covid)))

# SupF (endogenous break)
cat("SupF test (endogenous break detection):\n\n")
supf_a <- data.frame()
for (pair in fwd_pairs_a) {
  for (lag in c(1, 5, 10)) {
    res <- run_supf(df_post, pair[1], pair[2], lag)
    if (!is.null(res)) supf_a <- rbind(supf_a, res)
  }
}
print(supf_a, row.names = FALSE)
cat(sprintf("\nSupF significant at 5%%: %d / %d\n\n", sum(supf_a$p_value < 0.05), nrow(supf_a)))


# ── A4. Figures (Post-2018) ──

cat("================================================================\n")
cat("  A4: FIGURES (POST-2018)\n")
cat("================================================================\n\n")

# Clean plot data: remove NA but keep zeros (they are real)
df_plot_a <- df_post %>% filter(!is.na(SOFR_EFFR))

# Fig A1: SOFR-EFFR time series
p_a1 <- ggplot(df_plot_a, aes(x = Date, y = SOFR_EFFR * 100)) +
  geom_line(color = "steelblue", linewidth = 0.4) +
  geom_vline(xintercept = as.Date("2020-03-15"), linetype = "dashed",
             color = "red", linewidth = 0.6) +
  annotate("text", x = as.Date("2020-06-01"), y = max(df_plot_a$SOFR_EFFR*100, na.rm=T)*0.85,
           label = "COVID\n(Mar 2020)", color = "red", size = 3, hjust = 0) +
  labs(title = "SOFR-EFFR Spread (Post-2018)",
       subtitle = "Overnight secured-unsecured funding spread; zero = no arbitrage gap",
       x = NULL, y = "Basis points") +
  theme_thesis
ggsave("fig_a1_sofr_effr.pdf", p_a1, width = 10, height = 5)
cat("  Saved: fig_a1_sofr_effr.pdf\n")

# Fig A2: Credit spreads (clean — remove NA days, not zero-spread days)
df_credit_a <- df_post %>%
  filter(!is.na(Baa_spread) & !is.na(Aaa_spread)) %>%
  select(Date, Baa_spread, Aaa_spread) %>%
  tidyr::pivot_longer(-Date, names_to = "Variable", values_to = "Spread") %>%
  mutate(Variable = ifelse(Variable == "Baa_spread", "Baa", "Aaa"))

p_a2 <- ggplot(df_credit_a, aes(x = Date, y = Spread, color = Variable)) +
  geom_line(linewidth = 0.4) +
  geom_vline(xintercept = as.Date("2020-03-15"), linetype = "dashed",
             color = "red", linewidth = 0.6) +
  scale_color_manual(values = c("Baa" = "darkorange", "Aaa" = "steelblue")) +
  labs(title = "Moody's Corporate Credit Spreads (Post-2018)",
       subtitle = "Option-adjusted spreads over risk-free rate",
       x = NULL, y = "Spread (%)", color = NULL) +
  theme_thesis
ggsave("fig_a2_credit_spreads.pdf", p_a2, width = 10, height = 5)
cat("  Saved: fig_a2_credit_spreads.pdf\n")

# Fig A3: Reserves & TGA (clean)
df_qty_a <- df_post %>%
  filter(!is.na(Reserves)) %>%
  select(Date, Reserves, TGA) %>%
  mutate(Reserves = Reserves / 1e6, TGA = TGA / 1e6) %>%
  tidyr::pivot_longer(-Date, names_to = "Variable", values_to = "Trillions")

p_a3 <- ggplot(df_qty_a, aes(x = Date, y = Trillions, color = Variable)) +
  geom_line(linewidth = 0.5) +
  geom_vline(xintercept = as.Date("2020-03-15"), linetype = "dashed",
             color = "red", linewidth = 0.6) +
  scale_color_manual(values = c("Reserves" = "darkgreen", "TGA" = "purple")) +
  labs(title = "Federal Reserve Reserves and Treasury General Account (Post-2018)",
       subtitle = "Weekly H.4.1 data, forward-filled to daily",
       x = NULL, y = "Trillions USD", color = NULL) +
  theme_thesis
ggsave("fig_a3_reserves_tga.pdf", p_a3, width = 10, height = 5)
cat("  Saved: fig_a3_reserves_tga.pdf\n")

# Fig A4: ON RRP (clean)
df_onrrp <- df_post %>%
  filter(!is.na(ON_RRP))

p_a4 <- ggplot(df_onrrp, aes(x = Date, y = ON_RRP)) +
  geom_line(color = "darkred", linewidth = 0.5) +
  labs(title = "ON RRP Outstanding (Post-2018)",
       subtitle = "Overnight reverse repo operations, FRED RRPONTSYD (billions USD)",
       x = NULL, y = "Billions USD") +
  theme_thesis
ggsave("fig_a4_on_rrp.pdf", p_a4, width = 10, height = 5)
cat("  Saved: fig_a4_on_rrp.pdf\n")

# Fig A5: Rolling Granger p-values (Post-2018)
cat("  Computing rolling Granger p-values (post-2018)...\n")

roll_a_baa <- roll_granger(df_post, "SOFR_EFFR", "d_Baa", lag = 5, window = 400, step = 15)
roll_a_baa$Direction <- "SOFR-EFFR -> d_Baa (lag 5)"
roll_a_aaa <- roll_granger(df_post, "SOFR_EFFR", "d_Aaa", lag = 5, window = 400, step = 15)
roll_a_aaa$Direction <- "SOFR-EFFR -> d_Aaa (lag 5)"

roll_a <- bind_rows(roll_a_baa, roll_a_aaa)

if (nrow(roll_a) > 0) {
  p_a5 <- ggplot(roll_a, aes(x = Date, y = p_value, color = Direction)) +
    geom_line(linewidth = 0.6) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
    geom_vline(xintercept = as.Date("2020-03-15"), linetype = "dashed",
               color = "red", linewidth = 0.6) +
    annotate("text", x = as.Date("2020-06-01"), y = 0.95,
             label = "COVID", color = "red", size = 3, hjust = 0) +
    scale_color_manual(values = c("SOFR-EFFR -> d_Baa (lag 5)" = "darkorange",
                                   "SOFR-EFFR -> d_Aaa (lag 5)" = "steelblue")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Rolling Granger p-values: SOFR-EFFR -> Credit Spreads (Post-2018)",
         subtitle = "400-day rolling window, lag 5; black dashed = 5% significance",
         x = NULL, y = "p-value", color = NULL) +
    theme_thesis
  ggsave("fig_a5_rolling_granger.pdf", p_a5, width = 10, height = 5)
  cat("  Saved: fig_a5_rolling_granger.pdf\n")
}


###############################################################################
###############################################################################
#                     PART B — TED SPREAD SAMPLE (2003-2022)
###############################################################################
###############################################################################

cat("\n\n================================================================\n")
cat("  PART B: TED SPREAD SAMPLE (2003-2022)\n")
cat("================================================================\n\n")

BREAK_DATE_B <- as.Date("2015-01-01")

df_ted <- df %>% filter(!is.na(TED_spread), Date >= as.Date("2003-05-01"))
df_pre  <- df_ted %>% filter(Date < BREAK_DATE_B)
df_post_b <- df_ted %>% filter(Date >= BREAK_DATE_B)

cat(sprintf("  Full TED sample:  %s to %s (%d obs)\n",
            min(df_ted$Date), max(df_ted$Date), nrow(df_ted)))
cat(sprintf("  Pre-2015:         %s to %s (%d obs)\n",
            min(df_pre$Date), max(df_pre$Date), nrow(df_pre)))
cat(sprintf("  Post-2015:        %s to %s (%d obs)\n\n",
            min(df_post_b$Date), max(df_post_b$Date), nrow(df_post_b)))


# ── B1. Summary Statistics ──

b_vars <- c("TED_spread", "d_Reserves", "d_TGA", "d_Baa", "d_Aaa", "d_Baa_Aaa")

stats_b_full <- compute_stats(df_ted, b_vars, "Full (2003-2022)")
stats_b_pre  <- compute_stats(df_pre, b_vars, "Pre-2015")
stats_b_post <- compute_stats(df_post_b, b_vars, "Post-2015")
stats_b <- bind_rows(stats_b_full, stats_b_pre, stats_b_post)

cat("Summary statistics:\n")
print(stats_b %>% select(Period, Variable, N, Mean, SD, Skewness, Kurtosis),
      row.names = FALSE)
cat("\n")

# ADF tests
adf_b <- data.frame()
for (v in b_vars) {
  adf_b <- rbind(adf_b,
                  run_adf(df_ted[[v]], v, "Full"),
                  run_adf(df_pre[[v]], v, "Pre-2015"),
                  run_adf(df_post_b[[v]], v, "Post-2015"))
}
cat("ADF results:\n")
print(adf_b, row.names = FALSE)
cat("\n")


# ── B2. Granger Causality (Full + Subsamples) ──

cat("================================================================\n")
cat("  B2: GRANGER CAUSALITY (TED SPREAD)\n")
cat("================================================================\n\n")

fwd_pairs_b <- list(
  c("TED_spread", "d_Baa"), c("TED_spread", "d_Aaa"), c("TED_spread", "d_Baa_Aaa"),
  c("d_Reserves", "d_Baa"), c("d_Reserves", "d_Aaa"), c("d_Reserves", "d_Baa_Aaa"),
  c("d_TGA", "d_Baa"), c("d_TGA", "d_Aaa"), c("d_TGA", "d_Baa_Aaa")
)

rev_pairs_b <- list(
  c("d_Baa", "TED_spread"), c("d_Aaa", "TED_spread"),
  c("d_Baa", "d_Reserves"), c("d_Aaa", "d_Reserves"),
  c("d_Baa", "d_TGA"), c("d_Aaa", "d_TGA")
)

all_pairs_b <- c(fwd_pairs_b, rev_pairs_b)

# Run on three samples
run_all_pairs <- function(data, pairs, label) {
  res <- data.frame()
  for (pair in pairs) {
    res <- rbind(res, run_granger(data, pair[1], pair[2]))
  }
  res$Period <- label
  return(res)
}

gc_b_full <- run_all_pairs(df_ted, all_pairs_b, "Full (2003-2022)")
gc_b_pre  <- run_all_pairs(df_pre, all_pairs_b, "Pre-2015 (2003-2014)")
gc_b_post <- run_all_pairs(df_post_b, all_pairs_b, "Post-2015 (2015-2022)")

gc_b <- bind_rows(gc_b_full, gc_b_pre, gc_b_post)

# Print price channel comparison
cat("TED SPREAD -> CREDIT SPREADS (comparison):\n")
print(gc_b %>% filter(X == "TED_spread") %>%
        select(Period, Y, Lag, N, F_stat, p_value, Sig) %>%
        arrange(Y, Lag, Period), row.names = FALSE)
cat("\n")

# Print quantity comparison
cat("QUANTITIES -> CREDIT SPREADS (comparison):\n")
print(gc_b %>% filter(X %in% c("d_Reserves", "d_TGA"),
                       Y %in% c("d_Baa", "d_Aaa", "d_Baa_Aaa")) %>%
        select(Period, X, Y, Lag, N, F_stat, p_value, Sig) %>%
        arrange(X, Y, Lag, Period), row.names = FALSE)
cat("\n")


# ── B3. Structural Break Tests ──

cat("================================================================\n")
cat("  B3: STRUCTURAL BREAK TESTS (TED SAMPLE)\n")
cat("================================================================\n\n")

# Chow at Jan 2015
chow_b <- data.frame()
for (pair in all_pairs_b) {
  for (lag in c(1, 5, 10)) {
    res <- run_chow(df_ted, pair[1], pair[2], lag, "2015-01-01")
    if (!is.null(res)) chow_b <- rbind(chow_b, res)
  }
}

cat("Chow test at January 2015:\n\n")
print(chow_b, row.names = FALSE)
cat(sprintf("\nSignificant at 5%%: %d / %d\n\n", sum(chow_b$p_value < 0.05), nrow(chow_b)))

# SupF (endogenous break)
supf_b <- data.frame()
for (pair in fwd_pairs_b) {
  for (lag in c(1, 5, 10)) {
    res <- run_supf(df_ted, pair[1], pair[2], lag)
    if (!is.null(res)) supf_b <- rbind(supf_b, res)
  }
}

cat("SupF test (endogenous break):\n\n")
print(supf_b, row.names = FALSE)
cat(sprintf("\nSupF significant at 5%%: %d / %d\n\n", sum(supf_b$p_value < 0.05), nrow(supf_b)))

# ── Endogenous break-date identification via breakpoints() ──
# The SupF test confirms a break exists but does NOT identify the date.
# We use strucchange::breakpoints() to find the data-optimal break date
# for each key specification and verify it aligns with the a priori 2015 choice.

cat("Endogenous break-date identification (breakpoints):\n\n")

bp_results <- data.frame()
for (pair in fwd_pairs_b) {
  for (lag in c(1, 5, 10)) {
    xy <- df_ted %>%
      select(Date, all_of(c(pair[1], pair[2]))) %>%
      tidyr::drop_na() %>%
      arrange(Date)

    n <- nrow(xy)
    if (n < lag + 50) next

    Y <- xy[[pair[2]]][(lag+1):n]
    X_mat <- data.frame(row.names = 1:length(Y))
    for (i in 1:lag) {
      X_mat[[paste0("Y_lag", i)]] <- xy[[pair[2]]][(lag+1-i):(n-i)]
      X_mat[[paste0("X_lag", i)]] <- xy[[pair[1]]][(lag+1-i):(n-i)]
    }
    X_mat$Y <- Y
    dates_used <- xy$Date[(lag+1):n]

    tryCatch({
      bp <- breakpoints(Y ~ ., data = X_mat, breaks = 1)
      if (!is.na(bp$breakpoints[1])) {
        bp_date <- dates_used[bp$breakpoints[1]]
        bp_results <- rbind(bp_results, data.frame(
          Direction  = paste0(pair[1], " -> ", pair[2]),
          Lag        = lag,
          Break_date = as.character(bp_date),
          Break_year = format(bp_date, "%Y"),
          stringsAsFactors = FALSE
        ))
      }
    }, error = function(e) {
      cat(sprintf("  breakpoints error: %s -> %s lag %d: %s\n",
                  pair[1], pair[2], lag, e$message))
    })
  }
}

print(bp_results, row.names = FALSE)
cat(sprintf("\nBreak dates falling in 2014-2015: %d / %d\n",
            sum(bp_results$Break_year %in% c("2014", "2015")), nrow(bp_results)))
cat(sprintf("Break dates falling in 2013-2016: %d / %d\n\n",
            sum(bp_results$Break_year %in% c("2013", "2014", "2015", "2016")), nrow(bp_results)))

# Save endogenous breakpoint results
write.csv(bp_results, "results_endogenous_breakpoints.csv", row.names = FALSE)

# Generate LaTeX table
bp_clean <- bp_results
bp_clean$Direction <- gsub("_", "\\\\_", bp_clean$Direction)
bp_clean$Direction <- gsub("TED\\\\_spread", "TED", bp_clean$Direction)
bp_clean$Direction <- gsub("d\\\\_Baa\\\\_Aaa", "$\\\\Delta$(Baa--Aaa)", bp_clean$Direction)
bp_clean$Direction <- gsub("d\\\\_Baa", "$\\\\Delta$Baa", bp_clean$Direction)
bp_clean$Direction <- gsub("d\\\\_Aaa", "$\\\\Delta$Aaa", bp_clean$Direction)
bp_clean$Direction <- gsub("d\\\\_Reserves", "$\\\\Delta$Reserves", bp_clean$Direction)
bp_clean$Direction <- gsub("d\\\\_TGA", "$\\\\Delta$TGA", bp_clean$Direction)

bp_tex <- xtable(bp_clean[, c("Direction", "Lag", "Break_date")],
                 caption = "Endogenous Break Dates (Bai-Perron, Single Break) — TED Sample",
                 label = "tab:endogenous_breakpoints")
print(bp_tex, file = "table_endogenous_breakpoints.tex",
      include.rownames = FALSE, booktabs = TRUE,
      sanitize.text.function = identity,
      scalebox = 0.78)
cat("  Saved: table_endogenous_breakpoints.tex\n")
cat("  Saved: results_endogenous_breakpoints.csv\n\n")

# Chow at GFC peak (Sep 2008)
cat("Chow test at September 2008 (GFC peak):\n\n")
chow_b_gfc <- data.frame()
for (pair in fwd_pairs_b[1:3]) {  # TED spread pairs only
  for (lag in c(1, 5, 10)) {
    res <- run_chow(df_ted, pair[1], pair[2], lag, "2008-09-15")
    if (!is.null(res)) chow_b_gfc <- rbind(chow_b_gfc, res)
  }
}
print(chow_b_gfc, row.names = FALSE)
cat("\n")


# ── B4. Figures (TED Spread) ──

cat("================================================================\n")
cat("  B4: FIGURES (TED SPREAD)\n")
cat("================================================================\n\n")

# Fig B1: TED spread with both break dates
p_b1 <- ggplot(df_ted, aes(x = Date, y = TED_spread)) +
  geom_line(color = "steelblue", linewidth = 0.35) +
  geom_vline(xintercept = as.Date("2008-09-15"), linetype = "dotted",
             color = "darkred", linewidth = 0.7) +
  geom_vline(xintercept = BREAK_DATE_B, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  annotate("text", x = as.Date("2008-09-15") - 250, y = 4.2,
           label = "Lehman\n(Sep 2008)", color = "darkred", size = 3, hjust = 1) +
  annotate("text", x = BREAK_DATE_B + 120, y = 4.2,
           label = "Structural break\n(Jan 2015)", color = "red", size = 3, hjust = 0) +
  labs(title = "TED Spread (2003-2022)",
       subtitle = "LIBOR minus T-bill rate; measures unsecured interbank funding stress",
       x = NULL, y = "TED spread (%)") +
  theme_thesis
ggsave("fig_b1_ted_spread.pdf", p_b1, width = 10, height = 5)
cat("  Saved: fig_b1_ted_spread.pdf\n")

# Fig B2: Credit spreads full sample
df_credit_b <- df_ted %>%
  filter(!is.na(Baa_spread) & !is.na(Aaa_spread)) %>%
  select(Date, Baa_spread, Aaa_spread) %>%
  tidyr::pivot_longer(-Date, names_to = "Variable", values_to = "Spread") %>%
  mutate(Variable = ifelse(Variable == "Baa_spread", "Baa", "Aaa"))

p_b2 <- ggplot(df_credit_b, aes(x = Date, y = Spread, color = Variable)) +
  geom_line(linewidth = 0.35) +
  geom_vline(xintercept = as.Date("2008-09-15"), linetype = "dotted",
             color = "darkred", linewidth = 0.7) +
  geom_vline(xintercept = BREAK_DATE_B, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  scale_color_manual(values = c("Baa" = "darkorange", "Aaa" = "steelblue")) +
  labs(title = "Moody's Credit Spreads (2003-2022)",
       subtitle = "Option-adjusted spreads; GFC and 2015 structural break marked",
       x = NULL, y = "Spread (%)", color = NULL) +
  theme_thesis
ggsave("fig_b2_credit_spreads_ted.pdf", p_b2, width = 10, height = 5)
cat("  Saved: fig_b2_credit_spreads_ted.pdf\n")

# Fig B3: TED density comparison
df_dens_b <- bind_rows(
  df_pre %>% mutate(Period = "Pre-2015") %>% select(Period, TED_spread),
  df_post_b %>% mutate(Period = "Post-2015") %>% select(Period, TED_spread)
)

p_b3 <- ggplot(df_dens_b, aes(x = TED_spread, fill = Period)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Pre-2015" = "steelblue", "Post-2015" = "coral")) +
  labs(title = "Distribution of TED Spread: Pre vs. Post-2015",
       subtitle = "Pre-2015 has heavy right tail from GFC; post-2015 is compressed",
       x = "TED spread (%)", y = "Density", fill = NULL) +
  theme_thesis
ggsave("fig_b3_density_ted.pdf", p_b3, width = 8, height = 5)
cat("  Saved: fig_b3_density_ted.pdf\n")

# Fig B4: Rolling Granger p-values (TED sample)
cat("  Computing rolling Granger p-values (TED sample)...\n")

roll_b_baa <- roll_granger(df_ted, "TED_spread", "d_Baa", lag = 5, window = 500, step = 20)
roll_b_baa$Direction <- "TED -> d_Baa (lag 5)"
roll_b_aaa <- roll_granger(df_ted, "TED_spread", "d_Aaa", lag = 5, window = 500, step = 20)
roll_b_aaa$Direction <- "TED -> d_Aaa (lag 5)"
roll_b_diff <- roll_granger(df_ted, "TED_spread", "d_Baa_Aaa", lag = 5, window = 500, step = 20)
roll_b_diff$Direction <- "TED -> d_Baa_Aaa (lag 5)"

roll_b <- bind_rows(roll_b_baa, roll_b_aaa, roll_b_diff)

if (nrow(roll_b) > 0) {
  p_b4 <- ggplot(roll_b, aes(x = Date, y = p_value, color = Direction)) +
    geom_line(linewidth = 0.5) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
    geom_vline(xintercept = as.Date("2008-09-15"), linetype = "dotted",
               color = "darkred", linewidth = 0.6) +
    geom_vline(xintercept = BREAK_DATE_B, linetype = "dashed",
               color = "red", linewidth = 0.6) +
    scale_color_manual(values = c("TED -> d_Baa (lag 5)" = "darkorange",
                                   "TED -> d_Aaa (lag 5)" = "steelblue",
                                   "TED -> d_Baa_Aaa (lag 5)" = "darkgreen")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Rolling Granger p-values: TED Spread -> Credit Spreads (2003-2022)",
         subtitle = "500-day window, lag 5; red dashed = 2015 break, dotted = Lehman",
         x = NULL, y = "p-value", color = NULL) +
    theme_thesis
  ggsave("fig_b4_rolling_granger_ted.pdf", p_b4, width = 10, height = 5.5)
  cat("  Saved: fig_b4_rolling_granger_ted.pdf\n")
}

# Fig B5: Reserves with regime annotations
df_res_b <- df_ted %>% filter(!is.na(Reserves))
p_b5 <- ggplot(df_res_b, aes(x = Date, y = Reserves / 1e6)) +
  geom_line(color = "darkgreen", linewidth = 0.4) +
  geom_vline(xintercept = BREAK_DATE_B, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  annotate("rect", xmin = as.Date("2008-12-01"), xmax = as.Date("2014-10-01"),
           ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "blue") +
  annotate("text", x = as.Date("2011-06-01"), y = max(df_res_b$Reserves/1e6)*0.95,
           label = "QE era", color = "blue", size = 3) +
  labs(title = "Bank Reserves at Federal Reserve (2003-2022)",
       subtitle = "Transition from scarce reserves (corridor) to abundant reserves (floor system)",
       x = NULL, y = "Trillions USD") +
  theme_thesis
ggsave("fig_b5_reserves_ted.pdf", p_b5, width = 10, height = 5)
cat("  Saved: fig_b5_reserves_ted.pdf\n")

# Fig B6: F-stat comparison bar chart (lag 5)
bar_b <- gc_b %>%
  filter(X == "TED_spread", Lag == 5) %>%
  mutate(Period = factor(Period, levels = c("Full (2003-2022)",
                                             "Pre-2015 (2003-2014)",
                                             "Post-2015 (2015-2022)")))

p_b6 <- ggplot(bar_b, aes(x = Y, y = F_stat, fill = Period)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = ifelse(p_value < 0.05,
                                sprintf("p=%.3f%s", p_value, Sig),
                                sprintf("p=%.2f", p_value))),
            position = position_dodge(0.7), vjust = -0.5, size = 2.5) +
  scale_fill_manual(values = c("Full (2003-2022)" = "grey60",
                                "Pre-2015 (2003-2014)" = "steelblue",
                                "Post-2015 (2015-2022)" = "coral")) +
  labs(title = "Granger F-statistics: TED Spread -> Credit Spreads (Lag 5)",
       subtitle = "Comparison across full, pre-2015, and post-2015 subsamples",
       x = NULL, y = "F-statistic", fill = NULL) +
  theme_thesis
ggsave("fig_b6_fstat_comparison.pdf", p_b6, width = 9, height = 5.5)
cat("  Saved: fig_b6_fstat_comparison.pdf\n")


###############################################################################
###############################################################################
#                     PART C — LATEX TABLES
###############################################################################
###############################################################################

cat("\n\n================================================================\n")
cat("  PART C: LATEX TABLES\n")
cat("================================================================\n\n")

# ── Table 1: Summary stats post-2018 ──
t1 <- stats_a %>% select(Variable, N, Mean, SD, Skewness, Kurtosis)
colnames(t1) <- c("Variable", "$N$", "Mean", "SD", "Skew.", "Kurt.")
t1_xt <- xtable(t1,
  caption = "Summary Statistics: Post-2018 Sample (Stationary Transformations)",
  label = "tab:sumstats_post2018",
  digits = c(0, 0, 0, 4, 4, 2, 2))
print(t1_xt, file = "table_sumstats_post2018.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity)
cat("  Saved: table_sumstats_post2018.tex\n")

# ── Table 2: ADF post-2018 ──
t2 <- adf_a %>% select(Variable, ADF_stat, p_value, N, Stationary)
colnames(t2) <- c("Variable", "ADF stat", "$p$-value", "$N$", "Stationary")
t2_xt <- xtable(t2,
  caption = "ADF Unit Root Tests (Post-2018 Sample)",
  label = "tab:adf_post2018",
  digits = c(0, 0, 3, 4, 0, 0))
print(t2_xt, file = "table_adf_post2018.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity)
cat("  Saved: table_adf_post2018.tex\n")

# ── Table 3: Granger forward post-2018 ──
t3 <- gc_fwd_a %>%
  mutate(Direction = paste0(X, " $\\rightarrow$ ", Y)) %>%
  select(Direction, Lag, N, F_stat, p_value, Sig)
colnames(t3) <- c("Direction", "Lag", "$N$", "$F$-stat", "$p$-value", "")
t3_xt <- xtable(t3,
  caption = "Granger Causality: Funding $\\rightarrow$ Credit Spreads (Post-2018)",
  label = "tab:granger_fwd_post2018",
  digits = c(0, 0, 0, 0, 3, 4, 0))
print(t3_xt, file = "table_granger_fwd_post2018.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78)
cat("  Saved: table_granger_fwd_post2018.tex\n")

# ── Table 4: Granger reverse post-2018 ──
t4 <- gc_rev_a %>%
  mutate(Direction = paste0(X, " $\\rightarrow$ ", Y)) %>%
  select(Direction, Lag, N, F_stat, p_value, Sig)
colnames(t4) <- c("Direction", "Lag", "$N$", "$F$-stat", "$p$-value", "")
t4_xt <- xtable(t4,
  caption = "Reverse Causality: Credit Spreads $\\rightarrow$ Funding (Post-2018)",
  label = "tab:granger_rev_post2018",
  digits = c(0, 0, 0, 0, 3, 4, 0))
print(t4_xt, file = "table_granger_rev_post2018.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78)
cat("  Saved: table_granger_rev_post2018.tex\n")

# ── Table 5: Chow + SupF post-2018 ──
t5_chow <- chow_a %>%
  filter(grepl("SOFR_EFFR|d_Reserves|d_TGA", Direction) &
         grepl("d_Baa|d_Aaa|d_Baa_Aaa", Direction) &
         !grepl("d_Baa.*SOFR|d_Aaa.*SOFR|d_Baa.*d_Res|d_Aaa.*d_Res|d_Baa.*d_TGA|d_Aaa.*d_TGA", Direction))
t5 <- t5_chow %>% select(Direction, Lag, Break_date, Chow_F, p_value, Sig)
colnames(t5) <- c("Direction", "Lag", "Break date", "Chow $F$", "$p$-value", "")
t5_xt <- xtable(t5,
  caption = "Chow Test for Structural Break in Granger Regressions (Post-2018, Midpoint)",
  label = "tab:chow_post2018",
  digits = c(0, 0, 0, 0, 3, 4, 0))
print(t5_xt, file = "table_chow_post2018.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78)
cat("  Saved: table_chow_post2018.tex\n")

# ── Table 6: SupF post-2018 ──
t6 <- supf_a %>% select(Direction, Lag, SupF, p_value, Sig)
colnames(t6) <- c("Direction", "Lag", "SupF", "$p$-value", "")
t6_xt <- xtable(t6,
  caption = "Andrews SupF Test for Unknown Structural Break (Post-2018)",
  label = "tab:supf_post2018",
  digits = c(0, 0, 0, 3, 4, 0))
print(t6_xt, file = "table_supf_post2018.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78)
cat("  Saved: table_supf_post2018.tex\n")

# ── Table 7: TED Granger comparison (price channel) ──
t7 <- gc_b %>%
  filter(X == "TED_spread") %>%
  mutate(Direction = paste0("TED $\\rightarrow$ ", Y)) %>%
  select(Period, Direction, Lag, N, F_stat, p_value, Sig) %>%
  arrange(Direction, Lag, Period)
colnames(t7) <- c("Sample", "Direction", "Lag", "$N$", "$F$-stat", "$p$-value", "")
t7_xt <- xtable(t7,
  caption = "Granger Causality: TED Spread $\\rightarrow$ Credit Spreads --- Structural Break Comparison (2003--2022)",
  label = "tab:granger_ted_price",
  digits = c(0, 0, 0, 0, 0, 3, 4, 0))
print(t7_xt, file = "table_granger_ted_price.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.72)
cat("  Saved: table_granger_ted_price.tex\n")

# ── Table 8: TED Granger comparison (quantity channel) ──
t8 <- gc_b %>%
  filter(X %in% c("d_Reserves", "d_TGA"),
         Y %in% c("d_Baa", "d_Aaa", "d_Baa_Aaa")) %>%
  mutate(Direction = paste0(X, " $\\rightarrow$ ", Y)) %>%
  select(Period, Direction, Lag, N, F_stat, p_value, Sig) %>%
  arrange(Direction, Lag, Period)
colnames(t8) <- c("Sample", "Direction", "Lag", "$N$", "$F$-stat", "$p$-value", "")
t8_xt <- xtable(t8,
  caption = "Granger Causality: Quantities $\\rightarrow$ Credit Spreads --- Structural Break Comparison (2003--2022)",
  label = "tab:granger_ted_qty",
  digits = c(0, 0, 0, 0, 0, 3, 4, 0))
print(t8_xt, file = "table_granger_ted_qty.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.68)
cat("  Saved: table_granger_ted_qty.tex\n")

# ── Table 9: TED reverse causality ──
t9 <- gc_b %>%
  filter(X %in% c("d_Baa", "d_Aaa"),
         Y %in% c("TED_spread", "d_Reserves", "d_TGA")) %>%
  mutate(Direction = paste0(X, " $\\rightarrow$ ", Y)) %>%
  select(Period, Direction, Lag, N, F_stat, p_value, Sig) %>%
  arrange(Direction, Lag, Period)
colnames(t9) <- c("Sample", "Direction", "Lag", "$N$", "$F$-stat", "$p$-value", "")
t9_xt <- xtable(t9,
  caption = "Reverse Causality: Credit $\\rightarrow$ Funding --- Structural Break Comparison (2003--2022)",
  label = "tab:granger_ted_reverse",
  digits = c(0, 0, 0, 0, 0, 3, 4, 0))
print(t9_xt, file = "table_granger_ted_reverse.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.68)
cat("  Saved: table_granger_ted_reverse.tex\n")

# ── Table 10: Chow TED at 2015 ──
t10 <- chow_b %>%
  filter(grepl("TED|d_Reserves|d_TGA", Direction) &
         grepl("d_Baa|d_Aaa|d_Baa_Aaa", Direction) &
         !grepl("d_Baa.*TED|d_Aaa.*TED|d_Baa.*d_Res|d_Aaa.*d_Res|d_Baa.*d_TGA|d_Aaa.*d_TGA", Direction)) %>%
  select(Direction, Lag, Break_date, Chow_F, p_value, Sig)
colnames(t10) <- c("Direction", "Lag", "Break date", "Chow $F$", "$p$-value", "")
t10_xt <- xtable(t10,
  caption = "Chow Test for Structural Break at January 2015 (TED Spread Sample)",
  label = "tab:chow_ted_2015",
  digits = c(0, 0, 0, 0, 3, 4, 0))
print(t10_xt, file = "table_chow_ted_2015.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78)
cat("  Saved: table_chow_ted_2015.tex\n")

# ── Table 11: SupF TED ──
t11 <- supf_b %>% select(Direction, Lag, SupF, p_value, Sig)
colnames(t11) <- c("Direction", "Lag", "SupF", "$p$-value", "")
t11_xt <- xtable(t11,
  caption = "Andrews SupF Test for Unknown Structural Break (TED Spread Sample, 2003--2022)",
  label = "tab:supf_ted",
  digits = c(0, 0, 0, 3, 4, 0))
print(t11_xt, file = "table_supf_ted.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78)
cat("  Saved: table_supf_ted.tex\n")

# ── Table 12: Summary stats TED by regime ──
t12 <- stats_b %>%
  filter(Variable %in% c("TED_spread", "d_Baa", "d_Aaa", "d_Reserves")) %>%
  select(Period, Variable, N, Mean, SD, Skewness, Kurtosis)
colnames(t12) <- c("Period", "Variable", "$N$", "Mean", "SD", "Skew.", "Kurt.")
t12_xt <- xtable(t12,
  caption = "Summary Statistics by Subsample: Pre vs.\\ Post-2015 Structural Break",
  label = "tab:sumstats_ted_break",
  digits = c(0, 0, 0, 0, 4, 4, 2, 2))
print(t12_xt, file = "table_sumstats_ted_break.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78)
cat("  Saved: table_sumstats_ted_break.tex\n")

# ── Table 13: ADF TED by regime ──
t13 <- adf_b %>% select(Period, Variable, ADF_stat, p_value, N, Stationary)
colnames(t13) <- c("Period", "Variable", "ADF stat", "$p$-value", "$N$", "Stationary")
t13_xt <- xtable(t13,
  caption = "ADF Unit Root Tests by Subsample (TED Spread Sample)",
  label = "tab:adf_ted_break",
  digits = c(0, 0, 0, 3, 4, 0, 0))
print(t13_xt, file = "table_adf_ted_break.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78)
cat("  Saved: table_adf_ted_break.tex\n")


###############################################################################
#  SAVE ALL RESULTS
###############################################################################

write.csv(gc_a, "results_granger_post2018.csv", row.names = FALSE)
write.csv(gc_b, "results_granger_ted.csv", row.names = FALSE)
write.csv(chow_a, "results_chow_post2018.csv", row.names = FALSE)
write.csv(chow_b, "results_chow_ted.csv", row.names = FALSE)
write.csv(supf_a, "results_supf_post2018.csv", row.names = FALSE)
write.csv(supf_b, "results_supf_ted.csv", row.names = FALSE)
write.csv(stats_a, "results_sumstats_post2018.csv", row.names = FALSE)
write.csv(stats_b, "results_sumstats_ted.csv", row.names = FALSE)
write.csv(adf_a, "results_adf_post2018.csv", row.names = FALSE)
write.csv(adf_b, "results_adf_ted.csv", row.names = FALSE)

cat("\n\n================================================================\n")
cat("  ALL RESULTS SAVED\n")
cat("================================================================\n\n")

# Final summary
cat("================================================================\n")
cat("  FINAL SUMMARY\n")
cat("================================================================\n\n")

cat("  PART A (Post-2018, SOFR-EFFR):\n")
n_fwd_a <- sum(gc_fwd_a$p_value < 0.05, na.rm = TRUE)
n_rev_a <- sum(gc_rev_a$p_value < 0.05, na.rm = TRUE)
cat(sprintf("    Forward significant: %d/%d\n", n_fwd_a, nrow(gc_fwd_a)))
cat(sprintf("    Reverse significant: %d/%d\n", n_rev_a, nrow(gc_rev_a)))
cat(sprintf("    Chow (midpoint) significant: %d/%d\n",
            sum(chow_a$p_value < 0.05), nrow(chow_a)))
cat(sprintf("    SupF significant: %d/%d\n\n", sum(supf_a$p_value < 0.05), nrow(supf_a)))

cat("  PART B (TED spread, 2003-2022):\n")
for (per in c("Full (2003-2022)", "Pre-2015 (2003-2014)", "Post-2015 (2015-2022)")) {
  sub <- gc_b %>% filter(Period == per, !is.na(p_value))
  fwd <- sub %>% filter(X %in% c("TED_spread", "d_Reserves", "d_TGA"),
                          Y %in% c("d_Baa", "d_Aaa", "d_Baa_Aaa"))
  cat(sprintf("    %s: %d/%d forward significant\n", per,
              sum(fwd$p_value < 0.05), nrow(fwd)))
}
cat(sprintf("    Chow (2015) significant: %d/%d\n",
            sum(chow_b$p_value < 0.05), nrow(chow_b)))
cat(sprintf("    SupF significant: %d/%d\n\n", sum(supf_b$p_value < 0.05), nrow(supf_b)))

cat("================================================================\n")
cat("  STEP 2 (ENHANCED) COMPLETE\n")
cat("================================================================\n")
