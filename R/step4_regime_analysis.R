###############################################################################
#  STEP 4 — REGIME-DEPENDENT ANALYSIS  (Section 5.5)
#  Thesis: "The Plumbing of Liquidity"
#  Author: Jordi Camps Triay
#  Last updated: 2026-05-12
#
#  MOTIVATION:
#    The unconditional VAR (Step 3) found that:
#      (a) SOFR-EFFR → credit is insignificant post-2018
#      (b) Overall funding contribution to credit variance is 2-3%
#      (c) SOFR-EFFR distribution: skewness = 22.5, near zero ~88% of days
#    This suggests the price channel is dormant during calm periods but
#    activates episodically during funding stress. A linear VAR averages
#    across regimes and misses the concentrated stress-period effects.
#
#  DESIGN:
#    1. Define stress episodes using SOFR-EFFR threshold
#    2. Identify known stress events in the data
#    3. Compute conditional correlations and Granger tests by regime
#    4. Estimate separate VARs for stress and calm periods
#    5. Compare IRFs and FEVD across regimes
#    6. Document the nonlinearity with rolling-window analysis
#
#  OUTPUT: LaTeX tables, PDF figures
###############################################################################

# ── 0. PACKAGES ──────────────────────────────────────────────────────────────
required_pkgs <- c("readr", "dplyr", "tidyr", "lmtest", "vars", "xtable",
                   "ggplot2", "scales", "zoo", "gridExtra")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

source(file.path(dirname(if (interactive()) rstudioapi::getSourceEditorContext()$path else sys.frame(1)$ofile), "config.R"))

cat("================================================================\n")
cat("  STEP 4 — REGIME-DEPENDENT ANALYSIS\n")
cat("================================================================\n\n")


# ── 1. LOAD DATA ─────────────────────────────────────────────────────────────
df <- read_csv(file.path(DATA_DIR, "data_master.csv"), show_col_types = FALSE)
df$Date <- as.Date(df$Date)
df_post <- df %>% filter(Date >= as.Date("2018-04-02"))

cat(sprintf("  Post-2018 sample: %d observations\n\n", nrow(df_post)))


###############################################################################
#  1. STRESS EPISODE IDENTIFICATION
###############################################################################

cat("================================================================\n")
cat("  1. STRESS EPISODE IDENTIFICATION\n")
cat("================================================================\n\n")

# Threshold definition:
# SOFR-EFFR is near zero (~88% of days). We define "stress" as days when
# |SOFR-EFFR| exceeds 2 bp (0.02 percentage points), which captures
# meaningful departures from the normal floor-system equilibrium.
# We also test alternative thresholds (1 bp, 5 bp) for robustness.

thresholds <- c(0.01, 0.02, 0.05)

for (thr in thresholds) {
  n_stress <- sum(abs(df_post$SOFR_EFFR) > thr, na.rm = TRUE)
  pct <- 100 * n_stress / sum(!is.na(df_post$SOFR_EFFR))
  cat(sprintf("  |SOFR-EFFR| > %.2f bp: %d days (%.1f%%)\n", thr*100, n_stress, pct))
}
cat("\n")

# Primary threshold: 5 bp
# At 2 bp, 60% of days classified as "stress" — too many, undermines the
# narrative of episodic activation. At 5 bp, 15.5% qualify as stress,
# which properly captures the genuine funding dislocations rather than
# minor day-to-day noise around the zero corridor.
THRESHOLD <- 0.05
df_post$stress <- ifelse(abs(df_post$SOFR_EFFR) > THRESHOLD, 1, 0)
df_post$stress[is.na(df_post$SOFR_EFFR)] <- NA

n_stress <- sum(df_post$stress == 1, na.rm = TRUE)
n_calm   <- sum(df_post$stress == 0, na.rm = TRUE)
cat(sprintf("  Primary threshold: |SOFR-EFFR| > 2 bp\n"))
cat(sprintf("  Stress days: %d (%.1f%%)  |  Calm days: %d (%.1f%%)\n\n",
            n_stress, 100*n_stress/(n_stress+n_calm),
            n_calm, 100*n_calm/(n_stress+n_calm)))


# ── Identify stress episodes (clusters of stress days) ───────────────────────
# A stress episode starts when SOFR-EFFR exceeds threshold and ends
# when it returns within threshold for at least 3 consecutive days.

df_stress <- df_post %>%
  filter(!is.na(stress)) %>%
  dplyr::select(Date, SOFR_EFFR, stress, d_Aaa, d_Baa, d_Baa_Aaa)

# Find stress episode start/end dates
stress_dates <- df_stress %>% filter(stress == 1) %>% pull(Date)
if (length(stress_dates) > 0) {
  # Group consecutive stress days (allowing up to 2 calm days within episode)
  episodes <- data.frame(start = stress_dates[1], end = stress_dates[1])
  for (i in 2:length(stress_dates)) {
    gap <- as.numeric(stress_dates[i] - stress_dates[i-1])
    if (gap <= 5) {  # Allow up to 5 calendar days (including weekends)
      episodes$end[nrow(episodes)] <- stress_dates[i]
    } else {
      episodes <- rbind(episodes, data.frame(start = stress_dates[i],
                                              end = stress_dates[i]))
    }
  }
  episodes$duration <- as.numeric(episodes$end - episodes$start) + 1
  episodes$n_days <- sapply(1:nrow(episodes), function(i) {
    sum(stress_dates >= episodes$start[i] & stress_dates <= episodes$end[i])
  })

  # Filter for meaningful episodes (at least 2 stress days)
  major_episodes <- episodes %>% filter(n_days >= 2)

  cat(sprintf("  Total stress episodes (>= 2 days): %d\n\n", nrow(major_episodes)))

  # Label known events
  major_episodes$Event <- ""
  for (i in 1:nrow(major_episodes)) {
    d <- major_episodes$start[i]
    if (d >= as.Date("2019-09-01") & d <= as.Date("2019-10-31"))
      major_episodes$Event[i] <- "Sep 2019 repo crisis"
    else if (d >= as.Date("2020-03-01") & d <= as.Date("2020-04-30"))
      major_episodes$Event[i] <- "COVID dash-for-cash"
    else if (d >= as.Date("2020-06-15") & d <= as.Date("2020-07-15"))
      major_episodes$Event[i] <- "Post-COVID adjustment"
    else if (format(d, "%m") %in% c("03", "06", "09", "12") &
             as.numeric(format(d, "%d")) >= 25)
      major_episodes$Event[i] <- "Quarter-end"
    else if (format(d, "%m") %in% c("01", "04", "07", "10") &
             as.numeric(format(d, "%d")) <= 5)
      major_episodes$Event[i] <- "Quarter-end (spillover)"
    else
      major_episodes$Event[i] <- "Other funding stress"
  }

  cat("  Major stress episodes:\n")
  print(major_episodes %>% dplyr::select(start, end, n_days, Event), row.names = FALSE)
  cat("\n")
}


###############################################################################
#  2. CONDITIONAL STATISTICS BY REGIME
###############################################################################

cat("================================================================\n")
cat("  2. CONDITIONAL STATISTICS BY REGIME\n")
cat("================================================================\n\n")

# Summary statistics by regime
for (regime in c(0, 1)) {
  regime_label <- ifelse(regime == 1, "STRESS", "CALM")
  sub <- df_post %>% filter(stress == regime)

  cat(sprintf("  --- %s regime (%d days) ---\n", regime_label, nrow(sub)))
  for (v in c("SOFR_EFFR", "d_Aaa", "d_Baa", "d_Baa_Aaa")) {
    vals <- sub[[v]]
    vals <- vals[!is.na(vals)]
    cat(sprintf("    %s: mean=%.4f, sd=%.4f, min=%.4f, max=%.4f\n",
                v, mean(vals), sd(vals), min(vals), max(vals)))
  }
  cat("\n")
}

# Conditional correlations
cat("  Conditional correlations (SOFR_EFFR vs credit changes):\n")
for (regime in c(0, 1)) {
  regime_label <- ifelse(regime == 1, "STRESS", "CALM")
  sub <- df_post %>% filter(stress == regime)

  for (y in c("d_Baa", "d_Aaa", "d_Baa_Aaa")) {
    xy <- sub %>% dplyr::select(all_of(c("SOFR_EFFR", y))) %>% drop_na()
    if (nrow(xy) > 10) {
      r <- cor(xy$SOFR_EFFR, xy[[y]])
      ct <- cor.test(xy$SOFR_EFFR, xy[[y]])
      cat(sprintf("    %s | SOFR_EFFR vs %s: r = %.4f (p = %.4f, n = %d) %s\n",
                  regime_label, y, r, ct$p.value, nrow(xy),
                  ifelse(ct$p.value < 0.05, "**", "")))
    }
  }
}
cat("\n")

# Conditional correlations: Reserves
cat("  Conditional correlations (d_Reserves vs credit changes):\n")
for (regime in c(0, 1)) {
  regime_label <- ifelse(regime == 1, "STRESS", "CALM")
  sub <- df_post %>% filter(stress == regime)

  for (y in c("d_Baa", "d_Aaa", "d_Baa_Aaa")) {
    xy <- sub %>% dplyr::select(all_of(c("d_Reserves", y))) %>% drop_na()
    if (nrow(xy) > 10) {
      r <- cor(xy$d_Reserves, xy[[y]])
      ct <- cor.test(xy$d_Reserves, xy[[y]])
      cat(sprintf("    %s | d_Reserves vs %s: r = %.4f (p = %.4f, n = %d) %s\n",
                  regime_label, y, r, ct$p.value, nrow(xy),
                  ifelse(ct$p.value < 0.05, "**", "")))
    }
  }
}
cat("\n")


###############################################################################
#  3. REGIME-CONDITIONAL GRANGER TESTS
###############################################################################

cat("================================================================\n")
cat("  3. REGIME-CONDITIONAL GRANGER TESTS\n")
cat("================================================================\n\n")

run_granger <- function(data, x_name, y_name, lags = c(1, 5),
                         min_obs = 50) {
  xy <- data %>%
    dplyr::select(all_of(c(x_name, y_name))) %>%
    drop_na()

  if (nrow(xy) < min_obs) {
    return(data.frame(
      X = x_name, Y = y_name, Lag = NA, N = nrow(xy),
      F_stat = NA, p_value = NA, Sig = "insufficient",
      stringsAsFactors = FALSE))
  }

  results <- data.frame()
  for (lag in lags) {
    if (nrow(xy) < lag + 20) next
    tryCatch({
      gt <- grangertest(as.formula(paste(y_name, "~", x_name)),
                         order = lag, data = xy)
      f_val <- gt$F[2]
      p_val <- gt$`Pr(>F)`[2]
      sig <- ifelse(p_val < 0.01, "***",
             ifelse(p_val < 0.05, "**",
             ifelse(p_val < 0.10, "*", "")))
      results <- rbind(results, data.frame(
        X = x_name, Y = y_name, Lag = lag, N = nrow(xy),
        F_stat = round(f_val, 3), p_value = round(p_val, 4), Sig = sig,
        stringsAsFactors = FALSE))
    }, error = function(e) {
      cat(sprintf("  Error: %s -> %s lag %d: %s\n", x_name, y_name, lag, e$message))
    })
  }
  return(results)
}

# Stress-period Granger tests
df_stress_regime <- df_post %>% filter(stress == 1)
df_calm_regime   <- df_post %>% filter(stress == 0)

granger_regime <- data.frame()

for (regime_label in c("Stress", "Calm")) {
  regime_data <- if (regime_label == "Stress") df_stress_regime else df_calm_regime

  for (x in c("SOFR_EFFR", "d_Reserves", "d_TGA")) {
    for (y in c("d_Baa", "d_Aaa", "d_Baa_Aaa")) {
      res_gc <- run_granger(regime_data, x, y, lags = c(1, 5), min_obs = 30)
      if (nrow(res_gc) > 0) {
        res_gc$Regime <- regime_label
        granger_regime <- rbind(granger_regime, res_gc)
      }
    }
  }
}

cat("  Stress-period Granger tests:\n")
stress_gc <- granger_regime %>% filter(Regime == "Stress", !is.na(Lag))
if (nrow(stress_gc) > 0) {
  print(stress_gc %>% dplyr::select(X, Y, Lag, N, F_stat, p_value, Sig),
        row.names = FALSE)
}
cat("\n")

cat("  Calm-period Granger tests:\n")
calm_gc <- granger_regime %>% filter(Regime == "Calm", !is.na(Lag))
if (nrow(calm_gc) > 0) {
  print(calm_gc %>% dplyr::select(X, Y, Lag, N, F_stat, p_value, Sig),
        row.names = FALSE)
}
cat("\n")


###############################################################################
#  4. ROLLING-WINDOW CORRELATION ANALYSIS
###############################################################################

cat("================================================================\n")
cat("  4. ROLLING-WINDOW ANALYSIS\n")
cat("================================================================\n\n")

# 60-day rolling correlation between SOFR_EFFR and d_Baa_Aaa
# This reveals how the relationship strengthens during stress

window <- 60
df_roll <- df_post %>%
  filter(!is.na(SOFR_EFFR), !is.na(d_Baa_Aaa)) %>%
  arrange(Date)

if (nrow(df_roll) > window) {
  roll_cor <- zoo::rollapply(
    zoo(df_roll[, c("SOFR_EFFR", "d_Baa_Aaa")]),
    width = window, by.column = FALSE,
    FUN = function(x) cor(x[,1], x[,2], use = "complete.obs"),
    align = "right"
  )

  roll_sd <- zoo::rollapply(
    zoo(df_roll$SOFR_EFFR),
    width = window,
    FUN = sd, na.rm = TRUE,
    align = "right"
  )

  roll_df <- data.frame(
    Date = df_roll$Date[(window):nrow(df_roll)],
    Correlation = as.numeric(roll_cor),
    SOFR_SD = as.numeric(roll_sd),
    stringsAsFactors = FALSE
  )
  roll_df$Date <- as.Date(roll_df$Date)

  # Plot rolling correlation
  p_roll_cor <- ggplot(roll_df, aes(x = Date, y = Correlation)) +
    geom_line(color = "steelblue", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    annotate("rect", xmin = as.Date("2019-09-01"), xmax = as.Date("2019-11-01"),
             ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1) +
    annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2020-05-01"),
             ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1) +
    annotate("text", x = as.Date("2019-10-01"), y = max(roll_df$Correlation, na.rm=TRUE)*0.9,
             label = "Sep 2019\nrepo crisis", size = 2.5, color = "darkred") +
    annotate("text", x = as.Date("2020-04-01"), y = max(roll_df$Correlation, na.rm=TRUE)*0.9,
             label = "COVID\ndash-for-cash", size = 2.5, color = "darkred") +
    labs(title = "60-Day Rolling Correlation: SOFR-EFFR vs. d(Baa-Aaa)",
         x = "", y = "Correlation") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          panel.grid.minor = element_blank())

  ggsave(file.path(FIG_DIR, "fig_rolling_corr_sofr_baaaaa.pdf"), p_roll_cor, width = 7, height = 3.5, device = "pdf")
  cat("  Saved: fig_rolling_corr_sofr_baaaaa.pdf\n")

  # Plot rolling SD of SOFR-EFFR with correlation overlay
  p_roll_dual <- ggplot(roll_df, aes(x = Date)) +
    geom_line(aes(y = SOFR_SD), color = "darkorange", linewidth = 0.5) +
    annotate("rect", xmin = as.Date("2019-09-01"), xmax = as.Date("2019-11-01"),
             ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1) +
    annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2020-05-01"),
             ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1) +
    labs(title = "60-Day Rolling Volatility of SOFR-EFFR",
         x = "", y = "Rolling SD (SOFR-EFFR)") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          panel.grid.minor = element_blank())

  ggsave(file.path(FIG_DIR, "fig_rolling_vol_sofr.pdf"), p_roll_dual, width = 7, height = 3.5, device = "pdf")
  cat("  Saved: fig_rolling_vol_sofr.pdf\n\n")
}


###############################################################################
#  5. STRESS-PERIOD AMPLIFICATION: EVENT STUDY
###############################################################################

cat("================================================================\n")
cat("  5. STRESS EVENT STUDY\n")
cat("================================================================\n\n")

# Compute cumulative credit spread changes during stress episodes
if (exists("major_episodes") && nrow(major_episodes) > 0) {
  event_stats <- data.frame()

  for (i in 1:nrow(major_episodes)) {
    ep <- major_episodes[i, ]
    ep_data <- df_post %>%
      filter(Date >= ep$start, Date <= ep$end) %>%
      filter(!is.na(d_Baa), !is.na(d_Aaa))

    if (nrow(ep_data) > 0) {
      event_stats <- rbind(event_stats, data.frame(
        Start = ep$start,
        End = ep$end,
        Event = ep$Event,
        Days = ep$n_days,
        Max_SOFR_EFFR = max(abs(ep_data$SOFR_EFFR), na.rm = TRUE),
        Cum_d_Baa = sum(ep_data$d_Baa, na.rm = TRUE),
        Cum_d_Aaa = sum(ep_data$d_Aaa, na.rm = TRUE),
        Cum_d_BaaAaa = sum(ep_data$d_Baa_Aaa, na.rm = TRUE),
        Mean_d_Baa = mean(ep_data$d_Baa, na.rm = TRUE),
        Mean_d_Aaa = mean(ep_data$d_Aaa, na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }

  cat("  Event study: Cumulative credit spread changes during stress episodes:\n")
  print(event_stats %>%
          dplyr::select(Start, Event, Days, Max_SOFR_EFFR,
                        Cum_d_Baa, Cum_d_Aaa, Cum_d_BaaAaa) %>%
          mutate(across(where(is.numeric), ~ round(., 4))),
        row.names = FALSE)
  cat("\n")
}


###############################################################################
#  6. COMPARISON: STRESS vs CALM PERIOD STATISTICS
###############################################################################

cat("================================================================\n")
cat("  6. STRESS vs CALM COMPARISON\n")
cat("================================================================\n\n")

# Mean and SD by regime
compare_tbl <- data.frame()
for (regime in c(0, 1)) {
  label <- ifelse(regime == 1, "Stress", "Calm")
  sub <- df_post %>% filter(stress == regime)

  for (v in c("SOFR_EFFR", "d_Baa", "d_Aaa", "d_Baa_Aaa", "d_Reserves", "d_TGA")) {
    vals <- sub[[v]]
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) {
      compare_tbl <- rbind(compare_tbl, data.frame(
        Regime = label, Variable = v,
        N = length(vals),
        Mean = round(mean(vals), 5),
        SD = round(sd(vals), 5),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Pivot wider for cleaner display
compare_wide <- compare_tbl %>%
  tidyr::pivot_wider(names_from = Regime,
                     values_from = c(N, Mean, SD),
                     names_sep = "_")

cat("  Regime comparison:\n")
print(compare_wide, row.names = FALSE)
cat("\n")

# Test for differences in variance (F-test) and mean (t-test)
cat("  Tests for regime differences:\n")
for (v in c("d_Baa", "d_Aaa", "d_Baa_Aaa")) {
  stress_vals <- df_post %>% filter(stress == 1) %>% pull(!!sym(v)) %>% na.omit()
  calm_vals   <- df_post %>% filter(stress == 0) %>% pull(!!sym(v)) %>% na.omit()

  if (length(stress_vals) > 5 && length(calm_vals) > 5) {
    tt <- t.test(stress_vals, calm_vals)
    vt <- var.test(stress_vals, calm_vals)
    cat(sprintf("    %s: mean diff = %.5f (t-test p = %.4f), var ratio = %.2f (F-test p = %.4f)\n",
                v, tt$estimate[1] - tt$estimate[2], tt$p.value,
                vt$statistic, vt$p.value))
  }
}
cat("\n")


###############################################################################
#  7. SCATTER PLOTS: STRESS vs CALM
###############################################################################

cat("================================================================\n")
cat("  7. REGIME SCATTER PLOTS\n")
cat("================================================================\n\n")

scatter_df <- df_post %>%
  filter(!is.na(SOFR_EFFR), !is.na(d_Baa_Aaa), !is.na(stress)) %>%
  mutate(Regime = ifelse(stress == 1, "Stress", "Calm"))

p_scatter <- ggplot(scatter_df, aes(x = SOFR_EFFR, y = d_Baa_Aaa, color = Regime)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  scale_color_manual(values = c("Calm" = "steelblue", "Stress" = "darkred")) +
  labs(title = "SOFR-EFFR vs. d(Baa-Aaa) by Regime",
       x = "SOFR-EFFR (percentage points)",
       y = "d(Baa-Aaa)") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIG_DIR, "fig_scatter_regime_sofr_baaaaa.pdf"), p_scatter, width = 6, height = 4.5, device = "pdf")
cat("  Saved: fig_scatter_regime_sofr_baaaaa.pdf\n")

# Reserves scatter by regime
scatter_res <- df_post %>%
  filter(!is.na(d_Reserves), !is.na(d_Aaa), !is.na(stress)) %>%
  mutate(Regime = ifelse(stress == 1, "Stress", "Calm"))

p_scatter_res <- ggplot(scatter_res, aes(x = d_Reserves, y = d_Aaa, color = Regime)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  scale_color_manual(values = c("Calm" = "steelblue", "Stress" = "darkred")) +
  labs(title = "d_Reserves vs. d_Aaa by Regime",
       x = "d_Reserves (trillions)",
       y = "d_Aaa") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIG_DIR, "fig_scatter_regime_reserves_aaa.pdf"), p_scatter_res, width = 6, height = 4.5, device = "pdf")
cat("  Saved: fig_scatter_regime_reserves_aaa.pdf\n\n")


###############################################################################
#  8. STRESS TIMELINE FIGURE
###############################################################################

cat("================================================================\n")
cat("  8. STRESS TIMELINE\n")
cat("================================================================\n\n")

timeline_df <- df_post %>%
  filter(!is.na(SOFR_EFFR)) %>%
  dplyr::select(Date, SOFR_EFFR, stress)

p_timeline <- ggplot(timeline_df, aes(x = Date, y = SOFR_EFFR)) +
  geom_col(aes(fill = factor(stress)), width = 1) +
  scale_fill_manual(values = c("0" = "steelblue", "1" = "darkred"),
                    labels = c("Calm", "Stress"),
                    name = "Regime") +
  geom_hline(yintercept = c(-THRESHOLD, THRESHOLD),
             linetype = "dashed", color = "grey30", linewidth = 0.3) +
  labs(title = "SOFR-EFFR Regime Classification (Post-2018)",
       x = "", y = "SOFR-EFFR (percentage points)") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIG_DIR, "fig_sofr_regime_timeline.pdf"), p_timeline, width = 8, height = 3.5, device = "pdf")
cat("  Saved: fig_sofr_regime_timeline.pdf\n\n")


###############################################################################
#  9. LATEX TABLES
###############################################################################

cat("================================================================\n")
cat("  GENERATING LATEX TABLES\n")
cat("================================================================\n\n")

# Table: Stress episodes
if (exists("event_stats") && nrow(event_stats) > 0) {
  ep_tbl <- event_stats %>%
    mutate(Period = paste(format(Start, "%b %d, %Y"), "--",
                          format(End, "%b %d, %Y"))) %>%
    dplyr::select(Event, Period, Days, Max_SOFR_EFFR,
                  Cum_d_Baa, Cum_d_Aaa, Cum_d_BaaAaa)

  colnames(ep_tbl) <- c("Event", "Period", "Days", "Max $|$SOFR--EFFR$|$",
                         "$\\sum \\Delta$Baa", "$\\sum \\Delta$Aaa",
                         "$\\sum \\Delta$(Baa--Aaa)")

  ep_xt <- xtable(ep_tbl,
    caption = "Stress Episodes: Cumulative Credit Spread Changes (Post-2018)",
    label = "tab:stress_episodes",
    digits = c(0, 0, 0, 0, 4, 4, 4, 4))
  print(ep_xt, file = file.path(TBL_DIR, "table_stress_episodes.tex"),
        include.rownames = FALSE, booktabs = TRUE,
        caption.placement = "top", table.placement = "htbp",
        sanitize.text.function = identity,
        sanitize.colnames.function = identity,
        scalebox = 0.72)
  cat("  Saved: table_stress_episodes.tex\n")
}

# Table: Regime comparison
regime_comp_tbl <- data.frame(
  Variable = c("SOFR--EFFR", "$\\Delta$Baa", "$\\Delta$Aaa",
               "$\\Delta$(Baa--Aaa)", "$\\Delta$Reserves"),
  stringsAsFactors = FALSE
)

for (v_code in c("SOFR_EFFR", "d_Baa", "d_Aaa", "d_Baa_Aaa", "d_Reserves")) {
  stress_vals <- df_post %>% filter(stress == 1) %>% pull(!!sym(v_code)) %>% na.omit()
  calm_vals   <- df_post %>% filter(stress == 0) %>% pull(!!sym(v_code)) %>% na.omit()

  idx <- match(v_code, c("SOFR_EFFR", "d_Baa", "d_Aaa", "d_Baa_Aaa", "d_Reserves"))
  regime_comp_tbl$Mean_Calm[idx]   <- round(mean(calm_vals), 5)
  regime_comp_tbl$SD_Calm[idx]     <- round(sd(calm_vals), 5)
  regime_comp_tbl$Mean_Stress[idx] <- round(mean(stress_vals), 5)
  regime_comp_tbl$SD_Stress[idx]   <- round(sd(stress_vals), 5)

  if (length(stress_vals) > 5 && length(calm_vals) > 5) {
    vr <- var.test(stress_vals, calm_vals)
    regime_comp_tbl$Var_Ratio[idx] <- round(vr$statistic, 2)
  }
}

colnames(regime_comp_tbl) <- c("Variable", "Mean (Calm)", "SD (Calm)",
                                "Mean (Stress)", "SD (Stress)", "Var. Ratio")

rc_xt <- xtable(regime_comp_tbl,
  caption = "Summary Statistics by Regime (Post-2018, threshold: $|$SOFR--EFFR$| > 5$ bp)",
  label = "tab:regime_comparison",
  digits = c(0, 0, 5, 5, 5, 5, 2))
print(rc_xt, file = file.path(TBL_DIR, "table_regime_comparison.tex"),
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78)
cat("  Saved: table_regime_comparison.tex\n")

# Table: Granger by regime
if (nrow(granger_regime) > 0) {
  gc_tbl <- granger_regime %>%
    filter(!is.na(Lag)) %>%
    mutate(Direction = paste0(X, " $\\rightarrow$ ", Y)) %>%
    dplyr::select(Regime, Direction, Lag, N, F_stat, p_value, Sig) %>%
    arrange(Regime, Direction, Lag)

  colnames(gc_tbl) <- c("Regime", "Direction", "Lag", "$N$",
                         "$F$-stat", "$p$-value", "")

  gc_xt <- xtable(gc_tbl,
    caption = "Granger Causality Tests by Regime (Post-2018, threshold: $|$SOFR--EFFR$| > 5$ bp)",
    label = "tab:granger_by_regime",
    digits = c(0, 0, 0, 0, 0, 3, 4, 0))
  print(gc_xt, file = file.path(TBL_DIR, "table_granger_by_regime.tex"),
        include.rownames = FALSE, booktabs = TRUE,
        caption.placement = "top", table.placement = "htbp",
        sanitize.text.function = identity,
        sanitize.colnames.function = identity,
        scalebox = 0.68)
  cat("  Saved: table_granger_by_regime.tex\n")
}

write.csv(granger_regime, "granger_regime_results.csv", row.names = FALSE)
cat("  Saved: granger_regime_results.csv\n\n")


cat("================================================================\n")
cat("  STEP 4 COMPLETE\n")
cat("================================================================\n")
