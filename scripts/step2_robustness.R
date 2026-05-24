###############################################################################
#  STEP 2 — ROBUSTNESS CHECKS & ECONOMETRIC DIAGNOSTICS
#  Supplements step2_granger_enhanced.R
#
#  New computations:
#    1. AIC/BIC optimal lag selection for each bivariate Granger pair
#    2. Breusch-Godfrey serial correlation tests
#    3. ARCH-LM heteroskedasticity tests
#    4. HAC-robust (Newey-West) Granger F-tests
#    5. Toda-Yamamoto procedure for key results
#    6. COVID robustness (excluding Feb-May 2020)
#    7. Rolling Granger for quantity channel (reserves -> Aaa)
#
#  OUTPUT: LaTeX tables, PDF figures, CSV results
###############################################################################

# ── 0. PACKAGES ────────────────────────────────────────────────────────────
required_pkgs <- c("readr", "dplyr", "tidyr", "lmtest", "tseries",
                   "strucchange", "moments", "xtable", "ggplot2",
                   "scales", "zoo", "gridExtra", "sandwich")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

WD <- "/Users/jordi/Downloads/University/TFG/Data/Thesis/step2_granger_enhanced"
setwd(WD)

df <- read_csv("../data_master.csv", show_col_types = FALSE)
df$Date <- as.Date(df$Date)

theme_thesis <- theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# ── Post-2018 sample ──
df_post <- df %>%
  filter(Date >= as.Date("2018-04-03")) %>%
  select(Date, SOFR_EFFR, d_Baa, d_Aaa, d_Baa_Aaa,
         d_Reserves, d_TGA, d_ON_RRP) %>%
  tidyr::drop_na(SOFR_EFFR)

cat("================================================================\n")
cat("  ROBUSTNESS CHECKS & ECONOMETRIC DIAGNOSTICS\n")
cat("================================================================\n\n")
cat(sprintf("Post-2018 sample: %d observations (%s to %s)\n\n",
            nrow(df_post), min(df_post$Date), max(df_post$Date)))


###############################################################################
#  1. AIC/BIC OPTIMAL LAG SELECTION
###############################################################################
cat("── 1. AIC/BIC Lag Selection ──────────────────────────────────────\n")

select_lag_ic <- function(data, x_name, y_name, max_lag = 15) {
  xy <- data %>% select(all_of(c(x_name, y_name))) %>% tidyr::drop_na()
  n <- nrow(xy)
  if (n < max_lag + 30) return(data.frame(X=x_name, Y=y_name,
                                           AIC_lag=NA, BIC_lag=NA, N=n))
  aic_vals <- bic_vals <- numeric(max_lag)
  for (p in 1:max_lag) {
    Y <- xy[[y_name]][(p+1):n]
    X_mat <- data.frame(row.names = 1:length(Y))
    for (i in 1:p) {
      X_mat[[paste0("Y_lag", i)]] <- xy[[y_name]][(p+1-i):(n-i)]
      X_mat[[paste0("X_lag", i)]] <- xy[[x_name]][(p+1-i):(n-i)]
    }
    X_mat$Y <- Y
    fit <- lm(Y ~ ., data = X_mat)
    aic_vals[p] <- AIC(fit)
    bic_vals[p] <- BIC(fit)
  }
  data.frame(X = x_name, Y = y_name,
             AIC_lag = which.min(aic_vals), BIC_lag = which.min(bic_vals),
             N = n, stringsAsFactors = FALSE)
}

# Forward pairs (Post-2018)
fwd_pairs <- expand.grid(
  X = c("SOFR_EFFR", "d_Reserves", "d_TGA", "d_ON_RRP"),
  Y = c("d_Baa", "d_Aaa", "d_Baa_Aaa"),
  stringsAsFactors = FALSE
)
# Reverse pairs
rev_pairs <- expand.grid(
  X = c("d_Baa", "d_Aaa"),
  Y = c("SOFR_EFFR", "d_Reserves", "d_TGA"),
  stringsAsFactors = FALSE
)

lag_sel <- do.call(rbind, lapply(1:nrow(fwd_pairs), function(i) {
  select_lag_ic(df_post, fwd_pairs$X[i], fwd_pairs$Y[i])
}))
lag_sel_rev <- do.call(rbind, lapply(1:nrow(rev_pairs), function(i) {
  select_lag_ic(df_post, rev_pairs$X[i], rev_pairs$Y[i])
}))
lag_sel_all <- rbind(lag_sel, lag_sel_rev)

cat("AIC/BIC lag selection results:\n")
print(lag_sel_all)
write.csv(lag_sel_all, "results_lag_selection.csv", row.names = FALSE)

# Run Granger at AIC-selected lags for key pairs
run_granger_at_lag <- function(data, x_name, y_name, lag) {
  xy <- data %>% select(all_of(c(x_name, y_name))) %>% tidyr::drop_na()
  if (nrow(xy) < lag + 30) return(NULL)
  tryCatch({
    gt <- grangertest(as.formula(paste(y_name, "~", x_name)),
                       order = lag, data = xy)
    f_val <- gt$F[2]; p_val <- gt$`Pr(>F)`[2]
    sig <- ifelse(p_val < 0.01, "***", ifelse(p_val < 0.05, "**",
           ifelse(p_val < 0.10, "*", "")))
    data.frame(X = x_name, Y = y_name, Lag = lag,
               N = nrow(xy), F_stat = round(f_val, 3),
               p_value = round(p_val, 4), Sig = sig,
               stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}

# Key pairs at their AIC-selected lags
key_pairs <- lag_sel[lag_sel$X %in% c("d_Reserves", "d_TGA") &
                      lag_sel$Y %in% c("d_Aaa", "d_Baa_Aaa"), ]
gc_at_aic <- do.call(rbind, lapply(1:nrow(key_pairs), function(i) {
  run_granger_at_lag(df_post, key_pairs$X[i], key_pairs$Y[i],
                      key_pairs$AIC_lag[i])
}))
cat("\nGranger at AIC-selected lags:\n")
print(gc_at_aic)
write.csv(gc_at_aic, "results_granger_aic_lag.csv", row.names = FALSE)


###############################################################################
#  2-3. RESIDUAL DIAGNOSTICS (Breusch-Godfrey + ARCH-LM)
###############################################################################
cat("\n── 2-3. Residual Diagnostics ──────────────────────────────────────\n")

run_diagnostics <- function(data, x_name, y_name, lag, bg_order = 5,
                             arch_order = 5) {
  xy <- data %>% select(all_of(c(x_name, y_name))) %>% tidyr::drop_na()
  n <- nrow(xy)
  if (n < lag + 50) return(NULL)

  Y <- xy[[y_name]][(lag+1):n]
  X_mat <- data.frame(row.names = 1:length(Y))
  for (i in 1:lag) {
    X_mat[[paste0("Y_lag", i)]] <- xy[[y_name]][(lag+1-i):(n-i)]
    X_mat[[paste0("X_lag", i)]] <- xy[[x_name]][(lag+1-i):(n-i)]
  }
  X_mat$Y <- Y
  fit <- lm(Y ~ ., data = X_mat)

  # Breusch-Godfrey
  bg <- bgtest(fit, order = bg_order)

  # ARCH-LM (manual)
  resid2 <- residuals(fit)^2
  m <- length(resid2)
  arch_df <- data.frame(e2 = resid2[(arch_order+1):m])
  for (i in 1:arch_order) {
    arch_df[[paste0("e2_lag", i)]] <- resid2[(arch_order+1-i):(m-i)]
  }
  arch_fit <- lm(e2 ~ ., data = arch_df)
  arch_r2 <- summary(arch_fit)$r.squared
  arch_lm <- nrow(arch_df) * arch_r2
  arch_p <- 1 - pchisq(arch_lm, df = arch_order)

  data.frame(
    Direction = paste0(x_name, " -> ", y_name),
    Lag = lag,
    BG_stat = round(bg$statistic, 3),
    BG_p = round(bg$p.value, 4),
    BG_reject = ifelse(bg$p.value < 0.05, "Yes", "No"),
    ARCH_LM = round(arch_lm, 3),
    ARCH_p = round(arch_p, 4),
    ARCH_reject = ifelse(arch_p < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}

# Run on key significant regressions
diag_pairs <- list(
  c("d_Reserves", "d_Aaa", 1),
  c("d_Reserves", "d_Aaa", 5),
  c("d_Reserves", "d_Aaa", 10),
  c("d_Reserves", "d_Baa_Aaa", 1),
  c("d_Reserves", "d_Baa_Aaa", 5),
  c("d_TGA", "d_Baa_Aaa", 5),
  c("d_TGA", "d_Baa_Aaa", 10),
  c("SOFR_EFFR", "d_Baa", 1),
  c("SOFR_EFFR", "d_Baa", 5)
)

diag_results <- do.call(rbind, lapply(diag_pairs, function(p) {
  run_diagnostics(df_post, p[1], p[2], as.integer(p[3]))
}))

cat("Residual diagnostics:\n")
print(diag_results)
write.csv(diag_results, "results_residual_diagnostics.csv", row.names = FALSE)


###############################################################################
#  4. HAC-ROBUST (NEWEY-WEST) GRANGER F-TESTS
###############################################################################
cat("\n── 4. HAC-Robust Granger Tests ────────────────────────────────────\n")

run_granger_hac <- function(data, x_name, y_name, lag) {
  xy <- data %>% select(all_of(c(x_name, y_name))) %>% tidyr::drop_na()
  n <- nrow(xy)
  if (n < lag + 50) return(NULL)

  Y <- xy[[y_name]][(lag+1):n]
  X_mat <- data.frame(row.names = 1:length(Y))
  for (i in 1:lag) {
    X_mat[[paste0("Y_lag", i)]] <- xy[[y_name]][(lag+1-i):(n-i)]
    X_mat[[paste0("X_lag", i)]] <- xy[[x_name]][(lag+1-i):(n-i)]
  }
  X_mat$Y <- Y

  # Unrestricted
  fit_u <- lm(Y ~ ., data = X_mat)
  # Restricted (only Y lags)
  y_lag_cols <- grep("^Y_lag", names(X_mat), value = TRUE)
  fit_r <- lm(as.formula(paste("Y ~", paste(y_lag_cols, collapse = " + "))),
               data = X_mat)

  # Newey-West HAC covariance
  nw_vcov <- sandwich::NeweyWest(fit_u, lag = lag, prewhite = FALSE)

  # Wald test with HAC
  wald <- lmtest::waldtest(fit_r, fit_u, vcov = nw_vcov)
  f_val <- wald$F[2]; p_val <- wald$`Pr(>F)`[2]

  # Also get classical F for comparison
  gt <- grangertest(as.formula(paste(y_name, "~", x_name)),
                     order = lag, data = xy)
  f_class <- gt$F[2]; p_class <- gt$`Pr(>F)`[2]

  sig_hac <- ifelse(p_val < 0.01, "***", ifelse(p_val < 0.05, "**",
             ifelse(p_val < 0.10, "*", "")))
  sig_class <- ifelse(p_class < 0.01, "***", ifelse(p_class < 0.05, "**",
               ifelse(p_class < 0.10, "*", "")))

  data.frame(
    Direction = paste0(x_name, " -> ", y_name),
    Lag = lag,
    N = nrow(X_mat),
    F_classical = round(f_class, 3),
    p_classical = round(p_class, 4),
    Sig_classical = sig_class,
    F_HAC = round(f_val, 3),
    p_HAC = round(p_val, 4),
    Sig_HAC = sig_hac,
    stringsAsFactors = FALSE
  )
}

# Run on all key pairs
hac_pairs <- list(
  c("d_Reserves", "d_Aaa", 1),
  c("d_Reserves", "d_Aaa", 5),
  c("d_Reserves", "d_Aaa", 10),
  c("d_Reserves", "d_Baa_Aaa", 1),
  c("d_Reserves", "d_Baa_Aaa", 5),
  c("d_TGA", "d_Baa_Aaa", 5),
  c("d_TGA", "d_Baa_Aaa", 10),
  c("d_TGA", "d_Baa", 1),
  c("SOFR_EFFR", "d_Baa", 1),
  c("SOFR_EFFR", "d_Aaa", 5)
)

hac_results <- do.call(rbind, lapply(hac_pairs, function(p) {
  run_granger_hac(df_post, p[1], p[2], as.integer(p[3]))
}))

cat("HAC-robust vs classical Granger:\n")
print(hac_results)
write.csv(hac_results, "results_granger_hac.csv", row.names = FALSE)


###############################################################################
#  5. TODA-YAMAMOTO PROCEDURE
###############################################################################
cat("\n── 5. Toda-Yamamoto Procedure ─────────────────────────────────────\n")

run_toda_yamamoto <- function(data, x_name, y_name, lag, d_max = 1) {
  xy <- data %>% select(all_of(c(x_name, y_name))) %>% tidyr::drop_na()
  n <- nrow(xy)
  p_total <- lag + d_max
  if (n < p_total + 50) return(NULL)

  Y <- xy[[y_name]][(p_total+1):n]
  X_mat <- data.frame(row.names = 1:length(Y))
  for (i in 1:p_total) {
    X_mat[[paste0("Y_lag", i)]] <- xy[[y_name]][(p_total+1-i):(n-i)]
    X_mat[[paste0("X_lag", i)]] <- xy[[x_name]][(p_total+1-i):(n-i)]
  }
  X_mat$Y <- Y

  # Unrestricted: all p + d_max lags
  fit_u <- lm(Y ~ ., data = X_mat)

  # Restricted: drop only X_lag1..X_lag{p} (keep X_lag{p+1}..X_lag{p+d_max})
  x_test_cols <- paste0("X_lag", 1:lag)
  keep_cols <- setdiff(names(X_mat)[names(X_mat) != "Y"], x_test_cols)
  fit_r <- lm(as.formula(paste("Y ~", paste(keep_cols, collapse = " + "))),
               data = X_mat)

  wald <- anova(fit_r, fit_u)
  f_val <- wald$F[2]; p_val <- wald$`Pr(>F)`[2]

  # Chi-square version
  chi2 <- f_val * lag
  chi2_p <- 1 - pchisq(chi2, df = lag)

  sig <- ifelse(chi2_p < 0.01, "***", ifelse(chi2_p < 0.05, "**",
         ifelse(chi2_p < 0.10, "*", "")))

  data.frame(
    Direction = paste0(x_name, " -> ", y_name),
    p = lag, d_max = d_max, N = nrow(X_mat),
    Chi2 = round(chi2, 3), Chi2_p = round(chi2_p, 4),
    F_equiv = round(f_val, 3), F_p = round(p_val, 4),
    Sig = sig, stringsAsFactors = FALSE
  )
}

ty_pairs <- list(
  c("d_Reserves", "d_Aaa", 1),
  c("d_Reserves", "d_Aaa", 5),
  c("d_Reserves", "d_Aaa", 10),
  c("d_Reserves", "d_Baa_Aaa", 1),
  c("d_Reserves", "d_Baa_Aaa", 5),
  c("d_TGA", "d_Baa_Aaa", 5),
  c("d_TGA", "d_Baa_Aaa", 10)
)

ty_results <- do.call(rbind, lapply(ty_pairs, function(p) {
  run_toda_yamamoto(df_post, p[1], p[2], as.integer(p[3]))
}))

cat("Toda-Yamamoto results:\n")
print(ty_results)
write.csv(ty_results, "results_toda_yamamoto.csv", row.names = FALSE)


###############################################################################
#  6. COVID ROBUSTNESS
###############################################################################
cat("\n── 6. COVID Robustness ────────────────────────────────────────────\n")

df_post_nocovid <- df_post %>%
  filter(!(Date >= as.Date("2020-02-01") & Date <= as.Date("2020-05-31")))

cat(sprintf("Post-2018 excl. COVID (Feb-May 2020): %d obs (removed %d)\n",
            nrow(df_post_nocovid), nrow(df_post) - nrow(df_post_nocovid)))

run_granger_simple <- function(data, x_name, y_name, lag) {
  xy <- data %>% select(all_of(c(x_name, y_name))) %>% tidyr::drop_na()
  if (nrow(xy) < lag + 30) return(NULL)
  tryCatch({
    gt <- grangertest(as.formula(paste(y_name, "~", x_name)),
                       order = lag, data = xy)
    f_val <- gt$F[2]; p_val <- gt$`Pr(>F)`[2]
    sig <- ifelse(p_val < 0.01, "***", ifelse(p_val < 0.05, "**",
           ifelse(p_val < 0.10, "*", "")))
    data.frame(X = x_name, Y = y_name, Lag = lag,
               N = nrow(xy), F_stat = round(f_val, 3),
               p_value = round(p_val, 4), Sig = sig,
               stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}

covid_pairs <- list(
  c("d_Reserves", "d_Aaa", 1), c("d_Reserves", "d_Aaa", 5),
  c("d_Reserves", "d_Aaa", 10),
  c("d_Reserves", "d_Baa", 1), c("d_Reserves", "d_Baa", 5),
  c("d_Reserves", "d_Baa", 10),
  c("d_Reserves", "d_Baa_Aaa", 1), c("d_Reserves", "d_Baa_Aaa", 5),
  c("d_Reserves", "d_Baa_Aaa", 10),
  c("d_TGA", "d_Baa_Aaa", 1), c("d_TGA", "d_Baa_Aaa", 5),
  c("d_TGA", "d_Baa_Aaa", 10),
  c("d_TGA", "d_Baa", 1), c("d_TGA", "d_Baa", 5),
  c("d_TGA", "d_Baa", 10),
  c("SOFR_EFFR", "d_Baa", 1), c("SOFR_EFFR", "d_Baa", 5),
  c("SOFR_EFFR", "d_Baa", 10),
  c("SOFR_EFFR", "d_Aaa", 1), c("SOFR_EFFR", "d_Aaa", 5),
  c("SOFR_EFFR", "d_Aaa", 10)
)

# Full sample results for comparison
gc_full <- do.call(rbind, lapply(covid_pairs, function(p) {
  r <- run_granger_simple(df_post, p[1], p[2], as.integer(p[3]))
  if (!is.null(r)) r$Sample <- "Full"
  r
}))

gc_nocovid <- do.call(rbind, lapply(covid_pairs, function(p) {
  r <- run_granger_simple(df_post_nocovid, p[1], p[2], as.integer(p[3]))
  if (!is.null(r)) r$Sample <- "Excl. COVID"
  r
}))

covid_comparison <- rbind(gc_full, gc_nocovid)
cat("COVID robustness comparison:\n")
print(covid_comparison[covid_comparison$X == "d_Reserves" &
                         covid_comparison$Y == "d_Aaa", ])
write.csv(covid_comparison, "results_covid_robustness.csv", row.names = FALSE)


###############################################################################
#  7. ROLLING GRANGER — QUANTITY CHANNEL
###############################################################################
cat("\n── 7. Rolling Granger: Quantity Channel ───────────────────────────\n")

roll_granger <- function(data, x_name, y_name, lag, window = 400,
                          step = 15) {
  xy <- data %>%
    select(Date, all_of(c(x_name, y_name))) %>%
    tidyr::drop_na() %>%
    arrange(Date)
  n <- nrow(xy)
  if (n < window + lag) return(data.frame(Date=as.Date(character()),
                                            p_value=numeric()))
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

roll_res_aaa <- roll_granger(df_post, "d_Reserves", "d_Aaa",
                              lag = 5, window = 400, step = 15)
roll_res_aaa$Direction <- "Reserves -> Aaa (lag 5)"

roll_tga_diff <- roll_granger(df_post, "d_TGA", "d_Baa_Aaa",
                               lag = 5, window = 400, step = 15)
roll_tga_diff$Direction <- "TGA -> Baa-Aaa (lag 5)"

roll_qty <- rbind(roll_res_aaa, roll_tga_diff)

p_roll_qty <- ggplot(roll_qty, aes(x = Date, y = p_value,
                                     color = Direction)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black",
             linewidth = 0.5) +
  geom_vline(xintercept = as.Date("2020-03-11"), linetype = "dashed",
             color = "red", linewidth = 0.4) +
  annotate("text", x = as.Date("2020-04-15"), y = 0.98,
           label = "COVID", color = "red", size = 3, hjust = 0) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Reserves -> Aaa (lag 5)" = "#2E75B6",
                                  "TGA -> Baa-Aaa (lag 5)" = "#C55A11")) +
  labs(title = "Rolling Granger p-values: Quantity Channel (Post-2018)",
       subtitle = "400-day rolling window, lag 5; black dashed = 5% significance",
       x = NULL, y = "p-value", color = NULL) +
  theme_thesis

ggsave("fig_a6_rolling_granger_qty.pdf", p_roll_qty,
       width = 10, height = 5)
cat("Saved: fig_a6_rolling_granger_qty.pdf\n")


###############################################################################
#  8. GENERATE LaTeX TABLES
###############################################################################
cat("\n── 8. Generating LaTeX Tables ─────────────────────────────────────\n")

clean_name <- function(x) {
  x <- gsub("^d_Baa_Aaa$", "$\\\\Delta$(Baa--Aaa)", x)
  x <- gsub("d_Baa_Aaa", "$\\\\Delta$(Baa--Aaa)", x)
  x <- gsub("^d_ON_RRP$", "$\\\\Delta$ON RRP", x)
  x <- gsub("d_ON_RRP", "$\\\\Delta$ON RRP", x)
  x <- gsub("^d_Reserves$", "$\\\\Delta$Reserves", x)
  x <- gsub("d_Reserves", "$\\\\Delta$Reserves", x)
  x <- gsub("^SOFR_EFFR$", "SOFR--EFFR", x)
  x <- gsub("SOFR_EFFR", "SOFR--EFFR", x)
  x <- gsub("^TED_spread$", "TED", x)
  x <- gsub("^d_Baa$", "$\\\\Delta$Baa", x)
  x <- gsub("d_Baa", "$\\\\Delta$Baa", x)
  x <- gsub("^d_Aaa$", "$\\\\Delta$Aaa", x)
  x <- gsub("d_Aaa", "$\\\\Delta$Aaa", x)
  x <- gsub("^d_TGA$", "$\\\\Delta$TGA", x)
  x <- gsub("d_TGA", "$\\\\Delta$TGA", x)
  x <- gsub(" -> ", " $\\\\rightarrow$ ", x)
  x
}

# Table: Lag selection
lag_tbl <- lag_sel_all
lag_tbl$X <- clean_name(lag_tbl$X)
lag_tbl$Y <- clean_name(lag_tbl$Y)
colnames(lag_tbl) <- c("$X$", "$Y$", "AIC lag", "BIC lag", "$N$")
xt_lag <- xtable(lag_tbl,
  caption = "Information Criterion Lag Selection for Bivariate Granger Pairs (Post-2018)",
  label = "tab:lag_selection",
  digits = c(0, 0, 0, 0, 0, 0))
print(xt_lag, file = "table_lag_selection.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "H",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.82, comment = FALSE)

# Table: Residual diagnostics
diag_tbl <- diag_results
diag_tbl$Direction <- clean_name(diag_tbl$Direction)
colnames(diag_tbl) <- c("Direction", "Lag", "BG stat", "BG $p$",
                         "BG reject?", "ARCH-LM", "ARCH $p$", "ARCH reject?")
xt_diag <- xtable(diag_tbl,
  caption = "Residual Diagnostics: Breusch-Godfrey (order 5) and ARCH-LM (order 5) Tests",
  label = "tab:diagnostics",
  digits = c(0, 0, 0, 3, 4, 0, 3, 4, 0))
print(xt_diag, file = "table_residual_diagnostics.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "H",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78, comment = FALSE)

# Table: HAC-robust comparison
hac_tbl <- hac_results
hac_tbl$Direction <- clean_name(hac_tbl$Direction)
colnames(hac_tbl) <- c("Direction", "Lag", "$N$",
                         "$F$ (classical)", "$p$ (classical)", "",
                         "$F$ (HAC)", "$p$ (HAC)", " ")
xt_hac <- xtable(hac_tbl,
  caption = "Classical vs.\\ HAC-Robust (Newey-West) Granger F-Tests",
  label = "tab:hac",
  digits = c(0, 0, 0, 0, 3, 4, 0, 3, 4, 0))
print(xt_hac, file = "table_granger_hac.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "H",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.72, comment = FALSE)

# Table: Toda-Yamamoto
ty_tbl <- ty_results
ty_tbl$Direction <- clean_name(ty_tbl$Direction)
colnames(ty_tbl) <- c("Direction", "$p$", "$d_{\\max}$", "$N$",
                        "$\\chi^2$", "$\\chi^2$ $p$",
                        "$F$", "$F$ $p$", "")
xt_ty <- xtable(ty_tbl,
  caption = "Toda-Yamamoto Granger Non-Causality Tests ($d_{\\max} = 1$)",
  label = "tab:toda_yamamoto",
  digits = c(0, 0, 0, 0, 0, 3, 4, 3, 4, 0))
print(xt_ty, file = "table_toda_yamamoto.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "H",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.78, comment = FALSE)

# Table: COVID robustness
covid_tbl <- covid_comparison
covid_tbl$Direction <- clean_name(paste0(covid_tbl$X, " -> ", covid_tbl$Y))
covid_tbl <- covid_tbl[, c("Sample", "Direction", "Lag", "N",
                             "F_stat", "p_value", "Sig")]
colnames(covid_tbl) <- c("Sample", "Direction", "Lag", "$N$",
                           "$F$-stat", "$p$-value", "")
xt_covid <- xtable(covid_tbl,
  caption = "COVID Robustness: Granger Causality With and Without February--May 2020",
  label = "tab:covid",
  digits = c(0, 0, 0, 0, 0, 3, 4, 0))
print(xt_covid, file = "table_covid_robustness.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "H",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.70, comment = FALSE)

# Table: SOFR-EFFR only (Battery 1 — fixing table duplication)
fwd_all <- read.csv("results_granger_post2018.csv")
fwd_price <- fwd_all[fwd_all$Type == "Forward" & fwd_all$X == "SOFR_EFFR", ]
fwd_price$Direction <- paste0(clean_name(fwd_price$X), " $\\rightarrow$ ",
                               clean_name(fwd_price$Y))
fwd_price_tbl <- fwd_price[, c("Direction", "Lag", "N", "F_stat",
                                 "p_value", "Sig")]
colnames(fwd_price_tbl) <- c("Direction", "Lag", "$N$", "$F$-stat",
                               "$p$-value", "")
xt_price <- xtable(fwd_price_tbl,
  caption = "Battery 1: SOFR--EFFR $\\rightarrow$ Credit Spread Changes (Post-2018)",
  label = "tab:price_post2018",
  digits = c(0, 0, 0, 0, 3, 4, 0))
print(xt_price, file = "table_granger_price_post2018.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "H",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      scalebox = 0.82, comment = FALSE)

cat("\nAll tables and figures generated successfully.\n")
cat("================================================================\n")
