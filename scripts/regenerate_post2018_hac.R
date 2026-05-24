###############################################################################
#  Regenerate Post-2018 Granger tables with HAC-robust inference
#  Matches the format used in the Dollar analysis (Section 5.4)
#  Output: 3 LaTeX tables → Empirical Analysis/tables/
###############################################################################

library(readr)
library(dplyr)
library(lmtest)
library(sandwich)
library(xtable)

OUT_DIR <- "/Users/jordi/Downloads/University/TFG/Data/Thesis/Empirical Analysis/tables"

# ── Load & prepare ──────────────────────────────────────────────────────────

df <- read_csv("/Users/jordi/Downloads/University/TFG/Data/Thesis/data_master.csv",
               show_col_types = FALSE)
df$Date <- as.Date(df$Date)

df_post <- df %>%
  filter(Date >= as.Date("2018-04-03")) %>%
  select(Date, SOFR_EFFR, d_Baa, d_Aaa, d_Baa_Aaa,
         d_Reserves, d_TGA, d_ON_RRP) %>%
  tidyr::drop_na(SOFR_EFFR)

cat(sprintf("Post-2018 sample: %d obs (%s to %s)\n",
            nrow(df_post), min(df_post$Date), max(df_post$Date)))

# ── Granger function with classical + HAC ───────────────────────────────────

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

  fit_u <- lm(Y ~ ., data = X_mat)
  y_lag_cols <- grep("^Y_lag", names(X_mat), value = TRUE)
  fit_r <- lm(as.formula(paste("Y ~", paste(y_lag_cols, collapse = " + "))),
              data = X_mat)

  # Classical F-test
  wt_cl <- waldtest(fit_r, fit_u)
  f_cl <- wt_cl$F[2]; p_cl <- wt_cl$`Pr(>F)`[2]

  # HAC-robust F-test (Newey-West, same as dollar analysis)
  wt_hac <- waldtest(fit_r, fit_u, vcov = vcovHAC)
  f_hac <- wt_hac$F[2]; p_hac <- wt_hac$`Pr(>F)`[2]

  data.frame(
    X = x_name, Y = y_name, Lag = lag, N = nrow(X_mat),
    F_classical = f_cl, p_classical = p_cl,
    F_HAC = f_hac, p_HAC = p_hac,
    stringsAsFactors = FALSE
  )
}

sig_stars <- function(p) {
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*",
  ifelse(p < 0.10,  ".",
                     ""))))
}

# ── Label mapping ───────────────────────────────────────────────────────────

label_map <- c(
  "SOFR_EFFR"  = "SOFR--EFFR",
  "d_Baa"      = "$\\Delta$Baa",
  "d_Aaa"      = "$\\Delta$Aaa",
  "d_Baa_Aaa"  = "$\\Delta$(Baa--Aaa)",
  "d_Reserves" = "$\\Delta$Reserves",
  "d_TGA"      = "$\\Delta$TGA",
  "d_ON_RRP"   = "$\\Delta$ON RRP"
)

make_direction <- function(x, y) {
  paste0(label_map[x], " $\\rightarrow$ ", label_map[y])
}

# ── Define all pairs ───────────────────────────────────────────────────────

# Battery 1: Price channel (SOFR--EFFR only, for separate table)
price_pairs <- expand.grid(
  X = "SOFR_EFFR",
  Y = c("d_Baa", "d_Aaa", "d_Baa_Aaa"),
  Lag = c(1, 5, 10),
  stringsAsFactors = FALSE
)

# Battery 2: All forward (price + quantity)
fwd_x <- c("SOFR_EFFR", "d_Reserves", "d_TGA", "d_ON_RRP")
fwd_y <- c("d_Baa", "d_Aaa", "d_Baa_Aaa")
fwd_pairs <- expand.grid(X = fwd_x, Y = fwd_y, Lag = c(1, 5, 10),
                          stringsAsFactors = FALSE)
# Order: group by X then Y then Lag
fwd_pairs <- fwd_pairs %>%
  mutate(X_order = match(X, fwd_x), Y_order = match(Y, fwd_y)) %>%
  arrange(X_order, Y_order, Lag) %>%
  select(X, Y, Lag)

# Battery 3: Reverse causality
rev_x <- c("d_Baa", "d_Aaa")
rev_y_sofr <- "SOFR_EFFR"
rev_y_qty  <- c("d_Reserves", "d_TGA")

rev_pairs <- bind_rows(
  expand.grid(X = rev_x, Y = rev_y_sofr, Lag = c(1, 5, 10), stringsAsFactors = FALSE),
  expand.grid(X = rev_x, Y = rev_y_qty,  Lag = c(1, 5, 10), stringsAsFactors = FALSE)
) %>%
  mutate(X_order = match(X, rev_x),
         Y_order = match(Y, c("SOFR_EFFR", "d_Reserves", "d_TGA"))) %>%
  arrange(X_order, Y_order, Lag) %>%
  select(X, Y, Lag)

# ── Run all tests ──────────────────────────────────────────────────────────

cat("\n── Computing HAC-robust Granger tests for all post-2018 pairs ──\n")

run_batch <- function(pairs_df) {
  results <- do.call(rbind, lapply(1:nrow(pairs_df), function(i) {
    run_granger_hac(df_post, pairs_df$X[i], pairs_df$Y[i], pairs_df$Lag[i])
  }))
  results$Direction <- make_direction(results$X, results$Y)
  results$Sig_cl  <- sig_stars(results$p_classical)
  results$Sig_hac <- sig_stars(results$p_HAC)
  results
}

price_results <- run_batch(price_pairs)
fwd_results   <- run_batch(fwd_pairs)
rev_results   <- run_batch(rev_pairs)

cat("\n── Battery 1 (Price) ──\n")
print(price_results %>% select(Direction, Lag, N, F_classical, p_classical, F_HAC, p_HAC))

cat("\n── Battery 2 (Forward, all) ──\n")
print(fwd_results %>% select(Direction, Lag, N, F_classical, p_classical, F_HAC, p_HAC))

cat("\n── Battery 3 (Reverse) ──\n")
print(rev_results %>% select(Direction, Lag, N, F_classical, p_classical, F_HAC, p_HAC))

# ── Generate LaTeX tables ──────────────────────────────────────────────────

make_table <- function(results, caption, label, filename) {
  tbl <- results %>%
    select(Direction, Lag, N,
           F_classical, p_classical, Sig_cl,
           F_HAC, p_HAC, Sig_hac)

  tbl$F_classical <- round(tbl$F_classical, 3)
  tbl$p_classical <- round(tbl$p_classical, 4)
  tbl$F_HAC       <- round(tbl$F_HAC, 3)
  tbl$p_HAC       <- round(tbl$p_HAC, 4)

  colnames(tbl) <- c("Direction", "Lag", "$N$",
                      "$F$ (class.)", "$p$ (class.)", "",
                      "$F$ (HAC)", "$p$ (HAC)", "")

  xt <- xtable(tbl, caption = caption, label = label,
               digits = c(0, 0, 0, 0, 3, 4, 0, 3, 4, 0))

  outpath <- file.path(OUT_DIR, filename)
  print(xt, file = outpath,
        include.rownames = FALSE,
        sanitize.text.function = identity,
        booktabs = TRUE, floating = TRUE,
        table.placement = "H",
        comment = FALSE,
        scalebox = 0.72)

  cat(sprintf("Written: %s\n", outpath))
}

# Table 1: Battery 1 (Price channel only)
make_table(price_results,
           "Battery 1: SOFR--EFFR $\\rightarrow$ Credit Spread Changes (Post-2018)",
           "tab:price_post2018",
           "table_granger_price_post2018.tex")

# Table 2: Battery 2 (All forward: price + quantity)
make_table(fwd_results,
           "Granger Causality: Funding $\\rightarrow$ Credit Spreads (Post-2018)",
           "tab:granger_fwd_post2018",
           "table_granger_fwd_post2018.tex")

# Table 3: Battery 3 (Reverse)
make_table(rev_results,
           "Reverse Causality: Credit Spreads $\\rightarrow$ Funding (Post-2018)",
           "tab:granger_rev_post2018",
           "table_granger_rev_post2018.tex")

cat("\nDone. All three tables regenerated with HAC columns.\n")
