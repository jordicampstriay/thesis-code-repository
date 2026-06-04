
required_pkgs <- c("readr", "dplyr", "tidyr", "zoo", "tseries", "ggplot2",
                   "scales", "corrplot", "xtable", "lubridate", "gridExtra",
                   "grDevices", "knitr")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

WD <- "/Users/jordi/Downloads/University/TFG/Data/Thesis"
setwd(WD)

cat("================================================================\n")
cat("  STEP 1 — DATA PREPARATION\n")
cat("  Working directory:", WD, "\n")
cat("================================================================\n\n")

raw <- read_csv("raw_data_thesis.csv",
                na = c("", "NA", "-", "#VALUE!"),
                show_col_types = FALSE)

colnames(raw) <- c("Date", "SOFR", "GC_Repo", "EFFR", "TED_spread", "DXY",
                    "HYG", "LQD", "SHV", "HYG_SHV", "EMB",
                    "TGCR", "Liq_Swaps", "BGCR",
                    "ON_RRP", "TGA", "Fed_Assets",
                    "Reserves", "Baa_spread", "Aaa_spread",
                    "TLT", "IEF", "MOVE", "SOFR_EFFR", "Aaa_Baa")

raw$Date <- as.Date(raw$Date)
for (col in colnames(raw)[-1])
  raw[[col]] <- suppressWarnings(as.numeric(raw[[col]]))

raw <- raw %>% arrange(Date)

cat(sprintf("Raw data: %d rows, %d columns\n", nrow(raw), ncol(raw)))
cat(sprintf("Date range: %s to %s\n\n",
            min(raw$Date, na.rm = TRUE), max(raw$Date, na.rm = TRUE)))

# ── DATA AVAILABILITY AUDIT ──────────────────────────────────────────────────
avail_audit <- function(data, vars) {
  out <- data.frame(Variable = character(), First = character(),
                    Last = character(), N = integer(),
                    stringsAsFactors = FALSE)
  for (v in vars) {
    if (v %in% colnames(data)) {
      idx <- which(!is.na(data[[v]]))
      if (length(idx) > 0)
        out <- rbind(out, data.frame(
          Variable = v,
          First = as.character(data$Date[min(idx)]),
          Last  = as.character(data$Date[max(idx)]),
          N     = length(idx)))
    }
  }
  out
}

audit_vars <- c("SOFR", "GC_Repo", "EFFR", "SOFR_EFFR", "TED_spread",
                "Baa_spread", "Aaa_spread", "Aaa_Baa",
                "DXY", "HYG", "LQD", "SHV", "HYG_SHV", "EMB",
                "TLT", "IEF", "MOVE",
                "Reserves", "TGA", "Fed_Assets", "ON_RRP", "Liq_Swaps",
                "TGCR", "BGCR")

availability <- avail_audit(raw, audit_vars)
cat("=== DATA AVAILABILITY (raw) ===\n")
print(availability, row.names = FALSE)
cat("\n")

cat("================================================================\n")
cat("  SOFR PROXY ANALYSIS\n")
cat("================================================================\n\n")


n_sofr      <- sum(!is.na(raw$SOFR))
n_gc_repo   <- sum(!is.na(raw$GC_Repo))
n_sofr_effr <- sum(!is.na(raw$SOFR_EFFR))
n_gap       <- sum(is.na(raw$SOFR_EFFR) & !is.na(raw$EFFR))

cat(sprintf("  SOFR (actual):       %d obs  (Apr 2018 -- Apr 2026)\n", n_sofr))
cat(sprintf("  GC Repo proxy:       %d obs  (May 2003 -- Feb 2018)\n", n_gc_repo))
cat(sprintf("  SOFR-EFFR (combined):%d obs\n", n_sofr_effr))
cat(sprintf("  Gap (Mar 2018):      ~22 business days\n\n"))


both_avail <- raw %>%
  filter(!is.na(SOFR_EFFR) & !is.na(TED_spread))

if (nrow(both_avail) > 100) {
  rho_sofr_ted <- cor(both_avail$SOFR_EFFR, both_avail$TED_spread,
                      use = "complete.obs")
  cat(sprintf("  Correlation SOFR-EFFR vs TED spread: %.4f  (N=%d)\n",
              rho_sofr_ted, nrow(both_avail)))
  cat("  TED spread contains credit risk (LIBOR component) — not a pure\n")
  cat("  secured/unsecured funding measure. SOFR-EFFR (or proxy) is preferred\n")
  cat("  when available. TED spread discontinued Jan 2022 (LIBOR cessation).\n\n")
}


cat("  DECISION: Use SOFR-EFFR (with GC Repo proxy) as primary funding spread.\n")
cat("            Extended sample: May 2003 -- April 2026 (gap: March 2018).\n")
cat("            TED spread: supplementary robustness (2003--Jan 2022).\n\n")


cat("================================================================\n")
cat("  DATA CLEANING\n")
cat("================================================================\n\n")

df <- raw

# C1. Holiday zeros
holiday_mask <- !is.na(df$Baa_spread) & df$Baa_spread == 0 &
                !is.na(df$Aaa_spread) & df$Aaa_spread == 0
n_holidays <- sum(holiday_mask)

market_cols <- c("Baa_spread", "Aaa_spread", "Aaa_Baa",
                 "DXY", "HYG", "LQD", "SHV", "EMB", "TLT", "IEF",
                 "TED_spread", "MOVE", "HYG_SHV")
for (col in market_cols) {
  if (col %in% colnames(df)) df[[col]][holiday_mask] <- NA
}
cat(sprintf("  Holiday zeros set to NA: %d rows\n", n_holidays))

# C2. DXY = 0 on non-trading days
dxy_zero <- !is.na(df$DXY) & df$DXY == 0
n_dxy_zero <- sum(dxy_zero)
df$DXY[dxy_zero] <- NA
cat(sprintf("  Additional DXY zeros set to NA: %d rows\n", n_dxy_zero))


df$SOFR <- NULL; df$GC_Repo <- NULL; df$EFFR <- NULL
df$TGCR <- NULL; df$BGCR <- NULL
cat("  Dropped: SOFR, GC_Repo, EFFR (individual rate levels)\n")
cat("  Dropped: TGCR, BGCR (spread would have 97% zeros)\n")

cat(sprintf("\n  Remaining columns: %d\n", ncol(df)))
cat(sprintf("  %s\n\n", paste(colnames(df), collapse = ", ")))

weekly_cols <- c("Reserves", "TGA", "Fed_Assets", "Liq_Swaps")
for (col in weekly_cols) {
  n_before <- sum(is.na(df[[col]]))
  df[[col]] <- na.locf(df[[col]], na.rm = FALSE)
  n_after <- sum(is.na(df[[col]]))
  cat(sprintf("  Forward-filled %s: %d NA -> %d NA\n", col, n_before, n_after))
}
cat("\n")

cat("================================================================\n")
cat("  VARIABLE CONSTRUCTION\n")
cat("================================================================\n\n")

df$Baa_Aaa <- df$Baa_spread - df$Aaa_spread
df$Aaa_Baa <- NULL
cat("  Created: Baa_Aaa = Baa_spread - Aaa_spread\n")

df <- df %>% mutate(
  d_Reserves   = c(NA, diff(Reserves)),
  d_TGA        = c(NA, diff(TGA)),
  d_ON_RRP     = c(NA, diff(ON_RRP)),
  d_Fed_Assets = c(NA, diff(Fed_Assets))
)
cat("  Created: d_Reserves, d_TGA, d_ON_RRP, d_Fed_Assets\n")

df <- df %>% mutate(
  d_Baa     = c(NA, diff(Baa_spread)),
  d_Aaa     = c(NA, diff(Aaa_spread)),
  d_Baa_Aaa = c(NA, diff(Baa_Aaa))
)
cat("  Created: d_Baa, d_Aaa, d_Baa_Aaa\n")

df <- df %>% mutate(
  lr_HYG = c(NA, diff(log(HYG))) * 100,
  lr_LQD = c(NA, diff(log(LQD))) * 100,
  lr_DXY = c(NA, diff(log(DXY))) * 100,
  lr_TLT = c(NA, diff(log(TLT))) * 100,
  lr_IEF = c(NA, diff(log(IEF))) * 100,
  lr_SHV = c(NA, diff(log(SHV))) * 100,
  lr_EMB = c(NA, diff(log(EMB))) * 100
)
cat("  Created: lr_HYG, lr_LQD, lr_DXY, lr_TLT, lr_IEF, lr_SHV, lr_EMB\n")

df <- df %>% mutate(
  log_HYG_LQD   = log(HYG / LQD),
  d_log_HYG_LQD = c(NA, diff(log_HYG_LQD)),
  log_HYG_SHV   = log(HYG_SHV),
  d_log_HYG_SHV = c(NA, diff(log_HYG_SHV))
)
cat("  Created: log_HYG_LQD, d_log_HYG_LQD, log_HYG_SHV, d_log_HYG_SHV\n\n")


cat("================================================================\n")
cat("  SAMPLE DEFINITIONS\n")
cat("================================================================\n\n")

df_extended <- df %>% filter(Date >= as.Date("2003-05-01"))

df_post <- df %>% filter(Date >= as.Date("2018-04-02"))

df_ted <- df %>% filter(Date >= as.Date("2003-05-01"),
                         Date <= as.Date("2022-01-21"))

cat(sprintf("  Extended sample (SOFR-EFFR w/ proxy): %s to %s  (%d rows)\n",
            min(df_extended$Date), max(df_extended$Date), nrow(df_extended)))
cat(sprintf("    Valid SOFR-EFFR: %d\n", sum(!is.na(df_extended$SOFR_EFFR))))
cat(sprintf("  Post-2018 sample (actual SOFR):       %s to %s  (%d rows)\n",
            min(df_post$Date), max(df_post$Date), nrow(df_post)))
cat(sprintf("    Valid SOFR-EFFR: %d\n", sum(!is.na(df_post$SOFR_EFFR))))
cat(sprintf("  TED spread sample:                    %s to %s  (%d rows)\n",
            min(df_ted$Date), max(df_ted$Date), nrow(df_ted)))
cat(sprintf("    Valid TED_spread: %d\n\n", sum(!is.na(df_ted$TED_spread))))


write.csv(df, "data_master.csv", row.names = FALSE)
cat("  Saved: data_master.csv\n\n")

ccb <- read_csv("cross_currency_basis.csv",
                col_names = c("Date", "AUD", "CHF", "EUR", "GBP", "JPY"),
                skip = 1, show_col_types = FALSE)
ccb$Date <- as.Date(ccb$Date)
ccb <- ccb %>% arrange(Date)
cat(sprintf("  Cross-currency basis: %s to %s (%d obs)\n\n",
            min(ccb$Date), max(ccb$Date), nrow(ccb)))

theme_thesis <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(face = "bold", size = 13, margin = margin(b = 4)),
      plot.subtitle    = element_text(color = "grey40", size = 9, margin = margin(b = 8)),
      plot.caption     = element_text(size = 7.5, color = "grey50", hjust = 0, margin = margin(t = 8)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      axis.title       = element_text(size = 10, color = "grey30"),
      axis.text        = element_text(size = 9),
      legend.position  = "bottom",
      legend.text      = element_text(size = 9),
      plot.margin      = margin(12, 14, 8, 10)
    )
}

COL_NAVY   <- "#1B4F72"
COL_RED    <- "#922B21"
COL_GREEN  <- "#1E8449"
COL_ORANGE <- "#D4AC0D"
COL_PURPLE <- "#7B1FA2"
COL_GREY   <- "grey60"

rec_bands <- data.frame(
  start = as.Date(c("2007-12-01", "2020-02-01")),
  end   = as.Date(c("2009-06-01", "2020-04-01"))
)

stress_bands <- data.frame(
  start = as.Date(c("2018-12-01", "2019-09-16", "2020-03-09", "2022-09-01")),
  end   = as.Date(c("2018-12-31", "2019-10-15", "2020-04-15", "2022-12-31")),
  label = c("Q4 2018", "Sept 2019\nrepo", "COVID-19", "QT stress")
)

add_rec <- function(p, from = as.Date("2003-01-01")) {
  bands <- rec_bands[rec_bands$end >= from, ]
  if (nrow(bands) > 0)
    p <- p + geom_rect(data = bands,
                        aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                        fill = "grey80", alpha = 0.4, inherit.aes = FALSE)
  p
}

add_stress <- function(p) {
  p + geom_rect(data = stress_bands,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                fill = "#FFCCCC", alpha = 0.35, inherit.aes = FALSE)
}


cat("================================================================\n")
cat("  TABLE: DATA SOURCES\n")
cat("================================================================\n\n")

source_table <- data.frame(
  Variable = c(
    "SOFR--EFFR spread",
    "GC Repo--EFFR spread (proxy)",
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
    "SHV (iShares Short Treasury ETF)",
    "HYG/SHV ratio",
    "EMB (iShares EM Bond ETF)",
    "MOVE index",
    "Cross-currency basis (3M)"
  ),
  Source = c(
    "FRBNY; FRED (SOFR, EFFR)",
    "FRBNY Primary Dealer Survey",
    "FRED (TEDRATE)",
    "FRED (BAA10Y)",
    "FRED (AAA10Y)",
    "Constructed",
    "FRED (WRESBAL)",
    "FRED (WTREGEN)",
    "FRED (WALCL)",
    "FRED (RRPONTSYD)",
    "FRED (SWPT)",
    "FRED (DTWEXBGS)",
    "Bloomberg",
    "Bloomberg",
    "Bloomberg",
    "Constructed",
    "Bloomberg",
    "Bloomberg",
    "Bloomberg"
  ),
  Frequency = c(
    "Daily", "Daily", "Daily", "Daily", "Daily", "Daily",
    "Weekly (Wed)", "Weekly (Wed)", "Weekly (Wed)",
    "Daily", "Weekly (Wed)",
    "Daily", "Daily", "Daily", "Daily", "Daily", "Daily", "Daily", "Daily"
  ),
  Available = c(
    "Apr 2018--Apr 2026",
    "May 2003--Feb 2018",
    "May 2003--Jan 2022",
    "May 2003--Apr 2026",
    "May 2003--Apr 2026",
    "May 2003--Apr 2026",
    "May 2003--Apr 2026",
    "May 2003--Apr 2026",
    "May 2003--Apr 2026",
    "May 2003--Apr 2026",
    "May 2003--Apr 2026",
    "Jan 2006--Apr 2026",
    "Apr 2007--Apr 2026",
    "Jan 2005--Apr 2026",
    "Jan 2007--Apr 2026",
    "Apr 2007--Apr 2026",
    "Dec 2007--Apr 2026",
    "Mar 2021--Apr 2026",
    "May 2024--May 2026"
  ),
  Role = c(
    "Funding input (price)",
    "Funding input (proxy)",
    "Funding input (robustness)",
    "Credit output",
    "Credit output",
    "Credit output (collateral)",
    "Funding input (quantity)",
    "Funding input (exogenous)",
    "Funding input (quantity)",
    "Funding input (quantity)",
    "Descriptive only",
    "Risk appetite / global",
    "Risk appetite extension",
    "Risk appetite extension",
    "Risk appetite extension",
    "Risk appetite extension",
    "EM transmission",
    "Validation",
    "Descriptive (post-2024)"
  ),
  stringsAsFactors = FALSE
)

source_xt <- xtable(source_table,
  caption = "Data Sources and Variable Descriptions",
  label = "tab:data_sources",
  align = c("l", "l", "l", "c", "c", "l"))
print(source_xt, file = "table_data_sources.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top",
      sanitize.text.function = function(x) gsub("--", "\\\\textendash{}", x),
      table.placement = "htbp", scalebox = 0.72)
cat("  Saved: table_data_sources.tex\n\n")


avail_clean_vars <- c("SOFR_EFFR", "TED_spread",
                       "Baa_spread", "Aaa_spread", "Baa_Aaa",
                       "DXY", "HYG", "LQD", "SHV", "EMB",
                       "Reserves", "TGA", "Fed_Assets", "ON_RRP",
                       "Liq_Swaps", "MOVE", "HYG_SHV")
avail_clean <- avail_audit(df, avail_clean_vars)

avail_xt <- xtable(avail_clean,
  caption = "Data Availability After Cleaning",
  label = "tab:data_availability", digits = c(0, 0, 0, 0, 0))
print(avail_xt, file = "table_data_availability.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = function(x) gsub("_", "\\\\_", x))
cat("  Saved: table_data_availability.tex\n\n")



cat("================================================================\n")
cat("  SUMMARY STATISTICS\n")
cat("================================================================\n\n")

compute_summary <- function(data, vars) {
  out <- data.frame(Variable = character(), N = integer(),
                    Mean = numeric(), SD = numeric(),
                    Min = numeric(), Q25 = numeric(),
                    Median = numeric(), Q75 = numeric(),
                    Max = numeric(), stringsAsFactors = FALSE)
  for (v in vars) {
    if (v %in% colnames(data)) {
      x <- na.omit(data[[v]])
      if (length(x) > 0) {
        q <- quantile(x, c(0.25, 0.50, 0.75))
        out <- rbind(out, data.frame(
          Variable = v, N = length(x),
          Mean = round(mean(x), 4), SD = round(sd(x), 4),
          Min = round(min(x), 4), Q25 = round(q[1], 4),
          Median = round(q[2], 4), Q75 = round(q[3], 4),
          Max = round(max(x), 4)))
      }
    }
  }
  out
}

summ_vars <- c("SOFR_EFFR", "TED_spread",
               "Baa_spread", "Aaa_spread", "Baa_Aaa",
               "DXY", "HYG", "LQD", "SHV", "EMB", "MOVE",
               "Reserves", "TGA", "Fed_Assets", "ON_RRP")

summ_ext <- compute_summary(df_extended, summ_vars)
cat("--- Extended Sample ---\n")
print(summ_ext, row.names = FALSE)

summ_post <- compute_summary(df_post, summ_vars)
cat("\n--- Post-2018 Sample ---\n")
print(summ_post, row.names = FALSE)

for (s in list(
  list(tbl = summ_ext,  f = "table_summary_extended.tex",
       cap = "Summary Statistics --- Extended Sample (May 2003 -- April 2026)",
       lab = "tab:summary_ext"),
  list(tbl = summ_post, f = "table_summary_post2018.tex",
       cap = "Summary Statistics --- Post-2018 Sample (April 2018 -- April 2026)",
       lab = "tab:summary_post2018")
)) {
  xt <- xtable(s$tbl, caption = s$cap, label = s$lab,
               digits = c(0, 0, 0, 4, 4, 4, 4, 4, 4, 4))
  print(xt, file = s$f, include.rownames = FALSE, booktabs = TRUE,
        caption.placement = "top", table.placement = "htbp",
        sanitize.text.function = function(x) gsub("_", "\\\\_", x),
        scalebox = 0.78)
  cat(sprintf("  Saved: %s\n", s$f))
}
cat("\n")


cat("================================================================\n")
cat("  ADF UNIT ROOT TESTS\n")
cat("================================================================\n\n")

adf_vars <- c("SOFR_EFFR", "TED_spread",
               "Baa_spread", "Aaa_spread", "Baa_Aaa",
               "DXY", "HYG", "LQD", "SHV", "EMB",
               "Reserves", "TGA", "Fed_Assets", "ON_RRP", "MOVE")

adf_results <- data.frame(
  Variable = character(), Sample = character(), N = integer(),
  ADF_level = numeric(), p_level = numeric(), Decision_level = character(),
  ADF_diff = numeric(), p_diff = numeric(), Decision_diff = character(),
  stringsAsFactors = FALSE
)

for (v in adf_vars) {
  if (v == "MOVE") {
    x_src <- df_post
    samp <- "Post-2018"
  } else if (v == "SOFR_EFFR") {
    x_src <- df_extended
    samp <- "Extended"
  } else if (v == "TED_spread") {
    x_src <- df_ted
    samp <- "To Jan 2022"
  } else {
    x_src <- df_extended
    samp <- "Extended"
  }

  x <- x_src[[v]]
  x <- x[!is.na(x)]
  if (length(x) < 50) next

  tryCatch({
    adf_lev <- adf.test(x, alternative = "stationary")
    dx <- diff(x); dx <- dx[!is.na(dx)]
    adf_dif <- adf.test(dx, alternative = "stationary")

    adf_results <- rbind(adf_results, data.frame(
      Variable = v, Sample = samp, N = length(x),
      ADF_level  = round(adf_lev$statistic, 3),
      p_level    = round(adf_lev$p.value, 4),
      Decision_level = ifelse(adf_lev$p.value < 0.05, "Stationary", "Unit root"),
      ADF_diff   = round(adf_dif$statistic, 3),
      p_diff     = round(adf_dif$p.value, 4),
      Decision_diff  = ifelse(adf_dif$p.value < 0.05, "Stationary", "Unit root")))
  }, error = function(e) cat(sprintf("  ADF error %s: %s\n", v, e$message)))
}

print(adf_results, row.names = FALSE)

cat("\n  Stationary in levels: ",
    paste(adf_results$Variable[adf_results$Decision_level == "Stationary"],
          collapse = ", "), "\n")
cat("  Unit root in levels:  ",
    paste(adf_results$Variable[adf_results$Decision_level == "Unit root"],
          collapse = ", "), "\n")
cat("  All unit-root vars stationary after differencing: ",
    all(adf_results$Decision_diff[adf_results$Decision_level == "Unit root"] == "Stationary"),
    "\n\n")

adf_latex <- adf_results %>% select(-Sample)
colnames(adf_latex) <- c("Variable", "$N$", "ADF (Level)", "$p$-value",
                          "Decision", "ADF ($\\Delta$)", "$p$-value ", "Decision ")
adf_xt <- xtable(adf_latex,
  caption = "Augmented Dickey--Fuller Unit Root Tests",
  label = "tab:adf_tests", digits = c(0, 0, 0, 3, 4, 0, 3, 4, 0))
print(adf_xt, file = "table_adf_results.tex",
      include.rownames = FALSE, booktabs = TRUE,
      caption.placement = "top", table.placement = "htbp",
      sanitize.colnames.function = identity,
      sanitize.text.function = function(x) gsub("_", "\\\\_", x),
      scalebox = 0.80)
cat("  Saved: table_adf_results.tex\n\n")

cat("================================================================\n")
cat("  CORRELATION ANALYSIS\n")
cat("================================================================\n\n")

corr_fd_vars <- c("SOFR_EFFR", "d_Reserves", "d_TGA", "d_ON_RRP",
                   "d_Baa", "d_Aaa", "d_Baa_Aaa")
corr_fd <- df_post %>% select(all_of(corr_fd_vars)) %>% drop_na()
corr_fd_mat <- cor(corr_fd)

cat("--- Post-2018 First Differences ---\n")
print(round(corr_fd_mat, 3))

cat(sprintf("\n  SOFR_EFFR vs d_Baa: %.4f\n", corr_fd_mat["SOFR_EFFR", "d_Baa"]))
cat(sprintf("  SOFR_EFFR vs d_Aaa: %.4f\n", corr_fd_mat["SOFR_EFFR", "d_Aaa"]))
cat(sprintf("  d_Reserves vs d_TGA: %.4f\n", corr_fd_mat["d_Reserves", "d_TGA"]))
cat("  >> Near-zero contemporaneous correlations: transmission is lagged.\n\n")

corr_ext_vars <- c("SOFR_EFFR", "d_Reserves", "d_TGA",
                    "d_Baa", "d_Aaa", "d_Baa_Aaa")
corr_ext <- df_extended %>% select(all_of(corr_ext_vars)) %>% drop_na()
corr_ext_mat <- cor(corr_ext)

cat("--- Extended Sample (with SOFR proxy) ---\n")
print(round(corr_ext_mat, 3))
cat("\n")

corr_ted_vars <- c("TED_spread", "d_Reserves", "d_TGA",
                    "d_Baa", "d_Aaa", "d_Baa_Aaa")
corr_ted <- df_ted %>% select(all_of(corr_ted_vars)) %>% drop_na()
corr_ted_mat <- cor(corr_ted)

cat("--- TED Spread Sample ---\n")
print(round(corr_ted_mat, 3))
cat("\n")

corr_fd_nice <- round(corr_fd_mat, 3)
rownames(corr_fd_nice) <- colnames(corr_fd_nice) <- c(
  "SOFR-EFFR", "$\\Delta$Reserves", "$\\Delta$TGA", "$\\Delta$ON RRP",
  "$\\Delta$Baa", "$\\Delta$Aaa", "$\\Delta$(Baa-Aaa)")

corr_xt <- xtable(corr_fd_nice,
  caption = "Correlation Matrix --- First Differences (Post-2018 Sample)",
  label = "tab:corr_fd", digits = 3)
print(corr_xt, file = "table_correlation_fd.tex",
      booktabs = TRUE, caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      sanitize.rownames.function = identity, scalebox = 0.75)
cat("  Saved: table_correlation_fd.tex\n\n")

corr_ext_nice <- round(corr_ext_mat, 3)
rownames(corr_ext_nice) <- colnames(corr_ext_nice) <- c(
  "SOFR-EFFR (proxy)", "$\\Delta$Reserves", "$\\Delta$TGA",
  "$\\Delta$Baa", "$\\Delta$Aaa", "$\\Delta$(Baa-Aaa)")

corr_ext_xt <- xtable(corr_ext_nice,
  caption = "Correlation Matrix --- First Differences (Extended Sample with SOFR Proxy)",
  label = "tab:corr_ext", digits = 3)
print(corr_ext_xt, file = "table_correlation_extended.tex",
      booktabs = TRUE, caption.placement = "top", table.placement = "htbp",
      sanitize.text.function = identity,
      sanitize.colnames.function = identity,
      sanitize.rownames.function = identity, scalebox = 0.78)
cat("  Saved: table_correlation_extended.tex\n\n")


cat("================================================================\n")
cat("  MULTICOLLINEARITY DIAGNOSTIC\n")
cat("================================================================\n\n")

mc_vars <- c("SOFR_EFFR", "d_Reserves", "d_TGA")
mc_data <- df_post %>% select(all_of(mc_vars)) %>% drop_na()
mc_cor <- cor(mc_data)
cat("VAR funding inputs (post-2018):\n")
print(round(mc_cor, 4))
max_off <- max(abs(mc_cor[upper.tri(mc_cor)]))
cat(sprintf("\n  Max off-diagonal |r|: %.4f — %s\n\n",
            max_off, ifelse(max_off < 0.50, "no concern", "INVESTIGATE")))

cat("================================================================\n")
cat("  SOFR-EFFR DISTRIBUTION\n")
cat("================================================================\n\n")

se_clean <- na.omit(df_post$SOFR_EFFR)
n_near_zero <- sum(abs(se_clean) < 0.005)
cat(sprintf("  Total obs: %d\n", length(se_clean)))
cat(sprintf("  Near-zero (|x| < 0.5bp): %d (%.1f%%)\n",
            n_near_zero, n_near_zero / length(se_clean) * 100))
cat(sprintf("  Mean: %.4f | SD: %.4f | Skewness: %.2f\n",
            mean(se_clean), sd(se_clean),
            mean((se_clean - mean(se_clean))^3) / sd(se_clean)^3))
cat("  >> Motivates regime-dependent analysis (Section 5.5)\n\n")


cat("================================================================\n")
cat("  HEAT MAP FIGURES\n")
cat("================================================================\n\n")

heat_cols <- colorRampPalette(c("#B71C1C", "#FFCDD2", "white",
                                 "#BBDEFB", "#0D47A1"))(200)

pdf("fig_corr_fd_post2018.pdf", width = 8, height = 7)
par(mar = c(1, 1, 3, 1))
corrplot(corr_fd_mat, method = "color", type = "upper",
         addCoef.col = "black", number.cex = 0.75,
         tl.col = "black", tl.srt = 45, tl.cex = 0.85,
         col = heat_cols, title = "")
mtext("Correlation: First Differences (Post-2018)", side = 3, line = 1.5,
      cex = 1.1, font = 2)
dev.off()
cat("  Saved: fig_corr_fd_post2018.pdf\n")

pdf("fig_corr_fd_extended.pdf", width = 7.5, height = 6.5)
par(mar = c(1, 1, 3, 1))
corrplot(corr_ext_mat, method = "color", type = "upper",
         addCoef.col = "black", number.cex = 0.8,
         tl.col = "black", tl.srt = 45, tl.cex = 0.9,
         col = heat_cols, title = "")
mtext("Correlation: First Differences (Extended Sample)", side = 3, line = 1.5,
      cex = 1.1, font = 2)
dev.off()
cat("  Saved: fig_corr_fd_extended.pdf\n\n")


cat("================================================================\n")
cat("  TIME SERIES PLOTS\n")
cat("================================================================\n\n")

p1 <- ggplot(df_extended %>% filter(!is.na(SOFR_EFFR)),
             aes(x = Date, y = SOFR_EFFR * 100))
p1 <- add_rec(p1)
p1 <- p1 +
  geom_hline(yintercept = 0, color = COL_GREY, linewidth = 0.4, linetype = "dashed") +
  geom_vline(xintercept = as.Date("2018-04-02"), color = COL_ORANGE,
             linewidth = 0.5, linetype = "dotted") +
  annotate("text", x = as.Date("2018-06-01"), y = 250,
           label = "SOFR\npublication\nbegins", size = 2.8, color = COL_ORANGE,
           hjust = 0) +
  geom_line(color = COL_NAVY, linewidth = 0.3, alpha = 0.85) +
  labs(title = "SOFR-EFFR Spread (Extended Sample with GC Repo Proxy)",
       subtitle = "May 2003 -- April 2026. Pre-2018: GC Repo Rate minus EFFR. Shaded: NBER recessions.",
       x = NULL, y = "Basis points",
       caption = "Source: FRBNY, FRED. Proxy: Fed FEDS Notes (2019).") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0.01, 0)) +
  theme_thesis()

p2 <- ggplot(df_post %>% filter(!is.na(SOFR_EFFR)),
             aes(x = Date, y = SOFR_EFFR * 100))
p2 <- add_stress(p2)
p2 <- p2 +
  geom_hline(yintercept = 0, color = COL_GREY, linewidth = 0.4, linetype = "dashed") +
  geom_line(color = COL_NAVY, linewidth = 0.35, alpha = 0.85) +
  labs(title = "SOFR-EFFR Spread (Post-2018)",
       subtitle = "Shaded: stress episodes.",
       x = NULL, y = "Basis points",
       caption = "Source: FRBNY.") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0.01, 0)) +
  theme_thesis()

p3_data <- df_extended %>%
  select(Date, Baa_spread, Aaa_spread) %>%
  pivot_longer(-Date, names_to = "S", values_to = "V") %>% filter(!is.na(V))
p3 <- ggplot(p3_data, aes(x = Date, y = V, color = S))
p3 <- add_rec(p3)
p3 <- p3 +
  geom_line(linewidth = 0.3, alpha = 0.85) +
  scale_color_manual(values = c("Aaa_spread" = COL_NAVY, "Baa_spread" = COL_RED),
                     labels = c("Moody's Aaa", "Moody's Baa")) +
  labs(title = "Corporate Credit Spreads",
       subtitle = "Moody's Baa and Aaa over 10-year Treasury. Shaded: NBER recessions.",
       x = NULL, y = "Spread (pp)", color = NULL,
       caption = "Source: FRED (BAA10Y, AAA10Y).") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0.01, 0)) +
  theme_thesis()

p4 <- ggplot(df_extended %>% filter(!is.na(Baa_Aaa)),
             aes(x = Date, y = Baa_Aaa))
p4 <- add_rec(p4)
p4 <- p4 +
  geom_line(color = COL_GREEN, linewidth = 0.3, alpha = 0.85) +
  labs(title = "Baa -- Aaa Quality Differential",
       subtitle = "Positive = Baa wider than Aaa (normal). Higher = greater credit differentiation.",
       x = NULL, y = "Differential (pp)",
       caption = "Source: Constructed from FRED.") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0.01, 0)) +
  theme_thesis()

p5_data <- df_extended %>%
  filter(!is.na(Reserves)) %>%
  mutate(Reserves = Reserves / 1e6, TGA = TGA / 1e6) %>%
  select(Date, Reserves, TGA) %>%
  pivot_longer(-Date, names_to = "S", values_to = "V") %>% filter(!is.na(V))
p5 <- ggplot(p5_data, aes(x = Date, y = V, color = S))
p5 <- add_rec(p5)
p5 <- p5 +
  geom_line(linewidth = 0.4, alpha = 0.85) +
  scale_color_manual(values = c("Reserves" = COL_NAVY, "TGA" = COL_ORANGE)) +
  labs(title = "Reserve Balances and Treasury General Account",
       subtitle = "Weekly (forward-filled). Trillions USD.",
       x = NULL, y = "USD trillions", color = NULL,
       caption = "Source: FRED (WRESBAL, WTREGEN).") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0.01, 0)) +
  theme_thesis()

p6 <- ggplot(df_ted %>% filter(!is.na(TED_spread)),
             aes(x = Date, y = TED_spread * 100))
p6 <- add_rec(p6)
p6 <- p6 +
  geom_line(color = COL_RED, linewidth = 0.3, alpha = 0.85) +
  labs(title = "TED Spread",
       subtitle = "3M LIBOR minus 3M T-bill. Discontinued January 2022.",
       x = NULL, y = "Basis points",
       caption = "Source: FRED (TEDRATE).") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0.01, 0)) +
  theme_thesis()

p7_data <- df_extended %>%
  select(Date, HYG, LQD) %>%
  pivot_longer(-Date, names_to = "S", values_to = "V") %>% filter(!is.na(V))
p7 <- ggplot(p7_data, aes(x = Date, y = V, color = S))
p7 <- add_rec(p7)
p7 <- p7 +
  geom_line(linewidth = 0.3, alpha = 0.85) +
  scale_color_manual(values = c("HYG" = COL_RED, "LQD" = COL_NAVY)) +
  labs(title = "Bond ETF Prices: HYG and LQD",
       subtitle = "HYG = high yield, LQD = investment grade.",
       x = NULL, y = "Price (USD)", color = NULL,
       caption = "Source: Bloomberg.") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0.01, 0)) +
  theme_thesis()

p8 <- ggplot(df_extended %>% filter(!is.na(ON_RRP), Date >= as.Date("2013-01-01")),
             aes(x = Date, y = ON_RRP / 1e3)) +
  geom_line(color = COL_PURPLE, linewidth = 0.4, alpha = 0.85) +
  labs(title = "Overnight Reverse Repo Facility",
       subtitle = "Daily take-up, billions USD.",
       x = NULL, y = "USD billions",
       caption = "Source: FRED (RRPONTSYD).") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0.01, 0)) +
  theme_thesis()

p9 <- ggplot(df_extended %>% filter(!is.na(DXY)),
             aes(x = Date, y = DXY))
p9 <- add_rec(p9)
p9 <- p9 +
  geom_line(color = COL_GREEN, linewidth = 0.3, alpha = 0.85) +
  labs(title = "US Dollar Index (DXY)",
       subtitle = "Higher = stronger dollar = tighter global USD funding.",
       x = NULL, y = "Index",
       caption = "Source: FRED (DTWEXBGS).") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0.01, 0)) +
  theme_thesis()

p10_data <- df_post %>% filter(!is.na(SOFR_EFFR) & !is.na(Baa_spread))
p10 <- ggplot(p10_data)
p10 <- add_stress(p10)
p10 <- p10 +
  geom_line(aes(x = Date, y = SOFR_EFFR * 100, color = "SOFR-EFFR (bps)"),
            linewidth = 0.4, alpha = 0.85) +
  geom_line(aes(x = Date,
                y = (Baa_spread - mean(Baa_spread, na.rm = TRUE)) * 30,
                color = "Baa spread (demeaned, rescaled)"),
            linewidth = 0.4, alpha = 0.7) +
  scale_color_manual(values = c("SOFR-EFFR (bps)" = COL_NAVY,
                                 "Baa spread (demeaned, rescaled)" = COL_RED)) +
  labs(title = "Funding Spread vs Credit Spread (Post-2018)",
       subtitle = "Visual inspection for co-movement during stress.",
       x = NULL, y = "SOFR-EFFR (bps)", color = NULL,
       caption = "Baa spread demeaned and rescaled for visual comparison only.") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0.01, 0)) +
  theme_thesis()

# ── P11: Cross-currency basis (JPY, descriptive) ──
p11 <- ggplot(ccb %>% filter(!is.na(JPY)), aes(x = Date, y = JPY)) +
  geom_hline(yintercept = 0, color = COL_GREY, linewidth = 0.4, linetype = "dashed") +
  geom_line(color = COL_RED, linewidth = 0.5, alpha = 0.85) +
  labs(title = "3-Month JPY/USD Cross-Currency Basis",
       subtitle = "May 2024 -- May 2026. Negative = dollar funding premium.",
       x = NULL, y = "Basis points",
       caption = "Source: Bloomberg.") +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y", expand = c(0.01, 0)) +
  theme_thesis()

# Save all plots
plot_list <- list(
  list(p1,  "fig_ts_sofr_effr_extended.pdf", 10, 4.5),
  list(p2,  "fig_ts_sofr_effr_post2018.pdf", 10, 4.5),
  list(p3,  "fig_ts_credit_spreads.pdf",     10, 4.5),
  list(p4,  "fig_ts_baa_aaa_diff.pdf",       10, 4.5),
  list(p5,  "fig_ts_reserves_tga.pdf",       10, 4.5),
  list(p6,  "fig_ts_ted_spread.pdf",         10, 4.5),
  list(p7,  "fig_ts_etf_prices.pdf",         10, 4.5),
  list(p8,  "fig_ts_on_rrp.pdf",             10, 4.5),
  list(p9,  "fig_ts_dxy.pdf",                10, 4.5),
  list(p10, "fig_ts_funding_vs_credit.pdf",  10, 4.5),
  list(p11, "fig_ts_cross_currency_basis.pdf", 10, 4.5)
)

for (item in plot_list) {
  ggsave(item[[2]], item[[1]], width = item[[3]], height = item[[4]], dpi = 300)
  cat(sprintf("  Saved: %s\n", item[[2]]))
}
cat("\n")

df_vol <- df_extended %>%
  select(Date, d_Aaa, d_Baa) %>% filter(!is.na(d_Aaa) & !is.na(d_Baa))

df_vol <- df_vol %>% mutate(
  vol_Aaa = zoo::rollapply(d_Aaa, width = 60, FUN = sd, fill = NA, align = "right"),
  vol_Baa = zoo::rollapply(d_Baa, width = 60, FUN = sd, fill = NA, align = "right"),
  vol_ratio = vol_Aaa / vol_Baa)

p_vol <- ggplot(df_vol %>% filter(!is.na(vol_ratio)), aes(x = Date, y = vol_ratio))
p_vol <- add_rec(p_vol)
p_vol <- p_vol +
  geom_hline(yintercept = 1, color = COL_GREY, linewidth = 0.4, linetype = "dashed") +
  geom_line(color = COL_GREEN, linewidth = 0.4, alpha = 0.85) +
  annotate("rect", xmin = as.Date("2019-09-15"), xmax = as.Date("2019-10-15"),
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#E57373") +
  annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2020-04-15"),
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#E57373") +
  labs(title = "Volatility Ratio: Aaa / Baa (60-Day Rolling)",
       subtitle = "Above 1 = Aaa more volatile (signals collateral-channel stress).",
       x = NULL, y = expression(sigma[Aaa] / sigma[Baa]),
       caption = "Red bands: Sept 2019 repo, March 2020 COVID.") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0.01, 0)) +
  coord_cartesian(ylim = c(0, max(df_vol$vol_ratio, na.rm = TRUE) * 1.05)) +
  theme_thesis()

ggsave("fig_ts_vol_ratio.pdf", p_vol, width = 10, height = 4.5, dpi = 300)
cat("  Saved: fig_ts_vol_ratio.pdf\n\n")



cat("================================================================\n")
cat("  CLEANING REPORT\n")
cat("================================================================\n\n")

cat(sprintf("  Source: Thesis/RAW DATA.xlsx (May 7 2026)\n"))
cat(sprintf("  Output: data_master.csv (%d rows, %d columns)\n", nrow(df), ncol(df)))
cat(sprintf("  Holiday zeros: %d | DXY zeros: %d\n", n_holidays, n_dxy_zero))
cat("  Dropped: SOFR, GC_Repo, EFFR (levels), TGCR, BGCR\n")
cat("  Forward-filled: Reserves, TGA, Fed_Assets, Liq_Swaps\n")
cat("  New variables: SHV, HYG_SHV, EMB, Baa_Aaa,\n")
cat("    d_Reserves, d_TGA, d_ON_RRP, d_Fed_Assets,\n")
cat("    d_Baa, d_Aaa, d_Baa_Aaa, lr_HYG, lr_LQD, lr_DXY,\n")
cat("    lr_TLT, lr_IEF, lr_SHV, lr_EMB,\n")
cat("    log_HYG_LQD, d_log_HYG_LQD, log_HYG_SHV, d_log_HYG_SHV\n\n")

cat("  KEY FINDING — SOFR proxy:\n")
cat("    The ON Treasury GC Repo Rate provides a validated proxy for SOFR\n")
cat("    before April 2018 (Fed FEDS Notes, 2019). This extends the\n")
cat("    SOFR-EFFR spread to May 2003, with a ~1 month gap in March 2018.\n")
cat("    The extended sample covers the GFC — the most important stress\n")
cat("    episode for the theoretical framework.\n\n")

cat("  SAMPLE RECOMMENDATIONS:\n")
cat("    1. Extended sample (2003--2026, 5726 SOFR-EFFR obs): Granger tests,\n")
cat("       descriptive analysis, full VAR with proxy.\n")
cat("    2. Post-2018 sample (2018--2026, 2006 obs): Primary VAR with actual SOFR.\n")
cat("    3. TED sample (2003--2022, 4684 obs): Robustness specification.\n")
cat("    4. Cross-currency basis (2024--2026, 523 obs): Descriptive only.\n\n")

cat("================================================================\n")
cat("  STEP 1 COMPLETE\n")
cat("================================================================\n")
