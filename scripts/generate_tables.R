library(tseries)
library(xtable)
library(moments)

WD <- "/Users/jordi/Downloads/University/TFG/Data/Thesis/Empirical Analysis"
setwd(WD)

d <- read.csv("../data_master.csv")
d$Date <- as.Date(d$Date)

d53 <- d[d$Date >= as.Date("2018-04-02") & !is.na(d$SOFR_EFFR), ]
d54 <- d[!is.na(d$DXY) & !is.na(d$HYG) & !is.na(d$LQD) & !is.na(d$SHV) & !is.na(d$EMB), ]

run_adf <- function(v, label, data, transform, sample_label) {
  x <- na.omit(data[[v]])
  tt <- adf.test(x, alternative = "stationary")
  data.frame(
    Variable = label,
    Transform = transform,
    Sample = sample_label,
    N = length(x),
    ADF = round(tt$statistic, 2),
    p = ifelse(tt$p.value <= 0.01, "$< 0.01$", sprintf("%.3f", tt$p.value)),
    Stationary = ifelse(tt$p.value < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}

adf_rows <- list()

adf_rows[[1]]  <- run_adf("SOFR_EFFR", "SOFR--EFFR", d53, "Level", "Post-2018")
adf_rows[[2]]  <- run_adf("Baa_spread", "Baa spread", d53, "Level", "Post-2018")
adf_rows[[3]]  <- run_adf("Aaa_spread", "Aaa spread", d53, "Level", "Post-2018")
adf_rows[[4]]  <- run_adf("Baa_Aaa", "Baa--Aaa", d53, "Level", "Post-2018")
adf_rows[[5]]  <- run_adf("d_Baa", "$\\Delta$Baa", d53, "First diff.", "Post-2018")
adf_rows[[6]]  <- run_adf("d_Aaa", "$\\Delta$Aaa", d53, "First diff.", "Post-2018")
adf_rows[[7]]  <- run_adf("d_Baa_Aaa", "$\\Delta$(Baa--Aaa)", d53, "First diff.", "Post-2018")
adf_rows[[8]]  <- run_adf("d_Reserves", "$\\Delta$Reserves", d53, "First diff.", "Post-2018")
adf_rows[[9]]  <- run_adf("d_TGA", "$\\Delta$TGA", d53, "First diff.", "Post-2018")
adf_rows[[10]] <- run_adf("d_ON_RRP", "$\\Delta$ON RRP", d53, "First diff.", "Post-2018")

adf_rows[[11]] <- run_adf("lr_DXY", "$\\Delta\\log$ DXY", d54, "Log return", "2008--2026")
adf_rows[[12]] <- run_adf("d_log_HYG_LQD", "$\\Delta\\log$(HYG/LQD)", d54, "Log return", "2008--2026")
adf_rows[[13]] <- run_adf("d_log_HYG_SHV", "$\\Delta\\log$(HYG/SHV)", d54, "Log return", "2008--2026")
adf_rows[[14]] <- run_adf("lr_EMB", "$\\Delta\\log$ EMB", d54, "Log return", "2008--2026")

adf_df <- do.call(rbind, adf_rows)
cat("ADF results:\n")
print(adf_df)

xt <- xtable(adf_df,
  caption = "Augmented Dickey--Fuller Unit Root Tests",
  label = "tab:adf_all")
colnames(xt) <- c("Variable", "Transform.", "Sample", "$N$", "ADF stat", "$p$-value", "Stationary?")

print(xt, file = "tables/table_adf_all.tex",
      include.rownames = FALSE, sanitize.text.function = identity,
      booktabs = TRUE, floating = TRUE, table.placement = "H",
      scalebox = 0.80, comment = FALSE)

sources <- data.frame(
  Variable = c("SOFR", "EFFR", "Baa corporate yield", "Aaa corporate yield",
               "10-year Treasury yield", "Reserve Balances", "TGA Balance",
               "ON RRP Balance", "DXY (Dollar Index)",
               "HYG (High Yield ETF)", "LQD (Inv. Grade ETF)",
               "SHV (Short Treasury ETF)", "EMB (EM Bond ETF)"),
  Source = c("Federal Reserve Bank of New York",
             "Federal Reserve (H.15)",
             "FRED (Moody's)",
             "FRED (Moody's)",
             "FRED",
             "Federal Reserve (H.4.1)",
             "Federal Reserve (H.4.1)",
             "Federal Reserve Bank of New York",
             "Google Finance",
             "Google Finance",
             "Google Finance",
             "Google Finance",
             "Google Finance"),
  Frequency = c("Daily", "Daily", "Daily", "Daily", "Daily",
                "Weekly (fwd-filled)", "Weekly (fwd-filled)",
                "Daily", "Daily", "Daily", "Daily", "Daily", "Daily"),
  Start = c("Apr 2018", "Apr 2018", "May 2003", "May 2003", "May 2003",
            "May 2003", "May 2003", "Mar 2014", "Jan 2006",
            "Apr 2007", "Jan 2005", "Jan 2007", "Dec 2007"),
  stringsAsFactors = FALSE
)

xt2 <- xtable(sources,
  caption = "Data Sources and Availability",
  label = "tab:data_sources")

print(xt2, file = "tables/table_data_sources.tex",
      include.rownames = FALSE, sanitize.text.function = identity,
      booktabs = TRUE, floating = TRUE, table.placement = "H",
      scalebox = 0.82, comment = FALSE)

cor_vars <- c("SOFR_EFFR", "d_Baa", "d_Aaa", "d_Baa_Aaa",
              "d_Reserves", "d_TGA")
cor_labels <- c("SOFR--EFFR", "$\\Delta$Baa", "$\\Delta$Aaa",
                "$\\Delta$(Baa--Aaa)", "$\\Delta$Reserves", "$\\Delta$TGA")

cor_data <- d53[, cor_vars]
cor_data <- cor_data[complete.cases(cor_data), ]
cor_mat <- cor(cor_data)
rownames(cor_mat) <- colnames(cor_mat) <- cor_labels

xt3 <- xtable(cor_mat,
  caption = "Contemporaneous Correlation Matrix (Post-2018, Stationary Variables)",
  label = "tab:corr_post2018",
  digits = 3)

print(xt3, file = "tables/table_corr_post2018.tex",
      sanitize.text.function = identity,
      booktabs = TRUE, floating = TRUE, table.placement = "H",
      scalebox = 0.78, comment = FALSE)

pdf("figures/fig_sofr_distribution.pdf", width = 8, height = 3.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.5, 0.8, 0))

sofr <- na.omit(d53$SOFR_EFFR)

hist(sofr, breaks = 100, col = "steelblue", border = "white",
     main = "(a) SOFR--EFFR Spread Distribution",
     xlab = "Spread (percentage points)", ylab = "Frequency",
     xlim = c(-0.15, 0.5))
abline(v = 0, col = "red", lty = 2, lwd = 2)

plot(d53$Date, d53$SOFR_EFFR, type = "l", col = "steelblue", lwd = 0.8,
     main = "(b) SOFR-EFFR Spread (2018-2026)",
     xlab = "", ylab = "Spread (pp)")
abline(h = 0, col = "red", lty = 2)
abline(v = as.Date("2020-03-01"), col = "darkred", lty = 3)
text(as.Date("2020-03-01"), max(d53$SOFR_EFFR, na.rm=TRUE)*0.8, "COVID",
     col = "darkred", cex = 0.7, pos = 4)
dev.off()

ss_vars <- c("SOFR_EFFR", "d_Baa", "d_Aaa", "d_Baa_Aaa",
             "d_Reserves", "d_TGA", "d_ON_RRP")
ss_labels <- c("SOFR--EFFR", "$\\Delta$Baa", "$\\Delta$Aaa",
               "$\\Delta$(Baa--Aaa)", "$\\Delta$Reserves",
               "$\\Delta$TGA", "$\\Delta$ON RRP")

ss_rows <- lapply(seq_along(ss_vars), function(i) {
  x <- na.omit(d53[[ss_vars[i]]])
  data.frame(
    Variable = ss_labels[i],
    N = length(x),
    Mean = round(mean(x), 5),
    SD = round(sd(x), 5),
    Skew = round(skewness(x), 2),
    Kurt = round(kurtosis(x), 2),
    Min = round(min(x), 4),
    Max = round(max(x), 4),
    stringsAsFactors = FALSE
  )
})
ss_df <- do.call(rbind, ss_rows)

xt4 <- xtable(ss_df,
  caption = "Summary Statistics: Post-2018 Sample (Stationary Variables)",
  label = "tab:sumstats_post2018",
  digits = c(0, 0, 0, 5, 5, 2, 2, 4, 4))

print(xt4, file = "tables/table_sumstats_post2018.tex",
      include.rownames = FALSE, sanitize.text.function = identity,
      booktabs = TRUE, floating = TRUE, table.placement = "H",
      scalebox = 0.82, comment = FALSE)

cat("\nAll tables and figures generated.\n")
