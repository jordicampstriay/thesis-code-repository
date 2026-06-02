library(tidyverse)
library(lmtest)
library(sandwich)

select <- dplyr::select

WD <- "/Users/jordi/Downloads/University/TFG/Data/Thesis/Empirical Analysis"
setwd(WD)

d <- read.csv("../data_master.csv")
d$Date <- as.Date(d$Date)
d53 <- d[d$Date >= as.Date("2018-04-02") & !is.na(d$SOFR_EFFR), ]
d54 <- d[!is.na(d$DXY) & !is.na(d$HYG) & !is.na(d$LQD) & !is.na(d$SHV) & !is.na(d$EMB), ]

pdf("figures/fig_sofr_distribution.pdf", width = 8, height = 3.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.5, 0.8, 0))

sofr <- na.omit(d53$SOFR_EFFR)

hist(sofr, breaks = 100, col = "steelblue", border = "white",
     main = "(a) SOFR-EFFR Spread Distribution",
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
cat("Done: fig_sofr_distribution.pdf\n")

roll_granger <- function(data, x_var, y_var, lag_order = 5,
                         window = 400, step = 15) {
  dd <- data[complete.cases(data[, c(x_var, y_var)]), ]
  n  <- nrow(dd)
  dates <- c()
  pvals <- c()

  for (i in seq(window, n, by = step)) {
    idx <- (i - window + 1):i
    sub <- dd[idx, ]
    tryCatch({
      gt <- grangertest(sub[[y_var]] ~ sub[[x_var]], order = lag_order)
      pvals <- c(pvals, gt$`Pr(>F)`[2])
      dates <- c(dates, as.character(dd$Date[i]))
    }, error = function(e) {
      pvals <<- c(pvals, NA_real_)
      dates <<- c(dates, as.character(dd$Date[i]))
    })
  }
  data.frame(Date = as.Date(dates), p_value = pvals, stringsAsFactors = FALSE)
}

theme_thesis <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    legend.position = "bottom"
  )

cat("Computing rolling Granger: SOFR-EFFR -> credit...\n")

roll_baa <- roll_granger(d53, "SOFR_EFFR", "d_Baa")
roll_baa$Direction <- "SOFR-EFFR -> d_Baa (lag 5)"
cat("  SOFR->Baa:", sum(!is.na(roll_baa$p_value)), "valid points\n")

roll_aaa <- roll_granger(d53, "SOFR_EFFR", "d_Aaa")
roll_aaa$Direction <- "SOFR-EFFR -> d_Aaa (lag 5)"
cat("  SOFR->Aaa:", sum(!is.na(roll_aaa$p_value)), "valid points\n")

roll_a <- rbind(roll_baa, roll_aaa)

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
  labs(title = NULL, subtitle = NULL,
       x = NULL, y = "p-value", color = NULL) +
  theme_thesis
ggsave("figures/fig_a5_rolling_granger.pdf", p_a5, width = 10, height = 5)
cat("Done: fig_a5_rolling_granger.pdf\n")

cat("Computing rolling Granger: quantity channel...\n")

roll_res_aaa <- roll_granger(d53, "d_Reserves", "d_Aaa")
roll_res_aaa$Direction <- "Reserves -> Aaa (lag 5)"
cat("  Reserves->Aaa:", sum(!is.na(roll_res_aaa$p_value)), "valid points\n")

roll_tga_diff <- roll_granger(d53, "d_TGA", "d_Baa_Aaa")
roll_tga_diff$Direction <- "TGA -> Baa-Aaa (lag 5)"
cat("  TGA->Baa-Aaa:", sum(!is.na(roll_tga_diff$p_value)), "valid points\n")

roll_qty <- rbind(roll_res_aaa, roll_tga_diff)

p_qty <- ggplot(roll_qty, aes(x = Date, y = p_value, color = Direction)) +
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
  labs(title = NULL, subtitle = NULL,
       x = NULL, y = "p-value", color = NULL) +
  theme_thesis
ggsave("figures/fig_a6_rolling_granger_qty.pdf", p_qty, width = 10, height = 5)
cat("Done: fig_a6_rolling_granger_qty.pdf\n")

roll_cor_fn <- function(x, y, w) {
  n <- length(x)
  r <- rep(NA_real_, n)
  for (i in w:n) {
    idx <- (i - w + 1):i
    xi <- x[idx]; yi <- y[idx]
    ok <- complete.cases(xi, yi)
    if (sum(ok) >= 30) r[i] <- cor(xi[ok], yi[ok])
  }
  r
}

df_rc <- d54 %>%
  mutate(
    rc_HYG_LQD = roll_cor_fn(lr_DXY, d_log_HYG_LQD, 250),
    rc_HYG_SHV = roll_cor_fn(lr_DXY, d_log_HYG_SHV, 250),
    rc_EMB     = roll_cor_fn(lr_DXY, lr_EMB, 250)
  )

pdf("figures/fig_rolling_corr_dxy.pdf", width = 9, height = 5)
par(mar = c(4, 4, 1, 1), mgp = c(2.5, 0.8, 0))
plot(df_rc$Date, df_rc$rc_HYG_LQD, type = "n",
     ylim = c(-0.7, 0.3), xlab = "", ylab = "Correlation",
     main = "")
abline(h = 0, col = "grey60", lty = 2)
lines(df_rc$Date, df_rc$rc_HYG_LQD, col = "#2166AC", lwd = 1.5)
lines(df_rc$Date, df_rc$rc_HYG_SHV, col = "#D6604D", lwd = 1.5)
lines(df_rc$Date, df_rc$rc_EMB,     col = "#1B7837", lwd = 1.5)
abline(v = as.Date("2020-03-01"), col = "red", lty = 3)
text(as.Date("2020-03-01"), 0.25, "COVID", col = "red", cex = 0.8, pos = 4)
abline(v = as.Date("2008-09-15"), col = "darkred", lty = 3)
text(as.Date("2008-09-15"), 0.25, "Lehman", col = "darkred", cex = 0.8, pos = 4)
legend("bottomleft",
       legend = c("DXY vs HYG/LQD", "DXY vs HYG/SHV", "DXY vs EMB"),
       col = c("#2166AC", "#D6604D", "#1B7837"), lwd = 2, cex = 0.85,
       bg = "white")
dev.off()
cat("Done: fig_rolling_corr_dxy.pdf\n")

cat("Computing rolling Granger: DXY -> credit (long sample, ~3 min)...\n")

rg_hyglqd <- roll_granger(d54, "lr_DXY", "d_log_HYG_LQD")
rg_hyglqd$Direction <- "DXY -> HYG/LQD"
cat("  DXY->HYG/LQD:", sum(!is.na(rg_hyglqd$p_value)), "valid points\n")

rg_hygshv <- roll_granger(d54, "lr_DXY", "d_log_HYG_SHV")
rg_hygshv$Direction <- "DXY -> HYG/SHV"
cat("  DXY->HYG/SHV:", sum(!is.na(rg_hygshv$p_value)), "valid points\n")

rg_emb <- roll_granger(d54, "lr_DXY", "lr_EMB")
rg_emb$Direction <- "DXY -> EMB"
cat("  DXY->EMB:", sum(!is.na(rg_emb$p_value)), "valid points\n")

rg_all <- rbind(rg_hyglqd, rg_hygshv, rg_emb)

pdf("figures/fig_rolling_granger_dxy.pdf", width = 9, height = 5)
par(mar = c(4, 4, 1, 1), mgp = c(2.5, 0.8, 0))
plot(rg_hyglqd$Date, rg_hyglqd$p_value, type = "n",
     ylim = c(0, 1), xlab = "", ylab = "HAC-robust p-value",
     main = "")
abline(h = 0.05, col = "black", lty = 2)
lines(rg_hyglqd$Date, rg_hyglqd$p_value, col = "#2166AC", lwd = 1.3)
lines(rg_hygshv$Date, rg_hygshv$p_value, col = "#D6604D", lwd = 1.3)
lines(rg_emb$Date,    rg_emb$p_value,    col = "#1B7837", lwd = 1.3)
abline(v = as.Date("2020-03-01"), col = "red", lty = 3)
text(as.Date("2020-03-01"), 0.95, "COVID", col = "red", cex = 0.8, pos = 4)
abline(v = as.Date("2008-09-15"), col = "darkred", lty = 3)
text(as.Date("2008-09-15"), 0.95, "Lehman", col = "darkred", cex = 0.8, pos = 4)
legend("topright",
       legend = c("DXY -> HYG/LQD", "DXY -> HYG/SHV", "DXY -> EMB",
                  "5% significance"),
       col = c("#2166AC", "#D6604D", "#1B7837", "black"),
       lty = c(1, 1, 1, 2), lwd = c(2, 2, 2, 1), cex = 0.8, bg = "white")
dev.off()
cat("Done: fig_rolling_granger_dxy.pdf\n")

cat("\nAll figures fixed.\n")
