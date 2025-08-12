###

# Keto CTA Plaque data
library(tidyverse)

df <- read_csv("data/keto-cta-quant-and-semi-quant.csv")

glimpse(df)   # structure + column types
head(df)      # first rows

# Check for any NAs
sum(is.na(df))                     # total count
colSums(is.na(df))                  # count per column

# Summary statistics
summary(df)

# columns
colnames(df)

# create change scores for each Plaque metric
df <- df %>%
  mutate(
    delta_NCPV = V2_Non_Calcified_Plaque_Volume - V1_Non_Calcified_Plaque_Volume,
    delta_TPS  = V2_Total_Plaque_Score - V1_Total_Plaque_Score,
    delta_PAV  = V2_Percent_Atheroma_Volume - V1_Percent_Atheroma_Volume,
    delta_PAV_pct  = (V2_Percent_Atheroma_Volume - V1_Percent_Atheroma_Volume)*100,
    V1_PAV_pct = V1_Percent_Atheroma_Volume * 100
  )

# Helper functions
qstats2 <- function(x) round(quantile(x, probs = c(.25, .5, .75), na.rm = TRUE), 2)  # Q1, median, Q3
mean2   <- function(x) round(mean(x, na.rm = TRUE), 2)
sd2     <- function(x) round(sd(x, na.rm = TRUE), 2)
cmp_num <- function(a, b, tol) ifelse(is.na(a) | is.na(b), NA, abs(a - b) <= tol)

# ---- Observed stats ----
obs <- list(
  V1_NCPV = list(mean = mean2(df$V1_Non_Calcified_Plaque_Volume),
                 sd   = sd2(df$V1_Non_Calcified_Plaque_Volume),
                 q    = qstats2(df$V1_Non_Calcified_Plaque_Volume)),
  V1_CAC  = list(mean = mean2(df$V1_CAC),
                 sd   = sd2(df$V1_CAC),
                 q    = qstats2(df$V1_CAC)),
  V1_PAV  = list(mean = mean2(df$V1_PAV_pct),
                 sd   = sd2(df$V1_PAV_pct),
                 q    = qstats2(df$V1_PAV_pct)),
  V1_TPS  = list(mean = mean2(df$V1_Total_Plaque_Score),
                 sd   = sd2(df$V1_Total_Plaque_Score),
                 q    = qstats2(df$V1_Total_Plaque_Score)),
  delta_NCPV = list(q = qstats2(df$delta_NCPV)),
  delta_PAV  = list(q = qstats2(df$delta_PAV_pct))
)


# Expected Values (from text)
exp <- list(
  V1_NCPV = list(mean = 75.90, sd = 88.30, q = c(15.40, 44.00, 102.30)),
  V1_CAC  = list(mean = 50.30, sd =100.90, q = c( 0.00,  0.00,  54.00)),
  V1_PAV  = list(mean =  3.20, sd =  3.80, q = c( 0.50,  1.60,   4.90)),
  V1_TPS  = list(mean =  1.70, sd =  2.60, q = c( 0.00,  0.00,   2.00)),
  delta_NCPV = list(q = c( 9.30, 18.90, 47.00)),
  delta_PAV  = list(q = c( 0.30,  0.80,  1.70))
)


print_baseline <- function(name, obs_entry, exp_entry, tol = 0.20) {
  mean_ok <- cmp_num(obs_entry$mean, exp_entry$mean, tol)
  sd_ok   <- cmp_num(obs_entry$sd,   exp_entry$sd,   tol)
  q_ok    <- abs(obs_entry$q - exp_entry$q) <= tol
  
  msg <- if (all(mean_ok, sd_ok, q_ok)) "YES" else {
    diffs <- c()
    if (!isTRUE(mean_ok)) diffs <- c(diffs, "mean")
    if (!isTRUE(sd_ok))   diffs <- c(diffs, "SD")
    labs <- c("Q1","median","Q3")
    diffs <- c(diffs, labs[!q_ok])
    paste("NO ->", paste(diffs, collapse = ", "))
  }
  
  cat(sprintf(
    "\n%s\n  observed: mean ± SD = %.2f ± %.2f ; median (Q1–Q3) = %.2f (%.2f–%.2f)\n  expected: mean ± SD = %.2f ± %.2f ; median (Q1–Q3) = %.2f (%.2f–%.2f)\n  match? %s\n",
    name,
    obs_entry$mean, obs_entry$sd, obs_entry$q[2], obs_entry$q[1], obs_entry$q[3],
    exp_entry$mean, exp_entry$sd, exp_entry$q[2], exp_entry$q[1], exp_entry$q[3],
    msg
  ))
}

print_change <- function(name, obs_q, exp_q, tol = 0.2) {
  q_ok <- abs(obs_q - exp_q) <= tol
  msg <- if (all(q_ok)) "YES" else {
    labs <- c("Q1","median","Q3")
    paste("NO ->", paste(labs[!q_ok], collapse = ", "))
  }
  
  cat(sprintf(
    "\n%s\n  observed: median (Q1–Q3) = %.2f (%.2f–%.2f)\n  expected: median (Q1–Q3) = %.2f (%.2f–%.2f)\n  match? %s\n",
    name,
    obs_q[2], obs_q[1], obs_q[3],
    exp_q[2], exp_q[1], exp_q[3],
    msg
  ))
}

tol <- 0.25  # tolerance for comparisons after rounding

print_baseline("Baseline NCPV (mm^3)", obs$V1_NCPV, exp$V1_NCPV, tol)
print_baseline("Baseline CAC",         obs$V1_CAC,  exp$V1_CAC,  tol)
print_baseline("Baseline PAV (%)",     obs$V1_PAV,  exp$V1_PAV,  tol)
print_baseline("Baseline TPS",         obs$V1_TPS,  exp$V1_TPS,  tol)

print_change("Δ NCPV (mm^3)",          obs$delta_NCPV$q, exp$delta_NCPV$q, tol=0.5)
print_change("Δ PAV (%)",              obs$delta_PAV$q,  exp$delta_PAV$q,  tol)





#############################################

library(BayesFactor)
rscale_input = 0.2

bf <- regressionBF(delta_NCPV ~ V1_CAC, data = df, rscale = rscale_input, progress = FALSE)
BayesFactor::extractBF(bf)
# Posterior for parameters (gives means/SD/credible intervals)
post <- BayesFactor::posterior(bf[1], iterations = 50000, progress = FALSE)
summ <- t(apply(as.data.frame(post), 2, function(x)
  c(mean = mean(x), sd = sd(x),
    q2.5 = quantile(x, .025), q97.5 = quantile(x, .975))))
round(summ, 3)

bf <- regressionBF(delta_NCPV ~ V1_Non_Calcified_Plaque_Volume, data = df, rscale = rscale_input, progress = FALSE)
BayesFactor::extractBF(bf)
# Posterior for parameters (gives means/SD/credible intervals)
post <- BayesFactor::posterior(bf[1], iterations = 50000, progress = FALSE)
summ <- t(apply(as.data.frame(post), 2, function(x)
  c(mean = mean(x), sd = sd(x),
    q2.5 = quantile(x, .025), q97.5 = quantile(x, .975))))
round(summ, 3)

bf <- regressionBF(delta_NCPV ~ V1_Percent_Atheroma_Volume, data = df, rscale = rscale_input, progress = FALSE)
BayesFactor::extractBF(bf)
# Posterior for parameters (gives means/SD/credible intervals)
post <- BayesFactor::posterior(bf[1], iterations = 50000, progress = FALSE)
summ <- t(apply(as.data.frame(post), 2, function(x)
  c(mean = mean(x), sd = sd(x),
    q2.5 = quantile(x, .025), q97.5 = quantile(x, .975))))
round(summ, 3)

bf <- regressionBF(delta_NCPV ~ V1_Total_Plaque_Score, data = df, rscale = rscale_input, progress = FALSE)
BayesFactor::extractBF(bf)
# Posterior for parameters (gives means/SD/credible intervals)
post <- BayesFactor::posterior(bf[1], iterations = 50000, progress = FALSE)
summ <- t(apply(as.data.frame(post), 2, function(x)
  c(mean = mean(x), sd = sd(x),
    q2.5 = quantile(x, .025), q97.5 = quantile(x, .975))))
round(summ, 3)


# “Our analysis was carried out by two experts in data analysis, 
# and it was independently reviewed by an expert in statistics during the peer review process,”
# https://www.wired.com/story/how-one-trial-set-off-a-new-war-in-the-nutrition-world-keto-cholesterol-fat/

# https://citizensciencefoundation.org/keto-cta/

#Keto-CTA Open Access Supplemental
#The following table lists values from all 100 participants for both visits (V1 & V2) in two categories of analyses:
# Semi-Quantitative
# Total Plaque Score (TPS)
# Coronary Artery Calcification (CAC)
# Quantitative
# Non-Calcified Plaque Volume (NCPV)
# Total Calcified Plaque Volume (TCPV)
# Percent Atheroma Volume (PAV)
# Keto-CTA Table – Imaging Metrics

