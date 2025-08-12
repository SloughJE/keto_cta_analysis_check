# Modeling Analysis

# Assumption checks

library(tidyverse)
library(performance)   # model_performance()
library(lmtest)        # bptest(), resettest()

library(tidyverse)
library(performance)
library(lmtest)

assumption_report <- function(formula, data) {
  m   <- lm(formula, data = data)
  n   <- nobs(m)
  p   <- length(coef(m))               # params incl. intercept
  perf<- performance::model_performance(m)
  
  # tests
  sh  <- shapiro.test(residuals(m))
  bp  <- bptest(m)
  rst <- resettest(m, power = 2:3, type = "fitted")
  
  # influence diagnostics
  cd  <- cooks.distance(m)
  h   <- hatvalues(m)
  rs  <- rstudent(m)
  cut_cd <- 4/n
  
  tibble(
    model        = deparse(formula),
    n            = n,
    slope        = unname(coef(m)[2]),
    slope_p      = coef(summary(m))[2, "Pr(>|t|)"],
    r2           = perf$R2,
    rmse         = perf$RMSE,
    # assumption tests
    bp_p         = bp$p.value,         # homoskedasticity
    shapiro_p    = sh$p.value,         # normality
    reset_p      = rst$p.value,        # linearity/specification
    # influence summary
    cooks_cut    = cut_cd,
    cooks_n      = sum(cd > cut_cd, na.rm = TRUE),
    cooks_prop   = mean(cd > cut_cd, na.rm = TRUE),
    cooks_max    = max(cd, na.rm = TRUE),
    cooks_ge1    = any(cd >= 1, na.rm = TRUE),
    leverage_hi2 = sum(h > 2*p/n, na.rm = TRUE),
    leverage_hi3 = sum(h > 3*p/n, na.rm = TRUE),
    rstud_gt3    = sum(abs(rs) > 3, na.rm = TRUE)
  )
}

forms <- list(
  `ΔNCPV ~ CAC_bl`  = delta_NCPV ~ V1_CAC,
  `ΔNCPV ~ NCPV_bl` = delta_NCPV ~ V1_Non_Calcified_Plaque_Volume,
  `ΔNCPV ~ PAV_bl`  = delta_NCPV ~ V1_PAV_pct,
  `ΔNCPV ~ TPS_bl`  = delta_NCPV ~ V1_Total_Plaque_Score
)

assumption_tbl <- purrr::map_dfr(forms, assumption_report, data = df)
assumption_tbl
View(assumption_tbl)


iwalk(forms, function(f, nm) {
  m <- lm(f, data = df)
  p <- check_model(m)
  print(p)  # force draw (esp. in scripts / knits)
})

library(patchwork)
#library(qqplotr) # if want loess band on Quantile Dev chart

purrr::iwalk(forms, function(f, nm) {
  m  <- lm(f, data = df)
  cm <- performance::check_model(m)
  pw <- getS3method("plot", "see_check_model")(cm) +   # invoke method explicitly
    plot_annotation(
      title = nm,
      #subtitle = paste(deparse(f), collapse = " "),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
      )
    )
  
  ggsave(
    file = file.path("figures/diagnostics", paste0("check_", gsub("[^A-Za-z0-9]+","_", nm), ".png")),
    plot = pw, device = ragg::agg_png, width = 8, height = 10, units = "in", dpi = 600, bg = "white"
  )
})
