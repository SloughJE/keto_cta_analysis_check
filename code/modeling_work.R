####
# lm work


df <- read_csv("data/keto-cta-quant-and-semi-quant.csv")

df$V1_Percent_Atheroma_Volume = df$V1_Percent_Atheroma_Volume*100
df$V2_Percent_Atheroma_Volume = df$V2_Percent_Atheroma_Volume*100


m <- lm(V2_Non_Calcified_Plaque_Volume ~ V1_Non_Calcified_Plaque_Volume,
        data = df)
summary(m)

library(lmtest)
library(car)

# 0) Quick visual checks
op <- par(mfrow = c(2, 2)); plot(m); par(op)  # residuals vs fitted, QQ, scale-location, influence

# 1) Normality of residuals
shapiro.test(residuals(m))                    # use with care for large n

# 2) Homoscedasticity
bptest(m)                                     # Breusch–Pagan
ncvTest(m)                                    # car’s non-constant variance test

# 3) Independence (autocorrelation in residuals)
dwtest(m)                                     # Durbin–Watson

# 4) Linearity (component+residual plot)
crPlots(m, terms = ~ V1_Non_Calcified_Plaque_Volume)

# 5) Influence / outliers
cooksd <- cooks.distance(m)
which(cooksd > 4/length(cooksd))              # indices of notable cases
outlierTest(m)                                 # Bonferroni test for studentized residuals
# influencePlot(m)                             # optional interactive plot


### model fails normality and homoscedasticity
library(lmtest); library(sandwich)

m_log <- lm(log1p(V2_Non_Calcified_Plaque_Volume) ~ log1p(V1_Non_Calcified_Plaque_Volume),
            data = df, na.action = na.exclude)


par(mfrow=c(2,2)); plot(m_log); par(mfrow=c(1,1))
bptest(m_log); shapiro.test(residuals(m_log))

df$delta_NCPV = df$V2_Non_Calcified_Plaque_Volume - df$V1_Non_Calcified_Plaque_Volume

m <- lm(V2_Non_Calcified_Plaque_Volume ~ V1_Non_Calcified_Plaque_Volume, data = df)
coeftest(m, vcov = vcovHC(m, type = "HC3"))
summary(m)
car::crPlots(m)

library(boot)
f <- function(d,i) coef(lm(V2_Non_Calcified_Plaque_Volume ~ V1_Non_Calcified_Plaque_Volume, data = d[i,]))[2]
b <- boot(df, f, R = 2000)
boot.ci(b, index = 1, type = c("perc","bca"))

m <- lm(V2_Non_Calcified_Plaque_Volume ~ V1_Non_Calcified_Plaque_Volume, data = df)
lmtest::coeftest(m, vcov = sandwich::vcovHC(m, type = "HC3"))

m_ols <- lm(V2_Non_Calcified_Plaque_Volume ~ V1_Non_Calcified_Plaque_Volume, data = df)
plot(fitted(m_ols), abs(residuals(m_ols)))  # pattern of spread vs mean
car::spreadLevelPlot(m_ols)         


###############################################
##############################################
library(nlme)

#
m0  <- gls(V2_Non_Calcified_Plaque_Volume ~ V1_Non_Calcified_Plaque_Volume, data=df, method="ML")  # equal variance

mP  <- gls(V2_Non_Calcified_Plaque_Volume ~ V1_Non_Calcified_Plaque_Volume, data=df,
           weights = varPower(form = ~ fitted(.)), method="ML")        # Var ∝ |μ|^{2δ}

mCP <- gls(V2_Non_Calcified_Plaque_Volume ~ V1_Non_Calcified_Plaque_Volume, data=df,
           weights = varConstPower(form = ~ fitted(.)), method="ML")   # c^2 + |μ|^{2δ}

mE  <- gls(V2_Non_Calcified_Plaque_Volume ~ V1_Non_Calcified_Plaque_Volume, data=df,
           weights = varExp(form = ~ fitted(.)), method="ML")          # Var ∝ exp(2δ μ)

anova(m0, mP, mCP, mE)   # LR tests (same fixed-effects)
AIC(m0, mP, mCP, mE)

best <- update(mP, method="REML")   #  mP won
summary(best)
intervals(best)
coef(best)                 # slope, intercept


# Diagnostics with *normalized* residuals
E <- resid(best, type = "normalized")
F <- fitted(best)

plot(F, E); abline(h = 0, lty = 2)   # should look flat
qqnorm(E); qqline(E)                  # heavy tails check
plot(F, abs(E))                       # should be patternless now

# Compare to OLS+HC3 to show conclusions are stable:
m_lm <- lm(V2_Non_Calcified_Plaque_Volume ~ V1_Non_Calcified_Plaque_Volume, data=df)
lmtest::coeftest(m_lm, vcov = sandwich::vcovHC(m_lm, type="HC3"))

delta0 <- as.numeric(coef(m_lm$modelStruct$varStruct, unconstrained = FALSE))

library(boot)
# 2) Pairs bootstrap
boot_fun <- function(d, i) {
  di <- d[i, ]
  fit <- gls(V2_Non_Calcified_Plaque_Volume ~ V1_Non_Calcified_Plaque_Volume,
             data = di,
             weights = varPower(value = delta0, form = ~ fitted(.)),  # start at delta0
             method = "ML",
             control = glsControl(msMaxIter = 200))
  unname(coef(fit)["V1_Non_Calcified_Plaque_Volume"])
}

set.seed(1)
b <- boot(df, boot_fun, R = 2000)
boot.ci(b, type = c("perc", "bca"))

check_model(best)

## prediction grid for the x-axis
new <- data.frame(
  V1_Non_Calcified_Plaque_Volume =
    seq(min(df$V1_Non_Calcified_Plaque_Volume),
        max(df$V1_Non_Calcified_Plaque_Volume), length.out = 200)
)

new$fit <- predict(best, newdata = new)

plot(V2_Non_Calcified_Plaque_Volume ~ V1_Non_Calcified_Plaque_Volume, data = df,
     pch = 21)
lines(new$V1_Non_Calcified_Plaque_Volume, new$fit, lwd = 2)

Xnew <- model.matrix(~ V1_Non_Calcified_Plaque_Volume, data = new)
bV   <- vcov(best)
se   <- sqrt(diag(Xnew %*% bV %*% t(Xnew)))
crit <- qt(0.975, df = nrow(df) - length(coef(best)))
new$lwr <- new$fit - crit*se
new$upr <- new$fit + crit*se


library(ggplot2)

# scatter of actuals
ggplot(df, aes(x = V1_Non_Calcified_Plaque_Volume,
               y = V2_Non_Calcified_Plaque_Volume)) +
  geom_point(alpha = .6) +
  # mean CI band (from vcov of the GLS)
  geom_ribbon(data = new,
              aes(x = V1_Non_Calcified_Plaque_Volume, ymin = lwr, ymax = upr),
              alpha = .15, inherit.aes = FALSE) +
  # fitted GLS mean line
  geom_line(data = new,
            aes(x = V1_Non_Calcified_Plaque_Volume, y = fit),
            linewidth = 1, inherit.aes = FALSE)


# with multiple covariates
#best2 <- gls( V2_Non_Calcified_Plaque_Volume ~ V1_Non_Calcified_Plaque_Volume + age + sex + ApoB,
#data = df,
#weights = varPower(~ fitted(.)),  # or try varPower(~ V1_Non_Calcified_Plaque_Volume)
#method = "REML"
#)
#summary(best2)

gls_fit_plot = ggplot(df, aes(x = V1_Non_Calcified_Plaque_Volume,
               y = V2_Non_Calcified_Plaque_Volume)) +
  geom_point(alpha = .6, size = 2, colour = "grey30") +
  geom_ribbon(
    data = new,
    aes(x = V1_Non_Calcified_Plaque_Volume, ymin = lwr, ymax = upr),
    fill = "grey50", alpha = .15, inherit.aes = FALSE
  ) +
  geom_line(
    data = new,
    aes(x = V1_Non_Calcified_Plaque_Volume, y = fit),
    colour = "red", linewidth = 1, inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45", linewidth = .4) +
  labs(
    title = "NCPV follow-up vs baseline",
    subtitle = expression(
      "GLS (REML) with power-of-the-mean variance: "
      ~ Var(epsilon[i]) == sigma^2 * "|" * mu[i] * "|"^{2*delta}
    ),
    x = expression(paste("NCPV baseline (mm"^3, ")")),
    y = expression(paste("NCPV follow-up (mm"^3, ")")),
    caption = "Shaded band: 95% pointwise CI for the GLS conditional mean E[Y|X] (from vcov())"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )



save_plot_simple <- function(p, filepath, width, height, dpi) {
  ggplot2::ggsave(
    filename = filepath,
    device = ragg::agg_png,
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

save_plot_simple(
  gls_fit_plot,
  "figures/plots_plaque_metrics/gls_fit__ncpv2_plot.png",
  6, 6, 600
)

###
# R2 or equiv

y    <- df$V2_Non_Calcified_Plaque_Volume
yhat <- fitted(best)

# weights proportional to 1/Var(epsilon_i)^{1/2}; square them for Σ^{-1} up to a constant
w    <- as.numeric(varWeights(best$modelStruct$varStruct))
w2   <- w^2

ybar_w <- weighted.mean(y, w2)
TSS_w  <- sum(w2 * (y - ybar_w)^2)
RSS_w  <- sum(w2 * (y - yhat  )^2)
R2_GLS <- 1 - RSS_w/TSS_w
R2_GLS


R2_corr <- cor(y, yhat)^2
R2_corr


####################################################






library(nlme)
m_gls <- gls(V2_Non_Calcified_Plaque_Volume ~ V1_Non_Calcified_Plaque_Volume,
             data = df, weights = varPower(form = ~ fitted(.)))
summary(m_gls); intervals(m_gls)  # look at delta; residual plots should de-fan

# Normalized residuals and fitted
r <- resid(m_gls, type = "normalized")
f <- fitted(m_gls)

# Base R: fitted vs residuals + QQ
op <- par(mfrow = c(1, 2), mar = c(4,4,1,1))
plot(f, r, xlab = "Fitted", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

qqnorm(r, main = "QQ-plot (normalized resid)")
qqline(r)
par(op)

linearHypothesis(m_gls,
                 c("(Intercept) = 0",
                   "V1_Non_Calcified_Plaque_Volume = 1"))

library(car)
linearHypothesis(m_gls, "(Intercept) = 0")
linearHypothesis(m_gls, "V1_Non_Calcified_Plaque_Volume = 1")

library(nlme)
inf_fun <- getS3method("influence", "gls")
inf <- inf_fun(m_gls, do.coef = TRUE)
head(inf[order(-inf[,"cook.d"]), c("cook.d","hat")], 5)

#################
##############

scatter_bf_with_stats <- function(metric) {
  g <- df_wide %>% dplyr::filter(metric == !!metric,
                                 is.finite(Baseline), is.finite(`Follow-up`))
  fit <- lm(`Follow-up` ~ Baseline, data = g)
  lims <- range(c(g$Baseline, g$`Follow-up`), na.rm = TRUE)
  
  lines_df <- tibble::tibble(
    type = c("Identity (no change)", "Linear fit"),
    slope = c(1, coef(fit)[["Baseline"]]),
    intercept = c(0, coef(fit)[["(Intercept)"]])
  )
  
  ggplot2::ggplot(g, aes(Baseline, `Follow-up`)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_abline(data = lines_df,
                         aes(slope = slope, intercept = intercept, linetype = type), linewidth = 0.9) +
    ggplot2::scale_linetype_manual(values = c("dashed", "solid"), name = NULL) +
    ggplot2::coord_equal(xlim = lims, ylim = lims, expand = TRUE) +
    ggplot2::labs(x = "Baseline", y = "Follow-up",
                  title = paste0(label_for(metric), " — Baseline vs Follow-up")) +
    theme_sl
}
scatter_bf_with_stats("Non_Calcified_Plaque_Volume")
residual_vs_baseline("Non_Calcified_Plaque_Volume")


scatter_bf_with_stats <- function(metric,
                                  transform = c("auto","none","log1p"),
                                  fit_color = "steelblue",
                                  id_color  = "black",
                                  id_linetype = "dashed",
                                  point_alpha = 0.6) {
  transform <- match.arg(transform)
  if (transform == "auto") transform <- if (use_log1p(metric)) "log1p" else "none"
  
  g <- df_wide %>%
    dplyr::filter(metric == !!metric, is.finite(Baseline), is.finite(`Follow-up`))
  if (!nrow(g)) return(ggplot2::ggplot() + ggplot2::theme_void())
  
  # ICC (two time points)
  d_icc <- g %>%
    dplyr::select(id, Baseline, `Follow-up`) %>%
    tidyr::pivot_longer(-id, names_to = "time", values_to = "y") %>%
    tidyr::drop_na(y)
  a <- stats::aov(y ~ factor(id), data = d_icc)
  s <- summary(a)[[1]]
  msb <- s["factor(id)", "Mean Sq"]; msw <- s["Residuals", "Mean Sq"]; k <- 2
  icc <- (msb - msw) / (msb + (k - 1) * msw)
  
  # Fit + predictions (shown back on original scale if log1p)
  grid <- data.frame(Baseline = seq(min(g$Baseline), max(g$Baseline), length.out = 200))
  if (transform == "log1p") {
    fit <- stats::lm(log1p(`Follow-up`) ~ log1p(Baseline), data = g)
    preds_g    <- exp(stats::predict(fit, newdata = g))    - 1
    preds_grid <- exp(stats::predict(fit, newdata = grid)) - 1
    r2 <- stats::cor(g$`Follow-up`, preds_g, use = "complete.obs")^2
    subtitle <- sprintf("R²=%.2f, ICC=%.2f (fit on log1p; shown in original units)", r2, icc)
  } else {
    fit <- stats::lm(`Follow-up` ~ Baseline, data = g)
    preds_grid <- stats::predict(fit, newdata = grid)
    r2 <- summary(fit)$r.squared
    subtitle <- sprintf("R²=%.2f, ICC=%.2f", r2, icc)
  }
  grid$`Follow-up` <- preds_grid
  
  # Lines (identity + fit) with legend & custom colors/linetypes
  id_df <- grid; id_df$`Follow-up` <- id_df$Baseline; id_df$type <- "Identity (no change)"
  fit_df <- grid; fit_df$type <- "Model fit"
  lines_df <- rbind(fit_df, id_df)
  
  lims <- range(c(g$Baseline, g$`Follow-up`, grid$`Follow-up`), na.rm = TRUE)
  
  ggplot2::ggplot(g, ggplot2::aes(Baseline, `Follow-up`)) +
    ggplot2::geom_point(alpha = point_alpha) +
    ggplot2::geom_line(data = lines_df,
                       ggplot2::aes(y = `Follow-up`, color = type, linetype = type),
                       linewidth = 0.9) +
    ggplot2::scale_color_manual(values = c("Model fit" = fit_color,
                                           "Identity (no change)" = id_color), name = NULL) +
    ggplot2::scale_linetype_manual(values = c("Model fit" = "solid",
                                              "Identity (no change)" = id_linetype), name = NULL) +
    ggplot2::coord_equal(xlim = lims, ylim = lims, expand = TRUE) +
    ggplot2::labs(
      title = paste0(label_for(metric), " — Baseline vs Follow-up"),
      subtitle = subtitle, x = "Baseline", y = "Follow-up"
    ) +
    theme_sl +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE),
           linetype = guide_legend(nrow = 1, byrow = TRUE))
}

p <- scatter_bf_with_stats("Non_Calcified_Plaque_Volume", transform = "log1p",
                           fit_color = "#2C7FB8", id_color = "gray30")
p


lm