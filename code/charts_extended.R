#####################
#### CHARTS #########
#####################
# Keto CTA Plaque data
library(tidyverse)
library(stringr)
library(purrr)
library(gt)
library(ggplot2)
library(ggtext)

df <- read_csv("data/keto-cta-quant-and-semi-quant.csv")

df$V1_Percent_Atheroma_Volume = df$V1_Percent_Atheroma_Volume*100
df$V2_Percent_Atheroma_Volume = df$V2_Percent_Atheroma_Volume*100

metric_keys <- c(
  "Non_Calcified_Plaque_Volume",
  "Total_Calcified_Plaque_Volume",
  "Percent_Atheroma_Volume",
  "Total_Plaque_Score",
  "CAC"
)

metric_labels <- c(
  Non_Calcified_Plaque_Volume   = "Noncalcified plaque volume (mm³)",
  Total_Calcified_Plaque_Volume = "Calcified plaque volume (mm³)",
  Percent_Atheroma_Volume       = "Percent atheroma volume (%)",
  Total_Plaque_Score            = "Total plaque score",
  CAC                           = "Coronary artery calcium score"
)

metric_labels_short <- c(
  Non_Calcified_Plaque_Volume   = "NCPV (mm³)",
  Total_Calcified_Plaque_Volume = "TCPV (mm³)",
  Percent_Atheroma_Volume       = "PAV (%)",
  Total_Plaque_Score            = "TPS",
  CAC                           = "CAC"
)

# Prefer short label (if available), else fall back to long label
label_for <- function(metric) {
  lab <- if (exists("metric_labels_short")) metric_labels_short[[metric]] else NULL
  if (is.null(lab) || is.na(lab)) lab <- metric_labels[[metric]]
  unname(lab)
}


# tidy long format for V1/V2
df_long <- df %>%
  mutate(id = row_number()) %>%
  pivot_longer(
    cols = matches("^V[12]_"),
    names_to = c("visit","metric"),
    names_pattern = "(V[12])_(.*)",
    values_to = "value"
  ) %>%
  filter(metric %in% metric_keys) %>%
  mutate(
    visit = if_else(visit == "V1", "Baseline", "Follow-up"),
    metric_lab = metric_labels[metric]
  )

# wide for paired calcs per metric
df_wide <- df_long %>%
  select(id, metric, visit, value) %>%
  pivot_wider(names_from = visit, values_from = value) %>%
  mutate(delta = `Follow-up` - Baseline,
         mean12 = (`Follow-up` + Baseline)/2)

# theme (slide-friendly, unobtrusive)
theme_sl <- theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

#  Waterfall 
waterfall <- function(metric) {
  g <- df_wide %>% dplyr::filter(metric == !!metric) %>%
    dplyr::arrange(delta) %>%
    dplyr::mutate(rank = dplyr::row_number())
  ggplot(g, aes(x = rank, y = delta)) +
    geom_col() +
    geom_hline(yintercept = 0, linewidth = 0.8) +
    labs(x = "Participant (ordered by Δ)", y = "Δ (FU − BL)",
         title = paste0(label_for(metric), " — Waterfall of Δ")) +
    theme_sl
}

# If deltas are integer-ish (e.g., CAC), use 1-unit bins centered on integers.
delta_hist <- function(metric, symmetric_limits = FALSE) {
  g <- df_wide %>% dplyr::filter(metric == !!metric)
  x <- g$delta
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0) {
    return(ggplot() + theme_void() +
             ggtitle(paste0(label_for(metric), " — no data")))
  }
  
  # helper
  is_integerish <- function(v) all(abs(v - round(v)) < 1e-8)
  
  med <- median(x, na.rm = TRUE)
  q   <- as.numeric(quantile(x, c(.25, .75), na.rm = TRUE))
  
  p <- ggplot(g, aes(x = delta)) 
  
  if (is_integerish(x)) {
    p <- p + geom_histogram(binwidth = 1, boundary = 0.5, closed = "left",
                            color = "white")
  } else {
    rng <- range(x, na.rm = TRUE)
    iqr <- IQR(x, na.rm = TRUE)
    bw  <- 2 * iqr / (n^(1/3))
    if (!is.finite(bw) || bw <= 0) bw <- diff(rng) / max(10, ceiling(n^(1/3)))
    p <- p + geom_histogram(binwidth = bw, boundary = 0, closed = "left",
                            color = "white")
  }
  
  p <- p +
    geom_vline(xintercept = 0, linewidth = 0.8) +
    geom_vline(xintercept = med, linetype = "dashed") +
    labs(x = "Change (Follow-up − Baseline)", y = "Count",
         title = paste0(label_for(metric), " — Δ distribution")) +
    theme_sl +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02)))
  
  if (symmetric_limits) {
    m <- max(abs(range(x, na.rm = TRUE)), na.rm = TRUE)
    p <- p + coord_cartesian(xlim = c(-m, m))
  }
  
  p
}

# Baseline vs Follow-up scatter — identity only
scatter_bf <- function(metric) {
  g <- df_wide %>%
    dplyr::filter(metric == !!metric, is.finite(Baseline), is.finite(`Follow-up`))
  
  lims <- range(c(g$Baseline, g$`Follow-up`), na.rm = TRUE)
  
  p <- ggplot(g, aes(x = Baseline, y = `Follow-up`)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0) +
    coord_equal(xlim = lims, ylim = lims, expand = TRUE) +
    labs(
      x = "Baseline",
      y = "Follow-up",
      title = paste0(label_for(metric), " — Baseline vs Follow-up"),
      subtitle = if (use_log1p(metric)) "Axes on log1p scale" else NULL
    ) +
    theme_sl
  
  if (use_log1p(metric)) {
    nice_breaks <- function(l) {
      v <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
      v[v >= l[1] & v <= l[2]]
    }
    p <- p +
      scale_x_continuous(trans = "log1p", breaks = nice_breaks(lims)) +
      scale_y_continuous(trans = "log1p", breaks = nice_breaks(lims))
  }
  
  p
}

prop_positive <- function(metric, threshold = 0, y_max = 1, show_labels = FALSE) {
  g <- df_long %>%
    dplyr::filter(metric == !!metric) %>%
    dplyr::mutate(
      visit = factor(visit, levels = c("Baseline", "Follow-up")),
      pos   = value > threshold
    ) %>%
    dplyr::group_by(visit) %>%
    dplyr::summarise(prop = mean(pos, na.rm = TRUE), n = dplyr::n(), .groups = "drop")
  
  p <- ggplot(g, aes(x = visit, y = prop)) +
    geom_col(width = 0.7) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = seq(0, y_max, by = 0.2),
      labels = scales::percent_format()
    ) +
    labs(
      x = NULL, y = "Proportion",
      title = paste0(label_for(metric), " — Proportion > ", threshold)
    ) +
    theme_sl +
    scale_x_discrete(drop = FALSE)
  
  if (isTRUE(show_labels)) {
    p <- p +
      geom_text(aes(label = scales::percent(prop, accuracy = 1)),
                vjust = -0.3, size = 4) +
      coord_cartesian(ylim = c(0, y_max))
  }
  
  p
}

# Bland–Altman with robust/transform options + numeric annotations
bland_altman <- function(metric,
                         transform = if (use_log1p(metric)) "log1p" else "none",
                         robust = TRUE,
                         annotate = TRUE,
                         digits = NULL) {
  transform <- match.arg(transform, c("none","log1p"))
  g <- df_wide %>% dplyr::filter(metric == !!metric)
  
  if (transform == "log1p") {
    x <- log1p(g$`Follow-up`); y <- log1p(g$Baseline)
    df <- tibble::tibble(x_mean = (x + y) / 2, y_diff = x - y)
    x_lab <- "Mean of log1p(Baseline, Follow-up)"
    y_lab <- "Δ log1p(Follow-up) − log1p(Baseline)"
    sub   <- "Log1p Bland–Altman (handles zeros; stabilizes spread)"
    if (is.null(digits)) digits <- 2
  } else {
    df <- tibble::tibble(x_mean = g$mean12, y_diff = g$delta)
    x_lab <- "Mean of Baseline & Follow-up"
    y_lab <- "Δ (Follow-up − Baseline)"
    sub   <- NULL
    if (is.null(digits)) digits <- 1
  }
  
  df <- df %>% dplyr::filter(is.finite(x_mean), is.finite(y_diff))
  if (nrow(df) == 0) return(ggplot() + theme_void() +
                              ggtitle(paste0(label_for(metric), " — Bland–Altman (no data)")))
  
  if (robust) {
    center <- stats::median(df$y_diff, na.rm = TRUE)
    spread <- 1.4826 * stats::mad(df$y_diff, center = center, constant = 1, na.rm = TRUE)
  } else {
    center <- mean(df$y_diff, na.rm = TRUE)
    spread <- sd(df$y_diff, na.rm = TRUE)
  }
  loa <- center + c(-1, 1) * 1.96 * spread
  
  ggplot(df, aes(x = x_mean, y = y_diff)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = center) +
    geom_hline(yintercept = loa, linetype = "dashed") +
    labs(x = x_lab, y = y_lab,
         title = paste0(label_for(metric), " — Bland–Altman"),
         subtitle = sub) +
    theme_sl +
    if (isTRUE(annotate)) {
      xr <- range(df$x_mean, na.rm = TRUE); yr <- range(df$y_diff, na.rm = TRUE)
      xpos <- xr[1] + 0.02 * diff(xr); ypos <- yr[2] - 0.03 * diff(yr)
      outside <- sum(df$y_diff < loa[1] | df$y_diff > loa[2], na.rm = TRUE)
      fmt <- function(z) formatC(z, format = "f", digits = digits)
      annotate("text", x = xpos, y = ypos, hjust = 0, vjust = 1,
               label = paste0(
                 if (robust) "Center (median Δ): " else "Center (mean Δ): ", fmt(center), "\n",
                 "LoA: ", fmt(loa[1]), " to ", fmt(loa[2]), "\n",
                 "Outside LoA: ", outside
               ), size = 4)
    }
}


#use_log1p <- function(metric) metric %in% c("CAC")

# Examples:
# bland_altman("Non_Calcified_Plaque_Volume", transform = "log1p", robust = TRUE)
# bland_altman("Non_Calcified_Plaque_Volume", transform = "none",  robust = TRUE)


# Examples to render:
# 
# delta_hist("Non_Calcified_Plaque_Volume")
# delta_hist("CAC") # need to deal with separately
# delta_hist("Total_Calcified_Plaque_Volume")
# 
# waterfall("Non_Calcified_Plaque_Volume")
# waterfall("Total_Calcified_Plaque_Volume")
# 
# scatter_bf("Non_Calcified_Plaque_Volume")
# 
# 
# prop_positive("CAC", threshold = 0, show_labels = TRUE)
# prop_positive("Non_Calcified_Plaque_Volume", threshold = 0, show_labels = TRUE)


########## 


plot_two_timepoints <- function(df,
                                var1, var2,
                                y_lab = "Value",
                                title = "",
                                label1 = "Baseline",
                                label2 = "1 Year",
                                scale_fn = identity) { 
  
  scale_fn <- rlang::as_function(scale_fn)  # allow ~ .x * 100 or function(x) x*100
  
  # Long data
  dat <- df %>%
    transmute(
      id = row_number(),
      !!label1 := scale_fn({{ var1 }}),
      !!label2 := scale_fn({{ var2 }})
    ) %>%
    drop_na() %>%
    pivot_longer(dplyr::all_of(c(label1, label2)),
                 names_to = "Time", values_to = "value") %>%
    mutate(x = as.numeric(factor(Time, levels = c(label1, label2))))
  
  # Summaries
  s <- dat %>%
    group_by(x) %>%
    summarise(
      q1  = quantile(value, 0.25),
      med = median(value),
      q3  = quantile(value, 0.75),
      .groups = "drop"
    ) %>%
    arrange(x)
  
  # IQR band polygon
  iqr_band <- tibble(
    x = c(1, 2, 2, 1),
    y = c(s$q1[1], s$q1[2], s$q3[2], s$q3[1])
  )
  
  print(iqr_band)
  # Median segment
  segment_df <- tibble(
    x = 1, xend = 2,
    y = s$med[s$x == 1],
    yend = s$med[s$x == 2]
  )
  
  y_top <- max(10, ceiling(max(dat$value, na.rm = TRUE)))
  minor_y <- seq(0, y_top, by = 1)
  
  ggplot(dat, aes(x = x, y = value, group = id)) +
    geom_polygon(data = iqr_band, aes(x = x, y = y),
                 inherit.aes = FALSE, alpha = 0.2, fill = "red") +
    geom_line(alpha = 0.25, colour = "grey50", linewidth = 0.4) +
    geom_point(size = 1.2, alpha = 0.30) +
    geom_segment(data = segment_df,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, colour = "red", linewidth = 1.2) +
    scale_x_continuous(breaks = c(1, 2), labels = c(label1, label2),
                       expand = expansion(mult = 0.45)   # try 0.3–0.5
    ) +
    scale_y_continuous(limits = c(0, y_top), minor_breaks = minor_y) +
    labs(x = NULL, y = y_lab, title = title) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_markdown(hjust = 0.5),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(face = "bold"),
      axis.text.x  = element_text(face = "bold"),
      axis.title.y = element_markdown(face = "bold"),
      panel.border = element_blank(),
      panel.grid.major = element_line(color = "black", linewidth = 0.2),
      panel.grid.minor.x = element_blank()
    )
}




figure_ncpv = plot_two_timepoints(df,
                                V1_Non_Calcified_Plaque_Volume, V2_Non_Calcified_Plaque_Volume,
                                y_lab = "NCPV mm<sup>3</sup>", title = "ΔNCPV (mm<sup>3</sup>)",
                                label1 = "Baseline", label2 = "1 Year"
)

figure_tps = plot_two_timepoints(df,
                                 V1_Total_Plaque_Score, V2_Total_Plaque_Score,
                                 y_lab = "TPS", title = "ΔTPS ",
                                 label1 = "Baseline", label2 = "1 Year"
)


figure_cac = plot_two_timepoints(df,
                                 V1_CAC, V2_CAC,
                                 y_lab = "CAC", title = "ΔCAC ",
                                 label1 = "Baseline", label2 = "1 Year"
)

figure_pav = plot_two_timepoints(df,
                                 V1_Percent_Atheroma_Volume, V2_Percent_Atheroma_Volume,
                                 y_lab = "PAV %", title = "ΔPAV (percentage point change) ",
                                 label1 = "Baseline", label2 = "1 Year"
)

figure_tcpv = plot_two_timepoints(df,
                                 V1_Total_Calcified_Plaque_Volume, V2_Total_Calcified_Plaque_Volume,
                                 y_lab = "TCPV", title = "ΔTCPV",
                                 label1 = "Baseline", label2 = "1 Year"
)

figure_ncpv = plot_two_timepoints(df,
                                V1_Non_Calcified_Plaque_Volume, V2_Non_Calcified_Plaque_Volume,
                                y_lab = "NCPV mm<sup>3</sup>", title = "ΔNCPV (mm<sup>3</sup>)",
                                label1 = "Baseline", label2 = "1 Year"
)

figure_tps = plot_two_timepoints(df,
                                 V1_Total_Plaque_Score, V2_Total_Plaque_Score,
                                 y_lab = "TPS", title = "ΔTPS ",
                                 label1 = "Baseline", label2 = "1 Year"
)


figure_cac = plot_two_timepoints(df,
                                 V1_CAC, V2_CAC,
                                 y_lab = "CAC", title = "ΔCAC ",
                                 label1 = "Baseline", label2 = "1 Year"
)

figure_pav = plot_two_timepoints(df,
                                 V1_Percent_Atheroma_Volume, V2_Percent_Atheroma_Volume,
                                 y_lab = "PAV %", title = "ΔPAV (percentage point change) ",
                                 label1 = "Baseline", label2 = "1 Year"
)

figure_tcpv = plot_two_timepoints(df,
                                 V1_Total_Calcified_Plaque_Volume, V2_Total_Calcified_Plaque_Volume,
                                 y_lab = "TCPV", title = "ΔTCPV",
                                 label1 = "Baseline", label2 = "1 Year"
)

slope_plots <- list(
  Non_Calcified_Plaque_Volume   = figure_ncpv,
  Total_Calcified_Plaque_Volume = figure_tcpv,
  Percent_Atheroma_Volume       = figure_pav,
  Total_Plaque_Score            = figure_tps,
  CAC                           = figure_cac
)


use_log1p <- function(metric) metric %in% c("Non_Calcified_Plaque_Volume","Total_Calcified_Plaque_Volume","CAC")


library(ragg)    # for ragg::agg_png

# Where to save and in what formats
out_dir  <- "figures/plots_plaque_metrics"
formats  <- c("png")     # add "pdf" if you want vector too
dpi_default <- 600       # slide-friendly
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Short stub for filenames (prefers short labels if available)
get_stub <- function(metric) {
  # prefer short label if present; else the raw key
  lab <- if (exists("metric_labels_short") && !is.null(metric_labels_short[[metric]])) {
    unname(metric_labels_short[[metric]])
  } else {
    metric
  }
  
  # normalize to ASCII and make filename-safe
  lab <- gsub("³", "3", lab, fixed = TRUE)     # mm³ -> mm3
  lab <- iconv(lab, to = "ASCII//TRANSLIT")    # general accent removal
  lab <- gsub("%", "pct", lab, fixed = TRUE)   # % -> pct
  lab <- gsub("[^A-Za-z0-9]+", "_", lab)       # non-alnum -> _
  lab <- gsub("^_+|_+$", "", lab)              # trim underscores
  tolower(lab)
}

# Device-aware saver (PNG via ragg; optional PDF via cairo)
save_plot <- function(plot, file_stub, ext = "png",
                      width = 8, height = 5, dpi = dpi_default, bg = "white") {
  filepath <- file.path(out_dir, paste0(file_stub, ".", ext))
  if (tolower(ext) == "png") {
    ggsave(filename = filepath, plot = plot, device = ragg::agg_png,
           width = width, height = height, units = "in",
           dpi = dpi, bg = bg)
  } else if (tolower(ext) == "pdf") {
    ggsave(filename = filepath, plot = plot, device = grDevices::cairo_pdf,
           width = width, height = height, units = "in")
  } else {
    ggsave(filename = filepath, plot = plot,
           width = width, height = height, units = "in",
           dpi = dpi, bg = bg)
  }
  invisible(filepath)
}

# Build the five plots for one metric
plots_for_metric <- function(metric) {
  list(
    slope     = slope_plots[[metric]],          # your prebuilt figure_* object
    delta     = delta_hist(metric),
    waterfall = waterfall(metric),
    scatter   = scatter_bf(metric),
    prop      = prop_positive(metric)
    #bland_altman = bland_altman(metric)
  )
}

# Save all five for a given metric, with sensible sizes for slides
save_metric_plots <- function(metric) {
  stub <- get_stub(metric)
  p    <- plots_for_metric(metric)
  if (!inherits(p$slope, "ggplot")) stop("No slope plot found for metric: ", metric)
  
  sizes <- list(                         # width, height in inches
    slope     = c(6, 6),                # wide banner
    delta     = c(6, 4),
    waterfall = c(6, 4),
    scatter   = c(6, 6),
    prop      = c(6, 6)
    #bland_altman = c(6,6)
  )
  
  walk(names(p), function(name) {
    plt <- p[[name]]
    if (inherits(plt, "ggplot")) {
      wh <- sizes[[name]]
      walk(formats, ~ save_plot(
        plot = plt,
        file_stub = paste0(stub, "_", name),
        ext = .x,
        width = wh[1], height = wh[2],
        bg = "white"     # or "transparent"
      ))
    }
  })
}

# Run for all metrics
walk(metric_keys, save_metric_plots)


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

#scatter_bf("Non_Calcified_Plaque_Volume")

save_plot_simple(
  scatter_bf("Non_Calcified_Plaque_Volume"),
  "figures/plots_plaque_metrics/ncpv_scatter_non_log.png",
  6, 6, 600
)

#######

ecdf_shift <- function(metric) {
  d <- df_long |> dplyr::filter(metric == !!metric, is.finite(value))
  ggplot(d, aes(x = value, color = visit)) +
    stat_ecdf(geom = "step", linewidth = 0.8) +
    labs(x = label_for(metric), y = "ECDF",
         title = paste0(label_for(metric), " — distribution shift")) +
    theme_sl
}

ecdf_shift("Non_Calcified_Plaque_Volume")



# helper: short label -> remove "(...)" and spaces
clean_metric_lab <- function(metric) {
  lab <- label_for(metric)                   # uses short if available
  lab <- gsub("\\s*\\([^)]*\\)", "", lab)    # drop everything in (...)
  lab <- gsub("\\s+", "", lab)               # drop spaces
  lab
}

delta_corr_heatmap <- function(method = "spearman") {
  M <- df_wide |>
    dplyr::select(id, metric, delta) |>
    tidyr::pivot_wider(id_cols = id, names_from = metric, values_from = delta) |>
    dplyr::select(-id)
  
  C <- stats::cor(as.matrix(M), use = "pairwise.complete.obs", method = method)
  
  levs <- intersect(metric_keys, colnames(C))
  lab_map <- setNames(vapply(levs, clean_metric_lab, character(1)), levs)
  
  tibble::as_tibble(C, rownames = "m1") |>
    tidyr::pivot_longer(-m1, names_to = "m2", values_to = "r") |>
    dplyr::mutate(m1 = factor(m1, levels = levs),
                  m2 = factor(m2, levels = levs)) |>
    ggplot2::ggplot(ggplot2::aes(m1, m2, fill = r)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = scales::number(r, accuracy = 0.01))) +
    ggplot2::scale_fill_gradient2(limits = c(-1, 1)) +
    ggplot2::scale_x_discrete(labels = lab_map) +
    ggplot2::scale_y_discrete(labels = lab_map) +
    ggplot2::labs(x = NULL, y = NULL, title = "Correlation of Δ across metrics") +
    theme_sl
}

delta_corr_heatmap()

save_plot_simple(
  delta_corr_heatmap(),
  "figures/plots_plaque_metrics/corr_deltas.png",
  6, 6, 600
)


icc_two_timepoints <- function(metric) {
  d <- df_wide %>%
    dplyr::filter(metric == !!metric) %>%
    dplyr::select(id, Baseline, `Follow-up`) %>%
    tidyr::pivot_longer(-id, names_to = "time", values_to = "y") %>%
    tidyr::drop_na(y)
  
  if (nrow(d) == 0) return(data.frame(metric = metric, ICC = NA_real_))
  
  a <- stats::aov(y ~ factor(id), data = d)
  s <- summary(a)[[1]]
  ms_between <- s["factor(id)", "Mean Sq"]
  ms_within  <- s["Residuals",  "Mean Sq"]
  k <- 2
  icc <- (ms_between - ms_within) / (ms_between + (k - 1) * ms_within)
  
  data.frame(
    metric = metric,
    ICC = icc,
    MS_between = ms_between,
    MS_within = ms_within,
    k = k,
    n_subjects = dplyr::n_distinct(d$id)
  )
}

icc_two_timepoints("Non_Calcified_Plaque_Volume")

baseline_quartile_delta <- function(metric, n_bins = 4) {
  g <- df_wide %>%
    dplyr::filter(metric == !!metric, is.finite(Baseline), is.finite(delta)) %>%
    dplyr::mutate(
      bin = cut(
        Baseline,
        breaks = quantile(Baseline, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE),
        include.lowest = TRUE
      )
    ) %>%
    dplyr::group_by(bin) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::ungroup()
  
  ggplot2::ggplot(g, ggplot2::aes(x = bin, y = delta)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.8) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.5, size = 1.4) +
    ggplot2::labs(
      x = "Baseline quantile",
      y = "Δ (Follow-up − Baseline)",
      title = paste0(label_for(metric), " — Δ by baseline quantile")
    ) +
    theme_sl
}

baseline_quartile_delta("Non_Calcified_Plaque_Volume")




residual_vs_baseline <- function(metric) {
  g <- df_wide %>%
    dplyr::filter(metric == !!metric, is.finite(Baseline), is.finite(`Follow-up`))
  if (!nrow(g)) return(ggplot2::ggplot() + ggplot2::theme_void())
  
  fit <- stats::lm(`Follow-up` ~ Baseline, data = g)
  g$resid <- stats::resid(fit)
  
  ggplot2::ggplot(g, ggplot2::aes(Baseline, resid)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.8) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.8) +
    ggplot2::labs(
      title = paste0(label_for(metric), " — Residual (FU − ŷ|BL) vs Baseline"),
      x = "Baseline", y = "Residual"
    ) + theme_sl
}



scatter_bf_with_stats <- function(metric,
                                  transform = c("auto","none","log1p"),
                                  fit_color = "steelblue",
                                  id_color  = "black",
                                  id_linetype = "dashed",
                                  point_alpha = 0.6,
                                  show_eq = TRUE) {
  
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
    b <- stats::coef(fit); A <- exp(b[1]); B <- b[2]
    eq_text <- sprintf("y + 1 = %.3g · (1 + x)^{%.2f}", A, B)
    subtitle <- sprintf("R²=%.2f; ICC=%.2f; fit: log1p, shown in original units", r2, icc)
  } else {
    fit <- stats::lm(`Follow-up` ~ Baseline, data = g)
    preds_grid <- stats::predict(fit, newdata = grid)
    r2 <- summary(fit)$r.squared
    b <- stats::coef(fit)
    eq_text <- sprintf("y = %.2f + %.2f · x", b[1], b[2])
    subtitle <- sprintf("R²=%.2f; ICC=%.2f", r2, icc)
  }
  grid$`Follow-up` <- preds_grid
  
  # Tighter limits; small pad
  lims <- range(c(g$Baseline, g$`Follow-up`), na.rm = TRUE)
  pad <- diff(lims) * 0.03
  lims <- c(lims[1] - pad, lims[2] + pad)
  
  # Build lines df so both appear in legend (no rbind name errors)
  types <- c("Model fit", "Identity (no change)")
  
  fit_df <- tibble::tibble(
    Baseline   = grid$Baseline,
    `Follow-up`= grid$`Follow-up`,
    type       = factor("Model fit", levels = types)
  )
  
  id_df <- tibble::tibble(
    Baseline   = lims,
    `Follow-up`= lims,
    type       = factor("Identity (no change)", levels = types)
  )
  
  lines_df <- dplyr::bind_rows(fit_df, id_df)
  
  p <- ggplot2::ggplot(g, ggplot2::aes(Baseline, `Follow-up`)) +
    ggplot2::geom_point(alpha = point_alpha) +
    ggplot2::geom_line(
      data = lines_df,
      ggplot2::aes(y = `Follow-up`, color = type, linetype = type),
      linewidth = 0.9
    ) +
    ggplot2::scale_color_manual(
      values = c("Model fit" = fit_color, "Identity (no change)" = id_color),
      breaks = types, name = NULL, drop = FALSE
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Model fit" = "solid", "Identity (no change)" = id_linetype),
      breaks = types, name = NULL, drop = FALSE
    ) +
    ggplot2::scale_x_continuous(limits = lims, expand = ggplot2::expansion(mult = c(0, 0.01))) +
    ggplot2::scale_y_continuous(limits = lims, expand = ggplot2::expansion(mult = c(0, 0.01))) +
    ggplot2::coord_fixed(ratio = 1, clip = "on") +
    ggplot2::labs(
      title = paste0(label_for(metric), " — Baseline vs Follow-up"),
      subtitle = subtitle, x = "Baseline", y = "Follow-up"
    ) +
    theme_sl +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      plot.title.position = "plot",
      plot.subtitle = ggplot2::element_text(margin = ggplot2::margin(t = 2, b = 6)),
      plot.margin = ggplot2::margin(t = 6, r = 6, b = 6, l = 6),
      panel.spacing = grid::unit(0, "pt")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
      linetype = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
    )
  
  if (isTRUE(show_eq)) {
    eq_x <- lims[1] + 0.02 * diff(lims)
    eq_y <- lims[2] - 0.02 * diff(lims)
    p <- p + ggplot2::annotate("text", x = eq_x, y = eq_y,
                               label = eq_text, hjust = 0, vjust = 1, size = 3.5)
  }
  
  p
}


p <- scatter_bf_with_stats(
  "Non_Calcified_Plaque_Volume",
  transform = "log1p",
  fit_color = "#2C7FB8",
  id_color = "gray30",
  show_eq = FALSE
)

p

save_plot_simple(p, "figures/plots_plaque_metrics/ncpv_fu_vs_bl_log.png", 6, 6, 600)


