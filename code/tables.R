# Keto CTA Plaque data
library(tidyverse)
library(stringr)
library(purrr)
library(gt)

df <- read_csv("data/keto-cta-quant-and-semi-quant.csv")

df$V1_Percent_Atheroma_Volume = df$V1_Percent_Atheroma_Volume*100
df$V2_Percent_Atheroma_Volume = df$V2_Percent_Atheroma_Volume*100


# discover paired baseline/follow-up variables
metrics <- names(df) |>
  keep(~ str_starts(.x, "V1_")) |>
  str_remove("^V1_") |>
  keep(~ paste0("V2_", .x) %in% names(df))

# labels and precision rules 
label_map <- c(
  "Total_Plaque_Score"            = "Total plaque score",
  "CAC"                           = "Coronary artery calcium score",
  "Non_Calcified_Plaque_Volume"   = "Noncalcified plaque volume (mm³)",
  "Total_Calcified_Plaque_Volume" = "Calcified plaque volume (mm³)",
  "Percent_Atheroma_Volume"       = "Percent atheroma volume (%)",
  "PAV_pct"                       = "Percent atheroma volume (%)"
)

digits_for <- function(metric_key) 1

# Drop-in fix: replace your fmt_num() with this base-R version (no 'scales' needed)
fmt_num <- function(x, d) {
  if (length(x) == 0 || all(!is.finite(x))) return(NA_character_)
  formatC(round(x, d), format = "f", digits = d, big.mark = ",")
}


fmt_mean_sd <- function(x, d) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_character_)
  paste0(fmt_num(mean(x), d), "\u00A0\u00B1\u00A0", fmt_num(sd(x), d))
}

fmt_med_iqr <- function(x, d) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_character_)
  qs <- as.numeric(quantile(x, probs = c(.25, .5, .75), na.rm = TRUE))
  paste0(fmt_num(qs[2], d), " (", fmt_num(qs[1], d), "–", fmt_num(qs[3], d), ")")
}

paired_n <- function(a, b) sum(complete.cases(a, b))

# build summary table with precision by metric 
summary_tbl <- map_dfr(metrics, function(m) {
  v1 <- df[[paste0("V1_", m)]]
  v2 <- df[[paste0("V2_", m)]]
  ch <- v2 - v1
  d  <- digits_for(m)
  
  tibble::tibble(
    Metric = dplyr::recode(m, !!!label_map, .default = m),
    N = paired_n(v1, v2),
    `Median (Q1–Q3): Baseline`  = fmt_med_iqr(v1, d),
    `Median (Q1–Q3): Follow-up` = fmt_med_iqr(v2, d),
    `Median (Q1–Q3): Change`    = fmt_med_iqr(ch, d),
    `Mean ± SD: Baseline`       = fmt_mean_sd(v1, d),
    `Mean ± SD: Follow-up`      = fmt_mean_sd(v2, d),
    `Mean ± SD: Change`         = fmt_mean_sd(ch, d)
  )
}) %>%
  arrange(Metric)

desired_order <- c(
  "Noncalcified plaque volume (mm³)",
  "Coronary artery calcium score",
  "Percent atheroma volume (%)",
  "Total plaque score",
  "Calcified plaque volume (mm³)"
)

summary_tbl <- summary_tbl %>%
  dplyr::mutate(.ord = match(Metric, desired_order)) %>%
  dplyr::arrange(is.na(.ord), .ord) %>%
  dplyr::select(-.ord)

# presentation table (clean headers, striping, separators, divider, notes)
gt_tbl_plaque <- summary_tbl %>%
  gt() %>%
  # short child headers under the spanners
  cols_label(
    `Median (Q1–Q3): Baseline`  = "Baseline",
    `Median (Q1–Q3): Follow-up` = "Follow-up",
    `Median (Q1–Q3): Change`    = "Change",
    `Mean ± SD: Baseline`       = "Baseline",
    `Mean ± SD: Follow-up`      = "Follow-up",
    `Mean ± SD: Change`         = "Change"
  ) %>%
  tab_spanner(
    label = "Median (Q1–Q3)",
    columns = c(`Median (Q1–Q3): Baseline`,
                `Median (Q1–Q3): Follow-up`,
                `Median (Q1–Q3): Change`)
  ) %>%
  tab_spanner(
    label = "Mean ± SD",
    columns = c(`Mean ± SD: Baseline`,
                `Mean ± SD: Follow-up`,
                `Mean ± SD: Change`)
  ) %>%
  # alignment
  cols_align(align = "left", columns = Metric) %>%
  cols_align(
    align = "center",
    columns = c(
      N,
      `Median (Q1–Q3): Baseline`, `Median (Q1–Q3): Follow-up`, `Median (Q1–Q3): Change`,
      `Mean ± SD: Baseline`, `Mean ± SD: Follow-up`, `Mean ± SD: Change`
    )
  ) %>%
  # subtle row striping
  opt_row_striping() %>%
  tab_options(
    row.striping.background_color = "#E7EFF7",  # light blue stripe
    row.striping.include_stub = TRUE            # include the Metric stub column
  ) %>%
  # vertical divider between stat blocks (apply to header AND body)
  tab_style(
    style = cell_borders(sides = "left", color = "#D0D0D0", weight = px(2)),
    locations = list(
      cells_column_labels(columns = `Mean ± SD: Baseline`),
      cells_body(columns = `Mean ± SD: Baseline`)
    )
  ) %>%
  # prevent wrapping anywhere in the table
  tab_style(
    style = cell_text(whitespace = "nowrap"),
    locations = list(
      cells_body(),
      cells_column_labels(),
      cells_stub(),
      cells_stubhead(),
      cells_column_spanners()
    )
  ) %>%
  # heading (styled like a caption): white text, larger, no-wrap, dark bg
  tab_header(title = "Baseline, follow-up, and paired change in plaque metrics at 1 year") %>%
  tab_style(
    style = list(
      cell_text(color = "white", size = px(22), whitespace = "nowrap"),
      cell_fill(color = "#000000")
    ),
    locations = cells_title(groups = "title")
  ) %>%
  # footnotes
  tab_footnote(
    footnote = "Change = follow-up − baseline, summarized across participants.",
    locations = list(
      cells_column_labels(columns = `Median (Q1–Q3): Change`),
      cells_column_labels(columns = `Mean ± SD: Change`)
    )
  ) %>%
  tab_footnote(
    footnote = "N = number of paired scans per metric.",
    locations = cells_column_labels(columns = N)
  ) %>%
  tab_options(table.font.size = px(16), data_row.padding = px(6))

gt_tbl_plaque <- gt_tbl_plaque %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_labels()) %>%
  tab_options(data_row.padding = px(8))


###############################################
## actual percent change (standard version)  ##
###############################################

#  helper: count usable pairs for % change (non-missing, nonzero baseline)
paired_n_pct <- function(a, b) sum(is.finite(a) & is.finite(b) & a != 0)

pct_change_tbl <- map_dfr(metrics, function(m) {
  v1 <- df[[paste0("V1_", m)]]
  v2 <- df[[paste0("V2_", m)]]
  pc <- ifelse(is.finite(v1) & is.finite(v2) & v1 != 0, (v2 - v1) / v1 * 100, NA_real_)
  d  <- 1  # % change precision
  
  tibble::tibble(
    Metric = dplyr::recode(m, !!!label_map, .default = m),
    `N (nonzero baseline)` = paired_n_pct(v1, v2),
    `Median (Q1–Q3): % change` = fmt_med_iqr(pc, d),
    `Mean ± SD: % change`      = fmt_mean_sd(pc, d)
  )
}) %>%
  arrange(Metric) %>%
  mutate(.ord = match(Metric, desired_order)) %>%
  arrange(is.na(.ord), .ord) %>%
  select(-.ord)

gt_tbl_plaque_pct <- pct_change_tbl %>%
  gt() %>%
  cols_label(
    `Median (Q1–Q3): % change` = "Median (Q1–Q3)",
    `Mean ± SD: % change`      = "Mean ± SD"
  ) %>%
  tab_spanner(
    label  = "Percent change from baseline",
    columns = c(`Median (Q1–Q3): % change`, `Mean ± SD: % change`)
  ) %>%
  cols_align(align = "left", columns = Metric) %>%
  cols_align(align = "center",
             columns = c(`N (nonzero baseline)`,
                         `Median (Q1–Q3): % change`,
                         `Mean ± SD: % change`)) %>%
  opt_row_striping() %>%
  tab_options(
    row.striping.background_color = "#E7EFF7",
    row.striping.include_stub = TRUE,
    table.font.size = px(16),
    data_row.padding = px(8)
  ) %>%
  tab_style(
    style = cell_text(whitespace = "nowrap"),
    locations = list(
      cells_body(),
      cells_column_labels(),
      cells_stub(),
      cells_stubhead(),
      cells_column_spanners()
    )
  ) %>%
  tab_header(title = "Relative percent change from baseline to follow-up in plaque metrics at 1 year") %>%
  tab_style(
    style = list(
      cell_text(color = "white", size = px(22), whitespace = "nowrap"),
      cell_fill(color = "#000000")
    ),
    locations = cells_title(groups = "title")
  ) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_labels()) %>%
  tab_footnote(
    footnote = "Percent change = (follow-up − baseline) / baseline × 100; calculated per participant and summarized across participants.",
    locations = list(cells_column_spanners(spanners = "Percent change from baseline"))
  ) %>%
  tab_footnote(
    footnote = "N counts participants with a nonzero baseline for the metric.",
    locations = cells_column_labels(columns = `N (nonzero baseline)`)
  )


# ---- Zeros/incidence/resolution (counts) ----
zero_summary_tbl <- purrr::map_dfr(metrics, function(m) {
  v1 <- df[[paste0("V1_", m)]]
  v2 <- df[[paste0("V2_", m)]]
  keep <- complete.cases(v1, v2)
  
  tibble::tibble(
    Metric = dplyr::recode(m, !!!label_map, .default = m),
    N      = sum(keep),
    `Baseline = 0`   = fmt_num(sum(v1[keep] == 0), 0),
    `Follow-up = 0`  = fmt_num(sum(v2[keep] == 0), 0),
    `0 → >0`         = fmt_num(sum(v1[keep] == 0 & v2[keep] > 0), 0),   # incidence
    `>0 → 0`         = fmt_num(sum(v1[keep] > 0 & v2[keep] == 0), 0)    # resolution
  )
}) %>%
  dplyr::mutate(.ord = match(Metric, desired_order)) %>%
  dplyr::arrange(is.na(.ord), .ord) %>%
  dplyr::select(-.ord)

gt_tbl_zero <- zero_summary_tbl %>%
  gt() %>%
  cols_label(
    `Baseline = 0`  = "Baseline = 0",
    `Follow-up = 0` = "Follow-up = 0",
    `0 → >0`        = " 0 → >0 ",
    `>0 → 0`        = " >0 → 0 "
  ) %>%
  opt_row_striping() %>%
  cols_align("left", Metric) %>%
  cols_align("center", c(`N`, `Baseline = 0`, `Follow-up = 0`, `0 → >0`, `>0 → 0`)) %>%
  tab_header(title = "Zeros, incidence (0→>0), and resolution (>0→0) by metric") %>%
  tab_style(style = cell_text(color = "white", size = px(22)),
            locations = cells_title(groups = "title")) %>%
  tab_style(style = cell_fill(color = "#000000"),
            locations = cells_title(groups = "title")) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_labels()) %>%
  tab_options(row.striping.background_color = "#E7EFF7",
              row.striping.include_stub = TRUE,
              table.font.size = px(16),
              data_row.padding = px(8))

gt_tbl_zero <- gt_tbl_zero %>%
  tab_style(
    style = cell_borders(sides = "left", color = "#FFFFFF", weight = px(4)),
    locations = cells_column_labels(columns = c(`N`,`Baseline = 0`,`Follow-up = 0`, `0 → >0`, `>0 → 0`))
  )

# Square-root change (works when x ≥ 0)
sqrt_change_tbl <- map_dfr(metrics, function(m) {
  v1 <- df[[paste0("V1_", m)]]
  v2 <- df[[paste0("V2_", m)]]
  dval <- sqrt(pmax(v2,0)) - sqrt(pmax(v1,0))
  tibble(
    Metric = dplyr::recode(m, !!!label_map, .default = m),
    N = paired_n(v1, v2),
    `Median (Q1–Q3): Δ√x` = fmt_med_iqr(dval, 2),
    `Mean ± SD: Δ√x`      = fmt_mean_sd(dval, 2)
  )
})

# Log-ratio with offset c
log_ratio_tbl <- map_dfr(metrics, function(m) {
  v1 <- df[[paste0("V1_", m)]]
  v2 <- df[[paste0("V2_", m)]]
  c <- 1  # choose per metric if needed
  lr <- 100 * log((v2 + c) / (v1 + c))
  tibble(
    Metric = dplyr::recode(m, !!!label_map, .default = m),
    N = paired_n(v1, v2),
    `Median (Q1–Q3): 100·log((FU+c)/(BL+c))` = fmt_med_iqr(lr, 1),
    `Mean ± SD: 100·log((FU+c)/(BL+c))`      = fmt_mean_sd(lr, 1)
  )
})





min_pos <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (!length(x)) NA_real_ else min(x)
}

# metric-aware offset (half the smallest nonzero; sensible fallbacks)
offset_for <- function(v1, v2, metric) {
  mp <- min_pos(c(v1, v2))
  c <- 0.5 * mp
  if (!is.finite(c) || c == 0) {
    if (metric %in% c("Percent_Atheroma_Volume","PAV_pct")) 0.1 else 1
  } else c
}

pct_change_offset_tbl <- map_dfr(metrics, function(m) {
  v1 <- df[[paste0("V1_", m)]]
  v2 <- df[[paste0("V2_", m)]]
  c0 <- offset_for(v1, v2, m)
  
  r  <- (v2 + c0) / (v1 + c0)
  pc <- (r - 1) * 100   # back-transformed, % scale
  d  <- 1
  
  tibble::tibble(
    Metric = dplyr::recode(m, !!!label_map, .default = m),
    N = paired_n(v1, v2),
    `Median (Q1–Q3): % change` = fmt_med_iqr(pc, d),
    `Mean ± SD: % change`      = fmt_mean_sd(pc, d)
  )
}) %>%
  mutate(.ord = match(Metric, desired_order)) %>%
  arrange(is.na(.ord), .ord) %>%
  select(-.ord)

gt_tbl_pct_offset <- pct_change_offset_tbl %>%
  gt() %>%
  cols_label(
    `Median (Q1–Q3): % change` = "Median (Q1–Q3)",
    `Mean ± SD: % change`      = "Mean ± SD"
  ) %>%
  tab_spanner(
    label  = "Offset-adjusted percent change (includes baseline = 0)",
    columns = c(`Median (Q1–Q3): % change`, `Mean ± SD: % change`)
  ) %>%
  cols_align("left", Metric) %>%
  cols_align("center", c(N, `Median (Q1–Q3): % change`, `Mean ± SD: % change`)) %>%
  opt_row_striping() %>%
  tab_header(title = "Offset-adjusted percent change from baseline to follow-up") %>%
  tab_style(
    style = list(cell_text(color = "white", size = px(22), whitespace = "nowrap"),
                 cell_fill(color = "#000000")),
    locations = cells_title(groups = "title")
  ) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels()) %>%
  tab_options(row.striping.background_color = "#E7EFF7",
              row.striping.include_stub = TRUE,
              table.font.size = px(16),
              data_row.padding = px(8)) %>%
  tab_style(style = cell_text(whitespace = "nowrap"),
            locations = list(cells_body(), cells_column_labels(),
                             cells_stub(), cells_stubhead(), cells_column_spanners())) %>%
  tab_footnote(
    footnote = "Percent change = ((follow-up + c) / (baseline + c) − 1) × 100; c = half the smallest nonzero across visits for that metric (fallback c = 1; for percentages c = 0.1).",
    locations = list(cells_column_spanners(spanners = "Offset-adjusted percent change (includes baseline = 0)"))
  ) %>%
  tab_footnote(
    footnote = "Includes all paired participants, including baseline = 0.",
    locations = cells_column_labels(columns = N)
  )



###################
# --- helpers ---
tmean <- function(x, trim = 0.1) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x, trim = trim)
}

robust_offset_for <- function(v1, v2, metric) {
  nz <- c(v1, v2); nz <- nz[is.finite(nz) & nz > 0]
  # half the 5th percentile, with sensible floors
  if (!length(nz)) return(if (metric %in% c("Percent_Atheroma_Volume","PAV_pct")) 0.5 else if (metric == "CAC") 10 else 5)
  c_est <- 0.5 * as.numeric(quantile(nz, 0.05, na.rm = TRUE))
  if (metric %in% c("Percent_Atheroma_Volume","PAV_pct")) max(c_est, 0.5)
  else if (metric == "CAC") max(c_est, 10)
  else max(c_est, 5)
}

# --- hybrid % change table: two medians + trimmed mean ---
pct_change_hybrid_tbl <- purrr::map_dfr(metrics, function(m) {
  v1 <- df[[paste0("V1_", m)]]
  v2 <- df[[paste0("V2_", m)]]
  c0 <- robust_offset_for(v1, v2, m)
  
  pc_all <- ((v2 + c0) / (v1 + c0) - 1) * 100            # includes zeros
  keep_nz <- is.finite(v1) & is.finite(v2) & v1 > 0
  pc_nz  <- ifelse(keep_nz, (v2 - v1) / v1 * 100, NA_real_)
  
  d <- 1
  tibble::tibble(
    Metric = dplyr::recode(m, !!!label_map, .default = m),
    N = paired_n(v1, v2),
    `N baseline > 0` = sum(keep_nz, na.rm = TRUE),
    `Median % chg (overall, offset)` = fmt_med_iqr(pc_all, d),
    `Median % chg (BL>0)`            = fmt_med_iqr(pc_nz, d),
    `Trimmed mean % chg (overall, 10%)` = fmt_num(tmean(pc_all, 0.1), d)
  )
}) %>%
  dplyr::mutate(.ord = match(Metric, desired_order)) %>%
  dplyr::arrange(is.na(.ord), .ord) %>%
  dplyr::select(-.ord)

gt_tbl_pct_hybrid <- pct_change_hybrid_tbl %>%
  gt() %>%
  cols_label(
    `Median % chg (overall, offset)` = "Overall (offset)",
    `Median % chg (BL>0)`            = "Baseline > 0",
    `Trimmed mean % chg (overall, 10%)` = "Trimmed mean, 10%"
  ) %>%
  tab_spanner(
    label = "Median (Q1–Q3) % change",
    columns = c(`Median % chg (overall, offset)`, `Median % chg (BL>0)`)
  ) %>%
  cols_align("left", Metric) %>%
  cols_align(
    "center",
    c(
      N, `N baseline > 0`,
      `Median % chg (overall, offset)`,
      `Median % chg (BL>0)`,
      `Trimmed mean % chg (overall, 10%)`
    )
  ) %>%
  opt_row_striping() %>%
  tab_header(title = "Percent change from baseline to follow-up: medians and trimmed mean") %>%
  tab_style(
    style = list(cell_text(color = "white", size = px(22), whitespace = "nowrap"),
                 cell_fill(color = "#000000")),
    locations = cells_title(groups = "title")
  ) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels()) %>%
  tab_options(
    row.striping.background_color = "#E7EFF7",
    row.striping.include_stub = TRUE,
    table.font.size = px(16),
    data_row.padding = px(8)
  ) %>%
  tab_style(
    style = cell_text(whitespace = "nowrap"),
    locations = list(
      cells_body(), cells_column_labels(),
      cells_stub(), cells_stubhead(), cells_column_spanners()
    )
  ) %>%
  tab_footnote(
    footnote =
      "Overall % change = ((follow-up + c)/(baseline + c) - 1) × 100. BL>0 % change = ((follow-up - baseline)/baseline) × 100. c = 0.5 × the 5th percentile of non-zero values across visits for that metric, with floors: percentages 0.5; CAC 10; volumes 5. Positive values indicate increase vs baseline; negative indicate decrease.",
    locations = list(cells_column_spanners(spanners = "Median (Q1–Q3) % change"))
  ) %>%
  tab_footnote(
    footnote = "Trimmed mean removes 10% from each tail before averaging; computed on the Overall (offset) % change.",
    locations = cells_column_labels(columns = `Trimmed mean % chg (overall, 10%)`)
  ) %>%
  tab_footnote(
    footnote = "N = paired participants; N baseline > 0 = number with nonzero baseline.",
    locations = cells_column_labels(columns = c(N, `N baseline > 0`))
  )



####################
# Abbreviations for display (keep the table order)
abbr_map <- c(
  "Noncalcified plaque volume (mm³)" = "NCPV",
  "Coronary artery calcium score"    = "CAC",
  "Percent atheroma volume (%)"      = "PAV",
  "Total plaque score"               = "TPS",
  "Calcified plaque volume (mm³)"    = "TCPV"
)
desired_order_abbr <- unname(abbr_map[desired_order])

# Build transition counts
zero_flow_tbl <- purrr::map_dfr(metrics, function(m) {
  v1 <- df[[paste0("V1_", m)]]
  v2 <- df[[paste0("V2_", m)]]
  keep <- complete.cases(v1, v2)
  n <- sum(keep)
  
  tibble::tibble(
    Metric_full = dplyr::recode(m, !!!label_map, .default = m),
    N_paired = n,
    `Stayed 0`   = sum(v1[keep] == 0 & v2[keep] == 0),
    `0 → >0`     = sum(v1[keep] == 0 & v2[keep] > 0),
    `>0 → 0`     = sum(v1[keep] > 0 & v2[keep] == 0),
    `Stayed >0`  = sum(v1[keep] > 0 & v2[keep] > 0)
  )
}) %>%
  # order by your table, then abbreviate
  dplyr::mutate(Metric_full = factor(Metric_full, levels = desired_order),
                Metric = dplyr::recode(Metric_full, !!!abbr_map)) %>%
  tidyr::pivot_longer(
    cols = c(`Stayed 0`, `0 → >0`, `>0 → 0`, `Stayed >0`),
    names_to = "Transition", values_to = "Count"
  ) %>%
  dplyr::mutate(
    Metric = factor(Metric, levels = desired_order_abbr),
    Transition = factor(Transition,
                        levels = c("Stayed 0", "0 → >0", ">0 → 0", "Stayed >0")),
    Share = Count / N_paired  # 0–1 ensures bars never spill past the edge
  )

zero_flow_tbl <- zero_flow_tbl %>%
  dplyr::mutate(
    Transition = factor(
      Transition,
      levels = c("Stayed >0", "0 → >0", ">0 → 0", "Stayed 0")
    )
  )

zero_flow_tbl <- zero_flow_tbl %>%
  dplyr::mutate(
    Transition = factor(Transition, levels = c("Stayed >0", "0 → >0", ">0 → 0", "Stayed 0")),
    # 2) Reverse the metric order so the first (NCPV) appears at the TOP after coord_flip()
    Metric = factor(Metric, levels = rev(levels(Metric)))
  )

fill_cols <- c(
  "Stayed >0" = "#8C2D2D",  # muted dark red (worse)
  "0 → >0"    = "#B9832F",  # muted orange (worse)
  ">0 → 0"    = "#A6C8A0",  # muted dark green (good)
  "Stayed 0"  = "#2F6B4F"   # muted light green (good)
)


# Plot with COUNTS
zero_transition_plot = ggplot(zero_flow_tbl, aes(x = Metric, y = Count, fill = Transition)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, max(zero_flow_tbl$N_paired, na.rm = TRUE))
  ) +
  scale_fill_manual(
    values = fill_cols,
    breaks = c("Stayed >0", "0 → >0", ">0 → 0", "Stayed 0")  # legend order
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, reverse = TRUE)) +  # count labels; hide tiny slices
  geom_text(
    aes(label = ifelse(Count >= 3, Count, "")),
    position = position_stack(vjust = 0.5),
    color = "white", size = 4
  ) +
  labs(
    title = "Plaque Metrics: Zero transitions from baseline to follow-up",
    y = "Participants", x = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.caption = element_text(hjust = 0.5),
    plot.margin = margin(t = 10, r = 16, b = 10, l = 16)   # outer margin
    
  )


# save_plot_simple(
#  zero_transition_plot,
#  "figures/plots_plaque_metrics/zero_transition_plot.png",
#  10, 6, 1000
# )


# Build counts of change categories (exact, no buffer)
change_flow_tbl <- purrr::map_dfr(metrics, function(m) {
  v1 <- df[[paste0("V1_", m)]]
  v2 <- df[[paste0("V2_", m)]]
  keep <- complete.cases(v1, v2)
  dlt  <- v2[keep] - v1[keep]
  tibble::tibble(
    Metric_full = dplyr::recode(m, !!!label_map, .default = m),
    N_paired = length(dlt),
    Decreased = sum(dlt < 0, na.rm = TRUE),
    `No change` = sum(dlt == 0, na.rm = TRUE),
    Increased = sum(dlt > 0, na.rm = TRUE)
  )
}) %>%
  dplyr::mutate(
    Metric_full = factor(Metric_full, levels = desired_order),
    Metric = dplyr::recode(Metric_full, !!!abbr_map)
  ) %>%
  tidyr::pivot_longer(
    cols = c(Decreased, `No change`, Increased),
    names_to = "Change", values_to = "Count"
  ) %>%
  dplyr::mutate(
    # Left→right stack order: Increased (worse), No change, Decreased (good)
    Change = factor(Change, levels = c("Increased", "No change", "Decreased")),
    # Put NCPV at the TOP after coord_flip()
    Metric = factor(Metric, levels = rev(desired_order_abbr))
  )

# Muted, logical colors + readable label colors
fill_cols <- c(
  "Increased" = "#8C2D2D",  # muted dark red (worse)
  "No change" = "#B8B8B8",  # neutral gray
  "Decreased" = "#2F6B4F"   # muted dark green (good)
)
label_cols <- c(
  "Increased" = "white",
  "No change" = "black",
  "Decreased" = "white"
)

# Headroom so bars don't touch the right edge
ymax <- max(change_flow_tbl$N_paired, na.rm = TRUE)
pad  <- max(2, ceiling(0.02 * ymax))

change_category_plot = ggplot(change_flow_tbl, aes(x = Metric, y = Count, fill = Change)) +
  geom_col(width = 0.7) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, ymax ), expand = c(0, 0)) +
  scale_fill_manual(values = fill_cols, breaks = levels(change_flow_tbl$Change)) +
  # Count labels; hide tiny slices
  geom_text(
    aes(label = ifelse(Count >= 2, Count, ""), color = Change),
    position = position_stack(vjust = 0.5),
    size = 4, show.legend = FALSE
  ) +
  scale_color_manual(values = label_cols, guide = "none") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, reverse = TRUE)) +  # count labels; hide tiny slices
  labs(
    title = "Plaque Metrics: Change from baseline to follow-up by category",
    y = "Participants", x = NULL, fill = NULL,
    caption = "Exact categories: Decreased (<0), No change (=0), Increased (>0)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.caption = element_text(hjust = 0.5),
    plot.margin = margin(10, 16, 10, 16)
  )

# save_plot_simple(
#  change_category_plot,
#  "figures/plots_plaque_metrics/change_category_plot.png",
#  10, 6, 1000
# )
