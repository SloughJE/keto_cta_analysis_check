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
    Metric = recode(m, !!!label_map, .default = m),
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


