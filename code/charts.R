# Charts

# Keto CTA Plaque data
library(tidyverse)
library(ragg)
library(ggplot2)
library(plotly)
library(ggtext)

df <- read_csv("data/keto-cta-quant-and-semi-quant.csv")

# create change scores for each Plaque metric
df <- df %>%
  mutate(
    delta_NCPV = V2_Non_Calcified_Plaque_Volume - V1_Non_Calcified_Plaque_Volume,
    delta_TPS  = V2_Total_Plaque_Score - V1_Total_Plaque_Score,
    delta_PAV  = V2_Percent_Atheroma_Volume - V1_Percent_Atheroma_Volume,
    delta_PAV_pct  = (V2_Percent_Atheroma_Volume - V1_Percent_Atheroma_Volume)*100,
    V1_PAV_pct = V1_Percent_Atheroma_Volume * 100
  )

##########

#The median change in NCPV was 18.9 mm3 (IQR: 9.3-47.0 mm3) 

quantile(df$V1_Non_Calcified_Plaque_Volume, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
IQR(df$V1_Non_Calcified_Plaque_Volume, na.rm = TRUE)

quantile(df$V2_Non_Calcified_Plaque_Volume, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
IQR(df$V2_Non_Calcified_Plaque_Volume, na.rm = TRUE)

quantile(df$delta_NCPV, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
IQR(df$delta_NCPV, na.rm = TRUE)


# function for main chart

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


figure_1A = plot_two_timepoints(df,
    V1_Non_Calcified_Plaque_Volume, V2_Non_Calcified_Plaque_Volume,
    y_lab = "NCPV mm<sup>3</sup>", title = "ΔNCPV (mm<sup>3</sup>)",
    label1 = "Baseline", label2 = "1 Year"
)

# figure 1A

gg_figure_1A = ggplotly(figure_1A) %>%
  config(displayModeBar = FALSE)
gg_figure_1A

# Figure 1B

figure_1B = plot_two_timepoints(df,
                                V1_Percent_Atheroma_Volume, V2_Percent_Atheroma_Volume,
                                y_lab = "PAV %", title = "Proportional NCP and CP by PAV",
                                label1 = "Baseline", label2 = "1 Year",
                                scale_fn = ~ .x * 100 # convert decimal to %
                                
)

gg_figure_1B = ggplotly(figure_1B) %>%
  config(displayModeBar = FALSE)


#### plot 2F

figure_2F = ggplot(df, aes(x = V1_CAC, y = delta_TPS)) +
  geom_point(alpha = 0.8, size = 1.6, na.rm = TRUE) +
  stat_smooth(method = "lm", formula = y ~ x,
              se = TRUE, level = 0.95,
              color = "red", fill = "grey70", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(x = "Baseline CAC", y = "ΔTotal Plaque Score") +
  scale_y_continuous(breaks = -1:6, minor_breaks = NULL) +
  coord_cartesian(ylim = c(-1, 6)) +   # keeps points & fit intact
  theme_classic(base_size = 16) +
  theme(
    
  axis.text.y = element_text(face = "bold"),
  axis.text.x  = element_text(face = "bold"),
  )

figure_2F

#### Delta NCPV ~ baseline NCPV
ggplot(df, aes(x = V1_Non_Calcified_Plaque_Volume, y = delta_NCPV)) +
  geom_point(alpha = 0.8, size = 1.6, na.rm = TRUE) +
  stat_smooth(method = "lm", formula = y ~ x,
              se = TRUE, level = 0.95,
              color = "red", fill = "grey70", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(x = "Baseline NCPV", y = "ΔNCPV", title = "Baseline NCPV vs ΔNCPV") +
  #scale_y_continuous(breaks = -1:6, minor_breaks = NULL) +
  #coord_cartesian(ylim = c(-1, 6)) +   # keeps points & fit intact
  theme_classic(base_size = 16) +
  theme(
    
    axis.text.y = element_text(face = "bold"),
    axis.text.x  = element_text(face = "bold"),
  )

# save figs as png
save_png <- function(p, file, width = 6, height = 4, dpi = 600, bg = "white") {
  ggsave(file, plot = p, device = ragg::agg_png,
         width = width, height = height, units = "in",
         dpi = dpi, bg = bg)
}

save_png(figure_1A, "figures/Figure1A.png", width = 6, height = 5.86, dpi = 800)
save_png(figure_1B, "figures/Figure1B.png", width = 6, height = 6.23, dpi = 800)
save_png(figure_2F, "figures/Figure2F.png", width = 6, height = 6.26, dpi = 800)

#ggsave("figures/Figure1A.pdf", plot = figure_1A, device = "pdf",
#       width = 6, height = 12, units = "in")

ggsave("figures/Figure1A.svg", plot = figure_1A, device = svglite::svglite,
       width = 6, height = 12, units = "in")
ggsave("figures/Figure1B.svg", plot = figure_1B, device = svglite::svglite,
       width = 6, height = 12, units = "in")
ggsave("figures/Figure2F.svg", plot = figure_2F, device = svglite::svglite,
       width = 6, height = 6, units = "in")
