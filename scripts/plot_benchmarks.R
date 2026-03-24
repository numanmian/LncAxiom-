#!/usr/bin/env Rscript

# scripts/plot_benchmarks.R
# Publication‑ready benchmark plot with panel labels (i) to (iv)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(viridis)
  library(extrafont)
  library(glue)   # for flexible string interpolation
})

# Try to load Times New Roman (if available)
tryCatch({
  font_import(pattern = "times", prompt = FALSE)
  loadfonts(quiet = TRUE)
  message("Times New Roman font loaded.")
}, error = function(e) {
  message("Times New Roman not available, using default sans font.")
})

# -------------------------------------------------------------------
# 1. Find and combine all benchmark files
# -------------------------------------------------------------------
benchmark_files <- list.files(path = "benchmarks",
                               pattern = "\\.txt$",
                               recursive = TRUE,
                               full.names = TRUE)

if (length(benchmark_files) == 0) {
  stop("No benchmark files found. Run the pipeline first to generate benchmarks.")
}

benchmarks <- benchmark_files %>%
  set_names() %>%
  map_dfr(read_tsv, col_types = cols(.default = "d"), .id = "file") %>%
  mutate(
    rule = str_extract(file, "(?<=benchmarks/)[^/]+")
  )

# Replace zeros with a small positive number to avoid log(0)
plot_data <- benchmarks %>%
  select(file, rule, s, max_rss, io_in, io_out) %>%
  pivot_longer(cols = c(io_in, io_out, max_rss, s),
               names_to = "metric",
               values_to = "value") %>%
  mutate(
    value = ifelse(value == 0, 0.01, value),
    metric = case_when(
      metric == "io_in"   ~ "I/O (in), MB",
      metric == "io_out"  ~ "I/O (out), MB",
      metric == "max_rss" ~ "RSS (max), MB",
      metric == "s"       ~ "Run time, s"
    ),
    metric = factor(metric, levels = c("I/O (in), MB", "I/O (out), MB",
                                        "RSS (max), MB", "Run time, s"))
  )

# -------------------------------------------------------------------
# 2. Custom labeller to add Roman numerals (i), (ii), (iii), (iv)
#    to each facet panel.
# -------------------------------------------------------------------
# Get the unique metrics (they are already in the desired order)
metrics_levels <- levels(plot_data$metric)

# Create a named vector mapping each metric to its panel label
# We use Roman numerals: i, ii, iii, iv
panel_labels <- setNames(
  glue("({as.roman(seq_along(metrics_levels))}) {metrics_levels}"),
  metrics_levels
)

# A labeller function that uses this mapping
panel_labeller <- labeller(metric = panel_labels)

# -------------------------------------------------------------------
# 3. Create the plot
# -------------------------------------------------------------------
p <- ggplot(plot_data, aes(x = rule, y = value, fill = rule, color = rule)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, show.legend = FALSE) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.8, show.legend = FALSE) +
  facet_wrap(~metric, scales = "free_y", ncol = 2, labeller = panel_labeller) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.1),
    breaks = scales::trans_breaks("log10", function(x) 10^x)
  ) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.9) +
  scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.9) +
  labs(
    x = "Pipeline rule",
    y = "Value (log scale)",
    title = "Computational resource usage of LncPro pipeline steps"
  ) +
  theme_bw(base_size = 12, base_family = "Times New Roman") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# -------------------------------------------------------------------
# 4. Save output
# -------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

ggsave("figures/resource_usage.png", p, width = 12, height = 9, dpi = 600)
ggsave("figures/resource_usage.pdf", p, width = 12, height = 9, device = cairo_pdf)

cat("Benchmark plot saved to figures/resource_usage.png and .pdf (600 dpi)\n")
cat("Panel labels: (i) I/O in, (ii) I/O out, (iii) RSS, (iv) Run time\n")
