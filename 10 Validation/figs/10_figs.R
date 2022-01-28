# Figures for lecture 10 -- Validation

# Load packages----
library(tidyverse)
library(amt)
library(raster)
library(ragg)
library(here)

# Biased counts ----
set.seed(123456789)
count_data <- data.frame(
  # True abundance in site
  true = round(runif(n = 50, 10, 60)),
  # Detection probability in site
  p = rbeta(n = 50, shape1 = 35, shape2 = 25)) %>% 
  # Count
  mutate(count = rbinom(n = 50, size = true, prob = p))

(site_counts <- ggplot(count_data, aes(x = count, y = true)) +
  geom_point(size = 0.9) +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  coord_cartesian(xlim = c(0, 60),
                  ylim = c(0, 60)) +
    xlab("Population Index") +
    ylab("True Abundance") +
  theme_bw())

ggsave(here("10 Validation", "figs/site_counts.png"), plot = site_counts, 
       device = agg_png, width = 700, height = 500, units = "px", dpi = 200)
