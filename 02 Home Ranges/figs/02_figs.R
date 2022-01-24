# Figures for lecture 02 -- Home Ranges
# Making home range figures with same data as Fig. 2 in Signer and Fieberg 2021.

# Load packages----
library(tidyverse)
library(amt)
library(raster)
library(ragg)
library(patchwork)
library(here)

# Range vs Occurence ----

# Source function for BCRW
source("fun/bcrw.R")

# Centroids
centroids <- list("A01" = c("x" = 446589, "y" = 4625899),
                  "A02" = c("x" = 440284, "y" = 4618197),
                  "A03" = c("x" = 448796, "y" = 4613795),
                  "A04" = c("x" = 442367, "y" = 4627042)
)

# Set random seed
set.seed(20220126)

# Generate random walk for range distribution
rng <- lapply(centroids, function(cent) {
  # The basic BCRW
  x <- bcrw(start_loc = cent + rnorm(n = 2, mean = 0, sd = 500),
            centroid = cent,
            n_steps = 1000,
            sl_distr = c('shape' = 1, 'scale' = 300),
            rho = 0.50, # stronger correlation (was 0.25)
            beta = 0.15) # stronger bias (was 0.1)
  
  # Assign actual dates and times (start on June 20, 2021)
  x$date <- as.POSIXct("2021-06-20 00:00:00") +
    x$t * 60 * 60
  
  # Return
  return(x)
}) %>% 
  bind_rows(.id = "ID") %>% 
  # Keep just individual 4
  filter(ID == "A04") %>% 
  # Make track
  make_track(x, y, date, crs = 32612)

# Data for occurrence distribution
occ <- rng[1:300,]

# Fit KDE to occurrence
occ_kde <- occ %>% 
  hr_kde()

# Fit aKDE to occurrence
occ_akde <- occ %>% 
  hr_akde(model = fit_ctmm(occ, model = "ou"))

# Get polygons
kde_poly <- hr_isopleths(occ_kde) %>% 
  mutate(what = "Est", type = "KDE") %>% 
  dplyr::select(what, type)
akde_poly <- hr_isopleths(occ_akde) %>% 
  mutate(what = c("95% CI", "Est", "95% CI"), type = "aKDE") %>% 
  dplyr::select(what, type)
polys <- rbind(kde_poly, akde_poly)

# Make plots
(occ_plot <- ggplot() +
    geom_sf(data = kde_poly, color = "gray70", size = 1, fill = NA) +
    # Just to set the limits (100% transparent)
    geom_point(data = rng, aes(x = x_, y = y_), color = "#00000000") +
    geom_point(data = occ, aes(x = x_, y = y_), color = "blue", size = 1.25) +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw())

(occ_plot_long <- occ_plot +
    geom_point(data = rng, aes(x = x_, y = y_), color = "orange", size = 1) +
    geom_point(data = occ, aes(x = x_, y = y_), color = "blue", size = 1.25))

(rng_plot <- ggplot() +
    geom_sf(data = polys, aes(color = type, linetype = what), 
            size = 1, fill = NA) +
    geom_point(data = rng, aes(x = x_, y = y_), color = "orange", size = 1) +
    geom_point(data = occ, aes(x = x_, y = y_), color = "blue", size = 1.25) +
    xlab(NULL) +
    ylab(NULL) +
    scale_linetype_manual(name = "CI",
                          breaks = c("Est", "95% CI"), 
                          values = c("solid", "dashed")) +
    scale_color_manual(name = "Distribution",
                       breaks = c("aKDE", "KDE"),
                       labels = c("Range", "Occurrence"),
                       values = c("black", "gray70")) +
    theme_bw())

(dist_fig <- occ_plot + rng_plot)

ggsave(here::here("02 Home Ranges", "figs/occ_fig.png"), plot = occ_plot, 
       width = 800, height = 500, units = "px", dpi = 110)
ggsave(here::here("02 Home Ranges", "figs/occ_fig_long.png"), plot = occ_plot_long, 
       width = 800, height = 500, units = "px", dpi = 110)
ggsave(here::here("02 Home Ranges", "figs/rng_fig.png"), plot = rng_plot, 
       width = 1000, height = 500, units = "px", dpi = 110)
ggsave(here::here("02 Home Ranges", "figs/distr_fig.png"), plot = dist_fig, 
       width = 1000, height = 350, units = "px", dpi = 110)

# FISHER FIGURES ####
# Load data ----
dat <- amt_fisher %>% 
  filter(id == "F1")

# Fit home ranges ----
mcp <- hr_mcp(dat, levels = c(1, 0.95, 0.5))
locoh <- hr_locoh(dat, levels = seq(0.2, 1, by = 0.2),
                  type = "a", n = max(dist(dat[, c("x_", "y_")])))
kde <- hr_kde(dat, levels = c(0.95, 0.5))
akde <- hr_akde(dat, levels = c(0.95, 0.5), model = fit_ctmm(dat, model = "ou"))

# Figures ----

# ... locations ----
loc_fig <- dat %>% 
  ggplot(aes(x = x_, y = y_)) +
  geom_point(alpha = 0.3) +
  coord_sf(crs = 5070) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw()

ggsave(here::here("02 Home Ranges", "figs/locs.png"), plot = loc_fig, device = agg_png,
       width = 800, height = 750, units = "px", dpi = 150)

# ... MCP ----
mcp_fig <- ggplot() +
  geom_sf(data = mcp$mcp, aes(color = factor(level, 
                                             levels = c("1", "0.95", "0.5"),
                                             labels = c("100%", "95%", "50%"))), 
          fill = NA, size = 2) +
  geom_point(data = mcp$data, aes(x = x_, y = y_),
             color = "gray30", alpha = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  scale_color_brewer(name = "MCP Level",
                     palette = "Set1") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(colour = "black", size = 0.8))

ggsave(here::here("02 Home Ranges", "figs/mcp.png"), plot = mcp_fig, device = agg_png,
       width = 800, height = 800, units = "px", dpi = 150)

# ... LoCoH ----
locoh_fig <- ggplot() +
  geom_sf(data = locoh$locoh, aes(fill = fct_rev(factor(level)),
                                  color = fct_rev(factor(level))), 
          size = 1) +
  geom_point(data = locoh$data, aes(x = x_, y = y_),
             color = "gray", size = 0.5, alpha = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_brewer(name = "Level", palette = "Reds") +
  scale_color_brewer(name = "Level", palette = "Reds") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(colour = "black", size = 0.8))

ggsave(here::here("02 Home Ranges", "figs/locoh.png"), plot = locoh_fig, device = agg_png,
       width = 800, height = 800, units = "px", dpi = 150)

# ... KDE ----

# ... ... KDE UD ----
ud_fig <- kde$ud %>% 
  as.data.frame(xy = TRUE) %>% 
  ggplot(aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Probability\nDensity") +
  coord_sf(crs = 5070,
           xlim = range(dat$x_),
           ylim = range(dat$y_)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(colour = "black", size = 0.8))

ggsave(here::here("02 Home Ranges", "figs/ud.png"), plot = ud_fig, device = agg_png,
       width = 800, height = 800, units = "px", dpi = 150)

# ... ... KDE isopleths ----
kde_fig <- hr_isopleths(kde) %>% 
  # Bug where levels are reversed
  mutate(level = rev(level)) %>% 
  ggplot() +
  geom_sf(aes(color = factor(level, 
                             levels = c("0.95", "0.5"),
                             labels = c("95%", "50%"))), 
          fill = NA, size = 2) +
  geom_point(data = mcp$data, aes(x = x_, y = y_),
             color = "gray30", alpha = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  scale_color_brewer(name = "KDE Isopleth",
                     palette = "Set2") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(colour = "black", size = 0.8))

ggsave(here::here("02 Home Ranges", "figs/kde.png"), plot = kde_fig, device = agg_png,
       width = 800, height = 800, units = "px", dpi = 150)

# ... aKDE ----

# ... ... work around bug in hr_isopleth.akde() ----
hr_isopleths.akde <- function (x, conf.level = 0.95, ...) {
  checkmate::assert_number(conf.level, lower = 0, upper = 1)
  res <- ctmm::SpatialPolygonsDataFrame.UD(x$akde, level.UD = x$levels, 
                                           conf.level = conf.level)
  res1 <- sf::st_as_sf(res)
  res1 <- sf::st_transform(res1, akde$crs)
  
  ## Proposed fix
  # Extract level and what from raw names
  split_names <- strsplit(res1$name, " ", fixed = TRUE)
  
  level_perc <- sapply(split_names, getElement, 2)
  est <- sapply(split_names, getElement, 3)
  
  res1$level <- as.numeric(gsub(pattern = "%", replacement = "", 
                                x = level_perc, fixed = TRUE))/100
  res1$what <- ifelse(est == "low", paste0("lci (", conf.level, ")"),
                      ifelse(est == "est", "estimate",
                             paste0("uci (", conf.level, ")")))
  res1$name <- NULL
  row.names(res1) <- NULL
  ##
  
  res1$area = sf::st_area(res1)
  res1[, c("level", "what", "area", "geometry")]
}

# ... ... plot ----
akde_fig <- hr_isopleths(akde) %>% 
  filter(what == "estimate") %>% 
  mutate(level = factor(level, 
                        levels = c("0.95", "0.5"),
                        labels = c("95%", "50%"))) %>% 
  ggplot() +
  geom_sf(aes(color = level), 
          fill = NA, size = 2) +
  geom_point(data = akde$data, aes(x = x_, y = y_),
             color = "gray30", alpha = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  scale_color_brewer(name = "aKDE Isopleth",
                     palette = "Set2") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(colour = "black", size = 0.8))

ggsave(here::here("02 Home Ranges", "figs/akde.png"), plot = akde_fig, device = agg_png,
       width = 800, height = 800, units = "px", dpi = 150)
