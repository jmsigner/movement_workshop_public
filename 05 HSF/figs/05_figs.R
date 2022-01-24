# Figures for lecture 05 -- HSFs

# Load packages----
library(tidyverse)
library(amt)
library(raster)
library(patchwork)
library(ragg)
library(mvtnorm)
library(here)

# E-space vs G-space ----

# ... G-space ----
# Coordinates
coords <- expand.grid("x" = seq(-1000, 950, by = 50),
                      "y" = seq(-1000, 950, by = 50)) %>% 
  arrange(x, y) %>% 
  as.matrix()

dm <- as.matrix(dist(coords))
dm <- dm/max(dm)

# Covariance matrix
# Exponential covariance function
expcov <- function(dists, rho, sigma) {
  n <- dim(dists)[1]
  result <- matrix(nrow = n, ncol = n)
  sigma2 <- sigma^2
  result <- sigma2 * exp(-1 * dists/rho)
  return(result)
}


set.seed(20220105)
temp <- rmvnorm(n = 1,
                mean = rep(25, nrow(coords)),
                sigma = expcov(dm, 5, 2))
precip <- rmvnorm(n = 1,
                  mean = rep(450, nrow(coords)),
                  sigma = expcov(dm, 5, 20))

g_dat <- data.frame(x = coords[, 1] + 5000,
                    y = coords[, 2] + 5000,
                    temp = temp[1, ],
                    precip = precip[1, ])

(temp_plot <- ggplot(g_dat, aes(x = x, y = y, fill = temp)) +
    geom_raster() +
    coord_fixed(ratio = 1) +
    scale_fill_viridis_c(name = "Temp (°C)", option = "B") +
    xlab("Easting") +
    ylab("Northing") +
    ggtitle("2-D G-space") +
    theme_bw())

(precip_plot <- ggplot(g_dat, aes(x = x, y = y, fill = precip)) +
    geom_raster() +
    coord_fixed(ratio = 1) +
    scale_fill_viridis_c(name = "Precip (mm)", option = "D") +
    xlab("Easting") +
    ylab("Northing") +
    ggtitle("2-D G-space") +
    theme_bw() +
    # theme(legend.position = "bottom") +
    NULL
)

g_plot <- temp_plot + precip_plot

ggsave(here::here("05 HSF", "figs/g-space_temp.png"), plot = temp_plot, 
       width = 400, height = 300, units = "px", dpi = 110)
ggsave(here::here("05 HSF", "figs/g-space_precip.png"), plot = precip_plot, 
       width = 400, height = 300, units = "px", dpi = 110)

# ... E-space ----
e_dat <- g_dat %>%
  mutate(temp = round(temp, 1),
         precip = round(precip, 0)) %>% 
  mutate(temp = factor(temp, levels = seq(min(temp), max(temp), by = 0.1)),
         precip = factor(precip, levels = seq(min(precip), max(precip), by = 1))) %>% 
  group_by(temp, precip, .drop = FALSE) %>% 
  tally() %>% 
  arrange(temp, precip) %>% 
  mutate(temp = as.numeric(as.character(temp)),
         precip = as.numeric(as.character(precip)))

(e_plot <- ggplot(e_dat, aes(x = temp, y = precip, fill = n)) +
    geom_raster() +
    coord_fixed(ratio = 0.12) +
    scale_fill_viridis_c(name = "Frequency") +
    xlab("Temperature (°C)") +
    ylab("Precipitation (mm)") +
    ggtitle("2-D E-space") +
    theme_bw())

ggsave(here::here("05 HSF", "figs/e-space.png"), plot = e_plot, 
       width = 450, height = 300, units = "px", dpi = 110)

# Habitat defn ----
# Defn: "a point in e-space..."
hab_pts <- data.frame(
  temp = round(runif(20, min = min(temp) + 0.5, max = max(temp) - 0.5), 1),
  precip = round(runif(20, min = min(precip) + 2, max = max(precip) - 2)))

(hab_plot <- e_plot +
    geom_point(data = hab_pts, size = 2, shape = 21, stroke = 0.75,
               fill = "#CC0000", color = "#FFFFFF"))

ggsave(here::here("05 HSF", "figs/habitat.png"), plot = hab_plot, 
       width = 600, height = 400, units = "px", dpi = 120)

# Availability ----
(temp_a <- g_dat %>% 
   ggplot(aes(x = temp)) +
   geom_density(size = 1.2) +
   xlab("Temperature (°C)") +
   ylab(expression(f[a](x))) +
   coord_cartesian(xlim = range(temp),
                   ylim = c(0, 0.625)) +
   theme_bw() + 
   NULL)

(precip_a <- g_dat %>% 
    ggplot(aes(x = precip)) +
    geom_density(size = 1.2) +
    xlab("Precipitation (mm)") +
    ylab(expression(f[a](x))) +
    coord_cartesian(xlim = range(precip),
                    ylim = c(0, 0.1)) +
    theme_bw() + 
    NULL)

(avail_fig <- (e_plot + 
                 ggtitle(NULL) +
                 theme(legend.position = "none")) + 
  (temp_a) + 
  (precip_a))

ggsave(here::here("05 HSF", "figs/avail.png"), plot = avail_fig, 
       width = 800, height = 300, units = "px", dpi = 120)

# Use ----
(temp_u <- hab_pts %>% 
   ggplot(aes(x = temp)) +
   geom_density(size = 1.2) +
   xlab("Temperature (°C)") +
   ylab(expression(f[u](x))) +
   coord_cartesian(xlim = range(temp),
                   ylim = c(0, 0.625)) +
   theme_bw() + 
   NULL)

(precip_u <- hab_pts %>% 
    ggplot(aes(x = precip)) +
    geom_density(size = 1.2) +
    xlab("Precipitation (mm)") +
    ylab(expression(f[u](x))) +
    coord_cartesian(xlim = range(precip),
                    ylim = c(0, 0.1)) +
    theme_bw() + 
    NULL)

(use_fig <- (hab_plot + 
                 ggtitle(NULL) +
                 theme(legend.position = "none")) + 
    (temp_u) + 
    (precip_u))

ggsave(here::here("05 HSF", "figs/use.png"), plot = use_fig, 
       width = 800, height = 300, units = "px", dpi = 120)

# Selection ----
(temp_w <- bind_rows("u" = g_dat,
                     "a" = hab_pts,
                     .id = "f") %>% 
   ggplot(aes(x = temp, color = f)) +
   geom_density(size = 1.2) +
   xlab("Temperature (°C)") +
   ylab(expression(f(x))) +
   coord_cartesian(xlim = range(temp),
                   ylim = c(0, 0.625)) +
   theme_bw() + 
   NULL)

(precip_w <- bind_rows("u" = g_dat,
                     "a" = hab_pts,
                     .id = "f") %>% 
    ggplot(aes(x = precip, color = f)) +
    geom_density(size = 1.2) +
    xlab("Precipitation (mm)") +
    ylab(expression(f(x))) +
    coord_cartesian(xlim = range(precip),
                    ylim = c(0, 0.1)) +
    theme_bw() + 
    NULL)

(sel_fig <- (hab_plot + 
               ggtitle(NULL) +
               theme(legend.position = "none")) + 
    (temp_w) + 
    (precip_w) +
    plot_layout(guides = "collect")) 

ggsave(here::here("05 HSF", "figs/select.png"), plot = sel_fig, 
       width = 900, height = 280, units = "px", dpi = 120)

# Simple example pie chart ----
{
  agg_png(here::here("05 HSF", "figs/simple_pie.png"), width = 600, height = 400, units = "px",
          res = 150)
  par(mar = rep(0, 4))
  pie(x = c("Grassland" = 0.80, 
            "Forest" = 0.15, 
            "Wetland" = 0.05),
      col = c("wheat", "forestgreen", "steelblue1"))
  # box()
  dev.off()
}
