# Figures for lecture 08 -- iSSF 2

# Load packages----
library(tidyverse)
library(amt)
library(raster)
library(ragg)
library(here)
library(circular)

# Beta as function ----

b1 <- 2
b3 <- 0.5

(bstar_fig <- data.frame(x2 = seq(0, 10, length.out = 50)) %>% 
    mutate(bs1 = b1 + b3 * x2) %>% 
    ggplot(aes(x = x2, y = bs1)) +
    geom_line(size = 1) +
    xlab(expression(x[2])) +
    ylab(expression(beta[1]^"*")) +
    ylim(c(0, 8)) +
    theme_bw())

ggsave(here("08 iSSF 2", "figs/beta_star_1.png"), plot = bstar_fig, device = agg_png,
       width = 700, height = 325, units = "px", dpi = 200)


# Resource-risk interaction ----
# ... ex 1 ----
# Risk makes foraging unacceptable
b1 <- 3
b3 <- -1

(res_risk_fig1 <- data.frame(x2 = seq(0, 6, length.out = 10)) %>% 
    mutate(bs1 = b1 + b3 * x2) %>% 
    ggplot(aes(x = x2, y = bs1)) +
    geom_hline(yintercept = 0, linetype = "dashed", 
               size = 1, color = "gray50") +
    geom_line(size = 1) +
    xlab(expression(x[2])) +
    scale_x_continuous(name = "Risk",
                       breaks = c(0, 3, 6),
                       labels = c("Low", "Med", "High")) +
    scale_y_continuous(name = expression(beta[j]^"*"),
                       breaks = c(-2, 0, 2),
                       labels = c("-", 0, "+")) +
    ggtitle("log-RSS for Resource") +
    theme_bw())

ggsave(here("08 iSSF 2", "figs/res_risk_intxn1.png"), plot = res_risk_fig1, device = agg_png,
       width = 600, height = 450, units = "px", dpi = 150)


# ... ex 2 ----
# Forage makes risk acceptable
b1 <- -3
b3 <- 1

(res_risk_fig2 <- data.frame(x2 = seq(0, 6, length.out = 10)) %>% 
    mutate(bs1 = b1 + b3 * x2) %>% 
    ggplot(aes(x = x2, y = bs1)) +
    geom_hline(yintercept = 0, linetype = "dashed", 
               size = 1, color = "gray50") +
    geom_line(size = 1) +
    xlab(expression(x[2])) +
    scale_x_continuous(name = "Resource",
                       breaks = c(0, 3, 6),
                       labels = c("Low", "Med", "High")) +
    scale_y_continuous(name = expression(beta[j]^"*"),
                       breaks = c(-2, 0, 2),
                       labels = c("-", 0, "+")) +
    ggtitle("log-RSS for Risk") +
    theme_bw())

ggsave(here("08 iSSF 2", "figs/res_risk_intxn2.png"), plot = res_risk_fig2, device = agg_png,
       width = 600, height = 450, units = "px", dpi = 150)


# Resource-condition interaction ----
# Perspective of the condition
b3 <- -4.5
b4 <- 6
b5 <- -1

(res_cond_fig <- data.frame(x1 = seq(0, 6, length.out = 100)) %>% 
    mutate(bs3 = b3 + b4 * x1 + b5 * x1^2) %>% 
    ggplot(aes(x = x1, y = bs3)) +
    geom_hline(yintercept = 0, linetype = "dashed", 
               size = 1, color = "gray50") +
    geom_line(size = 1) +
    xlab(expression(x[2])) +
    scale_x_continuous(name = "Condition",
                       breaks = c(0, 3, 6),
                       labels = c("Low", "Med", "High")) +
    scale_y_continuous(name = expression(beta[j]^"*"),
                       breaks = c(-3, 0, 3),
                       labels = c("-", 0, "+")) +
    ggtitle("log-RSS for Resource") +
    theme_bw())

ggsave(here("08 iSSF 2", "figs/res_cond_intxn.png"), plot = res_cond_fig, device = agg_png,
       width = 600, height = 350, units = "px", dpi = 150)

# Risk-condition interaction ----
# Perspective of the condition
b3 <- 2
b4 <- -6
b5 <- 1

(risk_cond_fig <- data.frame(x1 = seq(0, 6, length.out = 100)) %>% 
    mutate(bs3 = b3 + b4 * x1 + b5 * x1^2) %>% 
    ggplot(aes(x = x1, y = bs3)) +
    geom_hline(yintercept = 0, linetype = "dashed", 
               size = 1, color = "gray50") +
    geom_line(size = 1) +
    xlab(expression(x[2])) +
    scale_x_continuous(name = "Condition",
                       breaks = c(0, 3, 6),
                       labels = c("Low", "Med", "High")) +
    scale_y_continuous(name = expression(beta[j]^"*"),
                       breaks = c(-4, -2, 0, 2),
                       labels = c("-", "", 0, "+")) +
    ggtitle("log-RSS for Risk") +
    theme_bw())

ggsave(here("08 iSSF 2", "figs/risk_cond_intxn.png"), plot = risk_cond_fig, device = agg_png,
       width = 600, height = 350, units = "px", dpi = 150)

# Resource-resource interaction ----
# ... ex 1 ----
# Resources are antagonistic
b1 <- 3.5
b3 <- -0.5

(res_antag_fig <- data.frame(x2 = seq(0, 6, length.out = 10)) %>% 
    mutate(bs1 = b1 + b3 * x2) %>% 
    ggplot(aes(x = x2, y = bs1)) +
    geom_hline(yintercept = 0, linetype = "dashed", 
               size = 1, color = "gray50") +
    geom_line(size = 1) +
    xlab(expression(x[2])) +
    scale_x_continuous(name = "Resource 2",
                       breaks = c(0, 3, 6),
                       labels = c("Low", "Med", "High")) +
    scale_y_continuous(name = expression(beta[j]^"*"),
                       breaks = c(0, 3),
                       labels = c(0, "+")) +
    ggtitle("log-RSS for Resource 1") +
    theme_bw())

ggsave(here("08 iSSF 2", "figs/res_antag.png"), plot = res_antag_fig, device = agg_png,
       width = 600, height = 450, units = "px", dpi = 150)


# ... ex 2 ----
# Resources are complementary
b1 <- 0.75
b3 <- 0.5

(res_comp_fig <- data.frame(x2 = seq(0, 6, length.out = 10)) %>% 
    mutate(bs1 = b1 + b3 * x2) %>% 
    ggplot(aes(x = x2, y = bs1)) +
    geom_hline(yintercept = 0, linetype = "dashed", 
               size = 1, color = "gray50") +
    geom_line(size = 1) +
    xlab(expression(x[2])) +
    scale_x_continuous(name = "Resource 2",
                       breaks = c(0, 3, 6),
                       labels = c("Low", "Med", "High")) +
    scale_y_continuous(name = expression(beta[j]^"*"),
                       breaks = c(0, 3),
                       labels = c(0, "+")) +
    ggtitle("log-RSS for Resource 1") +
    theme_bw())

ggsave(here("08 iSSF 2", "figs/res_comp.png"), plot = res_comp_fig, device = agg_png,
       width = 600, height = 450, units = "px", dpi = 150)

# Tentative step-length dist ----
shp0 <- 3
scl0 <- 100

(tent_sl_fig <- data.frame(sl = seq(1, 2000, length.out = 100)) %>% 
    mutate(y = dgamma(sl, shape = shp0, scale = scl0)) %>% 
    ggplot(aes(x = sl, y = y)) +
    geom_line(size = 1) +
    xlab("Step Length (m)") +
    ylab("Probability Density") +
    theme_bw())


ggsave(here("08 iSSF 2", "figs/tent_sl.png"), plot = tent_sl_fig, 
       device = agg_png, width = 600, height = 400, units = "px", dpi = 150)

# Updated sl ----
shp_update <- function(x, b3, b5) {
  b_log_sl_ <- b3 + b5 * x
  return(shp0 + b_log_sl_)
}

scl_update <- function(x, b2, b4) {
  b_sl_ = b2 + b4 * x
  return(1/((1/scl0)) - b_sl_)
}

# ... continuous x ----
x_cont <- seq(-2, 2, length.out = 10)
# Shape
b3_cont <- 3
b5_cont <- -2

# Scale
b2_cont <- 0
b4_cont <- -0.1

(updated_sl_cont <- expand.grid(sl = seq(1, 2000, length.out = 100),
                                x = x_cont) %>% 
    mutate(shp = shp_update(x, b3_cont, b5_cont),
           scl = scl_update(x, b2_cont, b4_cont),
           y = dgamma(sl, shape = shp, scale = scl)) %>% 
    ggplot(aes(x = sl, y = y, color = x, group = x)) +
    geom_line(size = 1) +
    xlab("Step length (m)") +
    ylab("Probability Density") +
    scale_color_gradient(name = "Vegetation\nDensity (SD)", 
                         low = "wheat", high = "forestgreen") +
    theme_bw())

ggsave(here("08 iSSF 2", "figs/updated_sl_cont.png"), plot = updated_sl_cont, 
       device = agg_png, width = 700, height = 350, units = "px", dpi = 150)


# ... categorical x ----
# Shape
shp_grass <- 6
shp_for <- 3
shp_wet <- 1.5

# Scale
scl_grass <- 100
scl_for <- 100
scl_wet <- 100

(updated_sl_cat <- expand.grid(sl = seq(0.1, 2000, length.out = 100),
                               x = c("Grassland", "Forest", "Wetland")) %>% 
    mutate(
      shp = case_when(
        x == "Grassland" ~ shp_grass,
        x == "Forest" ~ shp_for,
        x == "Wetland" ~ shp_wet
      ),
      scl = case_when(
        x == "Grassland" ~ scl_grass,
        x == "Forest" ~ scl_for,
        x == "Wetland" ~ scl_wet
      ),
      y = dgamma(sl, shape = shp, scale = scl)) %>% 
    ggplot(aes(x = sl, y = y, color = x, group = x)) +
    geom_line(size = 1) +
    xlab("Step length (m)") +
    ylab("Probability Density") +
    scale_color_manual(name = "Habitat\nType", 
                       breaks = c("Grassland", "Forest", "Wetland"),
                       values = c("wheat", "forestgreen", "steelblue1")) +
    theme_bw())

ggsave(here("08 iSSF 2", "figs/updated_sl_cat.png"), plot = updated_sl_cat, 
       device = agg_png, width = 700, height = 350, units = "px", dpi = 150)

# Updated von Mises ----

(updated_ta_cont <- expand.grid(x = seq(-pi, pi, length.out = 100),
                                k = seq(-1, 1, length.out = 10)) %>% 
   mutate(
     mu = case_when(
       k < 0 ~ pi,
       k > 0 ~ 0,
       TRUE ~ 0
     ),
     k_abs = abs(k),
     forage = -2 * k) %>% 
   rowwise() %>% 
   mutate(y = dvonmises(x = x, mu = mu, kappa = k_abs)) %>% 
   ggplot(aes(x = x, y = y, color = forage, group = forage)) +
   geom_line(size = 1) +
   xlab("Turn Angle (radians)") +
   ylab("Probability Density") +
   scale_color_gradient(name = "Forage\nQuality (SD)", 
                        low = "wheat", high = "forestgreen") +
   scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi), 
                      labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
   theme_bw())

ggsave(here("08 iSSF 2", "figs/updated_ta_cont.png"), plot = updated_ta_cont, 
       device = agg_png, width = 700, height = 300, units = "px", dpi = 150)

# Turn angle and speed ----

(updated_ta_sl <- expand.grid(x = seq(-pi, pi, length.out = 100),
                                k = seq(-1, 1, length.out = 10)) %>% 
   mutate(
     mu = case_when(
       k < 0 ~ pi,
       k > 0 ~ 0,
       TRUE ~ 0
     ),
     k_abs = abs(k),
     sl = (k + 1) * 500) %>% 
   rowwise() %>% 
   mutate(y = dvonmises(x = x, mu = mu, kappa = k_abs)) %>% 
   ggplot(aes(x = x, y = y, color = sl, group = sl)) +
   geom_line(size = 1) +
   xlab("Turn Angle (radians)") +
   ylab("Probability Density") +
   scale_color_gradient2(name = "Mean Step\nLength (m)",
                        midpoint = 500) +
   scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi), 
                      labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
   theme_bw())

ggsave(here("08 iSSF 2", "figs/updated_ta_sl.png"), plot = updated_ta_sl, 
       device = agg_png, width = 700, height = 400, units = "px", dpi = 150)
