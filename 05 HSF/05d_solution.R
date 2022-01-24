#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#-------------------Module 05 -- HSFs------------------X
#----------------Last updated 2021-01-22---------------X
#-------------------Exercise Solution------------------X
#######################################################X

# Using amt_fisher data

# Load packages ----
library(tidyverse)
library(amt)
library(lubridate)
library(raster)

# 1. Load data ----
# Location data as 'track_xyt'
dat <- amt_fisher %>% 
  # Subset to one individual
  filter(name == "Ricky T")

# Habitat data as list
hab <- amt_fisher_covar

# ... subset data ----
# Subset to 4 locations per day
# One approach is to find the time closest to the hour
# First, only keep hours 0, 6, 12, and 18

# Keep only the hour of the time by setting seconds and minutes to 0
dat$h <- dat$t_
minute(dat$h) <- second(dat$h) <- 0

# Subset
sub <- dat %>% 
  # Keep only hours 0, 6, 12, and 18
  filter(hour(t_) %in% c(0, 6, 12, 18)) %>% 
  # Get difference
  mutate(diff = difftime(t_, h, units = "mins")) %>% 
  # Group by hour
  group_by(h) %>% 
  # Pick the row which has the smallest diff
  filter(diff == min(diff)) %>% 
  # Get rid of h and diff
  ungroup() %>% 
  dplyr::select(-h, -diff)

# Check sampling rate
summarize_sampling_rate(sub)

# 2. Consider habitat variables ----
# Plot them
lapply(hab, function(x){
  plot(x, main = names(x))
})

# Land use is meant to be categorical, but is stored as integers
# Make sure we convert that to factor
# Some land use categories might be resources and some might be risks

# Elevation is a condition
# Make sure we use quadratic term
# Would be interesting to know what elevation is most selected

# Pop den is human population density
# This is probably a risk to fishers

# 3. Fit HSF ----
# Third-order = within home range
# Let's use a 95% KDE as our home range
kde <- sub %>% 
  hr_kde(levels = 0.95)
plot(kde)

# Where does it fall?
plot(hab[[2]])
plot(hr_isopleths(kde)$geometry, add = T)

# Sample random points within HR
rand <- random_points(kde, n = nrow(sub) * 100)
plot(rand)

# Create dataset for model
hsf_dat <- sub %>%
  # Prep observed locations
  mutate(case_ = TRUE) %>% 
  dplyr::select(case_, x_, y_) %>% 
  # Combine with random locations
  bind_rows(rand) %>% 
  # Attach environmental covariates
  extract_covariates(hab$landuse) %>% 
  extract_covariates(hab$elevation) %>% 
  extract_covariates(hab$popden) %>% 
  # Convert landuse to factor
  mutate(landuse = factor(landuse)) %>% 
  # Create weights for observed and available locations
  mutate(weight = case_when(
    case_ ~ 1,
    !case_ ~ 5000))

head(hsf_dat)

# Fit model 
hsf <- glm(case_ ~ landuse + elevation + I(elevation^2) + popden,
           data = hsf_dat, family = binomial, weights = weight)

# Check summary
summary(hsf)

# 4. RSS Figures ----

# ... land use ----
# Levels of land use
lu_levs <- levels(hsf_dat$landuse)

# x1
x1_lu <- data.frame(landuse = factor(lu_levs, levels = lu_levs),
                    elevation = mean(hsf_dat$elevation),
                    popden = mean(hsf_dat$popden))
# x2
x2_lu <- data.frame(landuse = factor("30", levels = lu_levs),
                   elevation = mean(hsf_dat$elevation),
                   popden = mean(hsf_dat$popden))

# log-RSS
log_rss_lu <- log_rss(hsf, x1_lu, x2_lu, ci = "se")

# Plot
log_rss_lu$df %>% 
  ggplot(aes(x = landuse_x1, y = log_rss, ymin = lwr, ymax = upr)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(width = 0.1) +
  geom_point(size = 2) +
  xlab("Land Use Code") +
  ylab("log-RSS vs. category 30") +
  coord_cartesian(ylim = c(-13, 5)) +
  theme_bw()

# Interpretation:
# Compared to category 30, Ricky T:
#   - Selects category 50 equally
#   - Avoids category 70
#   - Might avoid category 100, but some overlap with 0
#   - Prefers category 110
#   - Selects category 210 equally
# As for categories 120 and 140, they come with a huge amount of uncertainty,
# so it's hard for us to draw inference on them.

# ... elevation ----
# x1
x1_elev <- data.frame(landuse = factor("30", levels = lu_levs),
                    elevation = seq(60, 140, length.out = 100),
                    popden = mean(hsf_dat$popden))
# x2
x2_elev <- data.frame(landuse = factor("30", levels = lu_levs),
                    elevation = mean(hsf_dat$elevation),
                    popden = mean(hsf_dat$popden))

# log-RSS
log_rss_elev <- log_rss(hsf, x1_elev, x2_elev, ci = "se")

# Plot
log_rss_elev$df %>% 
  ggplot(aes(x = elevation_x1, y = log_rss, ymin = lwr, ymax = upr)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(linetype = "dashed", fill = "gray80", 
              color = "black", alpha = 0.3) +
  geom_line(size = 1) +
  xlab("Elevation") +
  ylab("log-RSS vs. mean elevation") +
  # coord_cartesian(ylim = c(-13, 5)) +
  theme_bw()

# Intepretation:
# Ricky T does not select for elevation. We hypothesized it was a condition,
# so we fit a parabola, which we expected to be concave down. It is actually
# concave-up, but the confidence envelope overlaps 0 for almost the entire
# range of elevation values. This is supported by the large p-values for 
# elevation and elevation^2 in the model summary.

# ... popden ----
# x1
x1_popden <- data.frame(landuse = factor("30", levels = lu_levs),
                      elevation = mean(hsf_dat$elevation),
                      popden = seq(0, 1500, length.out = 100))
# x2
x2_popden <- data.frame(landuse = factor("30", levels = lu_levs),
                      elevation = mean(hsf_dat$elevation),
                      popden = mean(hsf_dat$popden))

# log-RSS
log_rss_popden <- log_rss(hsf, x1_popden, x2_popden, ci = "se")

# Plot
log_rss_popden$df %>% 
  ggplot(aes(x = popden_x1, y = log_rss, ymin = lwr, ymax = upr)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(linetype = "dashed", fill = "gray80", 
              color = "black", alpha = 0.3) +
  geom_line(size = 1) +
  xlab("Population Density") +
  ylab("log-RSS vs. mean population density") +
  # coord_cartesian(ylim = c(-13, 5)) +
  theme_bw()

# Interpretation:
# Ricky T avoids human population density. We hypothesized that this habitat
# variable was a risk, and our results support that hypothesis.

# An additional metric of use might be to try some scenarios to see how
# much more Ricky T likes some densities over others.

# How much more time does Ricky T spend at popden = 100 than popden = 500?
x1a <- data.frame(landuse = factor("30", levels = lu_levs),
                  elevation = mean(hsf_dat$elevation),
                  popden = 100)
x2a <- data.frame(landuse = factor("30", levels = lu_levs),
                  elevation = mean(hsf_dat$elevation),
                  popden = 500)

log_rss(hsf, x1a, x2a)$df %>% 
  mutate(rss = exp(log_rss))

# We are about 2x as likely to find Ricky T at popden 100 than popden 500.
# I.e., we expect 2x the density of telemetry locations at x1a vs x2a.

# How much more time does Ricky T spend at popden = 100 than popden = 1500?
x1b <- data.frame(landuse = factor("30", levels = lu_levs),
                  elevation = mean(hsf_dat$elevation),
                  popden = 100)
x2b <- data.frame(landuse = factor("30", levels = lu_levs),
                  elevation = mean(hsf_dat$elevation),
                  popden = 1500)

log_rss(hsf, x1b, x2b)$df %>% 
  mutate(rss = exp(log_rss))

# We are about 10x as likely to find Ricky T at popden 100 than popden 1500.
# I.e., we expect 10x the density of telemetry locations at x1b vs x2b.
