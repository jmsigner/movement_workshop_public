#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#------------------Module 03 -- CTMMs------------------X
#----------------Last updated 2021-01-22---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

# Load packages ----
library(tidyverse)
library(raster)
library(amt)
library(sf)
library(crawl)

# Generate example data ----
# We will use the same BCRW as in the home ranges example to generate some
# sample data, and then we will randomly drop some data.

# Source function for BCRW
source("fun/bcrw.R")

# Once again, we'll generate data for 4 different individuals in the Bear River 
# mountains just east of Logan, UT. We'll assume our spatial coordinates are 
# in UTMs (zone 12).

centroids <- list("A01" = c("x" = 446589, "y" = 4625899),
                  "A02" = c("x" = 440284, "y" = 4618197),
                  "A03" = c("x" = 448796, "y" = 4613795),
                  "A04" = c("x" = 442367, "y" = 4627042)
)

# Set random seed
set.seed(20220126)

# Generate random walk
dat_all <- lapply(centroids, function(cent) {
  # The basic BCRW
  x <- bcrw(start_loc = cent + rnorm(n = 2, mean = 0, sd = 500),
            centroid = cent,
            n_steps = 1000,
            sl_distr = c('shape' = 1, 'scale' = 300),
            rho = 0.50, # stronger correlation (was 0.25)
            beta = 0.15) # stronger bias (was 0.1)
  
  # Assign actual dates and times (start on August 1, 2021)
  x$date <- as.POSIXct("2021-08-01 00:00:00") +
    x$t * 60 * 60
  
  # Missing GPS locations are typically autocorrelated in time.
  # Loop through the locations to decide which to drop.
  x$drop <- 0
  for (i in 2:nrow(x)) {
    x$drop[i] <- ifelse(x$drop[i-1] == 0,
                        rbinom(n = 1, size = 1, prob = 0.2),
                        rbinom(n = 1, size = 1, prob = 0.9))
  }
  
  # But keep first and last no matter what
  x$drop[which(x$t == 1)] <- 0
  x$drop[which(x$t == 1000)] <- 0
  
  # Return
  return(x)
}) %>% 
  bind_rows(.id = "ID")

# Drop those flagged locations
dat <- dat_all %>% 
  filter(drop == 0)

# How many locations for each ID?
dat %>% 
  group_by(ID) %>% 
  tally()

# Plot with good locations
ggplot(dat, aes(x = x, y = y, color = ID, group = ID)) +
  geom_point() +
  theme_bw()

# Plot dropped data for A04
a04_all <- dat_all %>% 
  filter(ID == "A04")

a04 <- dat %>% 
  filter(ID == "A04")

palette(c("#000000FF", "#00000030"))
plot(a04_all$x, a04_all$y, cex = 0.7,
     pch = 16, col = factor(a04_all$drop))
lines(a04_all$x, a04_all$y, col = 2)
lines(a04$x, a04$y, col = 1)
palette("default")

# Let's make a simple comparison between 95% MCPs for the "true" and partial
# datasets.
(mcp_area_true <- a04_all %>% 
    make_track(x, y, date, crs = 5070) %>% 
    hr_mcp() %>% 
    hr_area(units = TRUE) %>% 
    mutate(area = units::set_units(area, "km2")) %>% 
    mutate(dataset = "True"))

(mcp_area_obs <- a04 %>% 
    make_track(x, y, date, crs = 5070) %>% 
    hr_mcp() %>% 
    hr_area(units = TRUE) %>% 
    mutate(area = units::set_units(area, "km2")) %>% 
    mutate(dataset = "Observed"))

# This gives us a simple metric to compare later.

# Imputation with CRAWL ----
# Useful code to reference when using 'crawl::crwMLE()'
# https://jmlondon.github.io/crawl-workshop/crawl-practical.html#fitting-with-crawlcrwmle

# CRAWL wants our location data as a spatial object in R.
a04_sf <- a04 %>% 
  dplyr::select(ID, x, y, date) %>% 
  st_as_sf(coords = c("x", "y"), crs = 32612)

a04_all_sf <- a04_all %>% 
  dplyr::select(ID, x, y, date) %>% 
  st_as_sf(coords = c("x", "y"), crs = 32612)

# Now we fit the continuous-time CRW
mod <- crwMLE(data = a04_sf, Time.name = "date")
# Printing the object shows the summary
print(mod)

# Now we can predict a movement track
pred <- crwPredict(mod, predTime = "1 hour")

# Check the result
head(pred)

# The columns 'x' and 'y' are our original data. You can see our missing
# fixes are NAs.
# The columns 'mu.x' and 'mu.y' are the mean estimated xy-coordinates

# We have functions for converting these to 'sf' objects.
pts <- crw_as_sf(pred, ftype = "POINT")
lns <- crw_as_sf(pred, ftype = "LINESTRING")

# Plot
par(mar = rep(0, 4))
plot(lns$geometry)
plot(pts$geometry, add = T, pch = 16)

# Notice that the mean point location simply falls along a straight line between
# our observed points. However, there is uncertainty around the mean, which
# we may wish to include. For this, we want to simulate rather than predict.

# Simulate
set.seed(1)
sim <- crwSimulator(mod, predTime = "1 hour")
sim_tracks <- crwPostIS(sim, fullPost = FALSE)

sim_pts <- crw_as_sf(sim_tracks, ftype = "POINT")
sim_lns <- crw_as_sf(sim_tracks, ftype = "LINESTRING")

# Plot
par(mar = rep(0, 4))
plot(a04_all_sf$geometry, pch = 16, col = 3)
plot(sim_lns$geometry, add = T)
plot(sim_pts$geometry, add = T, pch = 16, col = factor(sim_pts$locType))
legend("topright", pch = 16, col = 1:3, 
       legend = c("Observed", "Simulated", "True"))

# This is just one realization of a stochastic process.

# Let's see how the area of our 95% MCP compares between observed, simulated,
# and true hourly locations.

(mcp_area_sim <- st_coordinates(sim_pts) %>% 
    make_track(X, Y, crs = 5070) %>% 
    hr_mcp() %>% 
    hr_area(units = TRUE) %>% 
    mutate(area = units::set_units(area, "km2")) %>% 
    mutate(dataset = "Simulated"))

rbind(mcp_area_true, mcp_area_obs, mcp_area_sim)

# You can see that our single simulation does a better job reflecting the
# true size of the 95% MCP than the observed points alone do.

# Multiple imputation ----
# If we want to reflect the uncertainty in the simulated points, we need to
# draw many realizations from the underlying stochastic process.

# The process of multiple imputation involves drawing many realizations from
# the stochastic process, calculating the statistic(s) of interest, and then
# pooling the estimates.

# Specifically, in our case, we need to:
#   1. Fit a model
#   2. For many iterations:
#       a. Simulate from the model
#       b. Calculate area of the 95% MCP
#   3. Calculate mean and quantiles of 95% MCP area

# We've already done 1, so now we just need to do 2 and 3.

# It would streamline this process if we had a function to do it all once.

# For demonstration purposes, we want the tracks, too, so our function should
# return a list with both the MCP area and the LINESTRING for each track

sim_iter <- function(model) {
  # Setup simulation
  sim <- crwSimulator(model, predTime = "1 hour")
  # Get a posterior sample
  sim_tracks <- crwPostIS(sim, fullPost = FALSE)
  # Convert to sf points
  pts <- crw_as_sf(sim_tracks, ftype = "POINT")
  # Convert to sf lines
  lns <- crw_as_sf(sim_tracks, ftype = "LINESTRING")
  # Convert to track_xy and fit 95% MCP
  mcp <- st_coordinates(pts) %>% 
    make_track(X, Y, crs = 5070) %>% 
    hr_mcp() %>% 
    hr_area(units = TRUE) %>% 
    mutate(area = units::set_units(area, "km2")) %>% 
    mutate(dataset = "Multiple Imputation")
  
  # List to return
  res <- list(line = lns,
              mcp = mcp)
  return(res)
}

# Here's how it works for one iteration:
sim_iter(mod)

# Now let's replicate that 100 times (takes about 3 minutes).
system.time({
  mult_imp <- replicate(100, sim_iter(model = mod), simplify = FALSE)
})

# Each iteration is one element of a list.
# Each one of those list elements is also a list, with two elements.
# The first is the lines, the second is the MCP area.

# Let's extract each component.
line_list <- lapply(mult_imp, function(x){
  return(x$line)
})

mcp_list <- lapply(mult_imp, function(x){
  return(x$mcp)
})

# And now combine the lists into a single data.frame
lines <- do.call(rbind, line_list)
mcps <- do.call(rbind, mcp_list)

# First, we'll summarize the MCPs. We want an estimate of the area and a
# 95% confidence interval.
(mcp_area_mi <- mcps %>% 
    group_by(level, what, dataset) %>% 
    summarize(mean = mean(area),
              lwr = quantile(area, 0.025),
              upr = quantile(area, 0.975)))

# Break that apart and compare with what we did earlier
mcp_area_mean <- mcp_area_mi %>% 
  mutate(dataset = "MI mean")%>% 
  dplyr::select(level, what, area = mean, dataset) 
mcp_area_lwr <- mcp_area_mi %>% 
  mutate(dataset = "MI lwr")%>% 
  dplyr::select(level, what, area = lwr, dataset) 
mcp_area_upr <- mcp_area_mi %>% 
  mutate(dataset = "MI upr")%>% 
  dplyr::select(level, what, area = upr, dataset) 

rbind(mcp_area_true,
      mcp_area_obs,
      mcp_area_sim,
      mcp_area_lwr,
      mcp_area_mean,
      mcp_area_upr)

# Hopefully you can see that multiple imputation is doing a much better job
# capturing the true statistic than the observed data or a single imputation
# alone.

# What does the multiple imputation look like?

# Let's plot our lines
# Plot
{
  par(mar = rep(0, 4))
  plot(lines$geometry, col = 1)
  plot(a04_all_sf$geometry, add = T, pch = 16, col = 3)
  plot(a04_sf$geometry, add = T, pch = 16, col = 2)
  plot(lines$geometry, add = T, col = "#00000033")
  legend("topright", pch = c(NA, 16, 16), lty = c(1, NA, NA), col = 1:3, 
         legend = c("Imputed", "Observed", "True"))
}


