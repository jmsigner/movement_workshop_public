#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#----------------Module 08 -- iSSF pt 2----------------X
#----------------Last updated 2021-01-21---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

# Load packages ----
library(tidyverse)
library(raster)
library(amt)
library(lubridate)
library(circular)

# Generate example data ----
# We need to generate some example telemetry data to fit our iSSF.
# Generating points from an iSSF is trickier than generating points
# from an HSF because of (1) the dependence on the start location  
# (i.e., it is a time series) and (2) the addition of the movement model.

# For this walkthrough, we're interested in a fairly complex model, i.e.,
# we want to see how interactions work with both habitat and movement
# variables. 

# We'll stick with one individual for now.

# ... habitat variables ----
# We'll use the same habitat layers we generated for module 5.
hab <- stack("05 HSF/geo/habitat.tif")
names(hab) <- c("forage", "temp", "predator", "cover")
plot(hab)

# ... movement-free habitat kernel ----
# Now we have to choose betas that will give us our movement-free habitat
# selection kernel.

# We'll start with the same values for beta from module 5. We still interpret
# each beta as the log-RSS for a one-unit change in the covariate. But 
# now we need to recognize that the ratio of densities is the density
# of *steps* ending in that habitat, GIVEN THE START LOCATION.

beta_forage = log(5)/500
beta_pred = log(0.25)/5
beta_temp2 = -1 * log(2)/36
beta_temp = beta_temp2 * -26

# Now let's add a bit of complexity. Let's say that our animal is nocturnal,
# and they spend the night foraging in grasslands and the day resting in
# forests.

# Let's stick with grassland as the reference category as it was in module 05.
# Specifically, we'll use grassland during daytime as the reference.
# We'll pick some betas such that:
#   Daytime:
#     Forest 10x grassland
#     Wetland 1/3 x grassland
#   Nighttime:
#     Forest 1/3x grassland 
#     Wetland 1/5x grassland 

# To be clear, we want the time of day at the end of the step to affect the
# habitat selected at the end of the step.

beta_forest_day = log(10)
beta_wetland_day = log(1/3)
beta_forest_night = log(1/3)
beta_wetland_night = log(1/5)

# Since our landscape is small we can calculate the movement-free habitat
# selection kernel for the entire landscape. It only depends on the time of
# day, so we can simply calculate one for day and one for night.
hab_vals <- as.data.frame(hab) %>% 
  # Convert cover to factor
  mutate(cover = factor(cover,
                        levels = 1:3,
                        labels = c("grassland", "forest", "wetland"))) %>% 
  # Calculate g(x) for day and night
  mutate(g_day =
           # forage
           beta_forage * forage +
           # two terms for temperature
           beta_temp * temp +
           beta_temp2 * temp^2 +
           # predator density
           beta_pred * predator +
           # landcover
           beta_forest_day * (cover == "forest") +
           beta_wetland_day * (cover == "wetland"),
         g_night =
           # forage
           beta_forage * forage +
           # two terms for temperature
           beta_temp * temp +
           beta_temp2 * temp^2 +
           # predator density
           beta_pred * predator +
           # landcover
           beta_forest_night * (cover == "forest") +
           beta_wetland_night * (cover == "wetland")) %>% 
  # Exponentiate to get w(x)
  mutate(w_day = exp(g_day),
         w_night = exp(g_night)) %>% 
  # Normalize so they sum to 1
  mutate(w_prime_day = w_day/sum(w_day),
         w_prime_night = w_night/sum(w_night))

# Place in raster
# Use hab[[1]] as template
hab_kern_day <- hab_kern_night <- hab[[1]]
# Insert values
values(hab_kern_day) <- hab_vals$w_prime_day
values(hab_kern_night) <- hab_vals$w_prime_night
# Plot
plot(stack(hab_kern_day, hab_kern_night))

# ... selection-free movement kernel ----

# Now we have to specify the movement part of the model. 

# Let's assume our steps come from a gamma distribution and our turn
# angles come from a von Mises distribution.

# As we said above, our animal is nocturnal, so we want a step-length 
# distribution with short steps during the day and long steps at night.
# I.e., steps that start during day are short and steps that start during
# night are long.

# Also note that our landscape is only 2 km x 2 km, so our steps should be
# quite short to stay within the landscape.

# Let's choose these shape and scale parameters.
shp_day <- 2
scl_day <- 25

shp_night <- 4
scl_night <- 75

# Note that the mean of the gamma is given by shape * scale
# Day
shp_day * scl_day
# Night
shp_night * scl_night

# Plot the distributions
expand.grid(sl = seq(0.1, 1000, length.out = 100),
            time = c("day", "night")) %>% 
  mutate(
    shp = case_when(
      time == "day" ~ shp_day,
      time == "night" ~ shp_night),
    scl = case_when(
      time == "day" ~ scl_day,
      time == "night" ~ scl_night
    ),
    y = dgamma(sl, shape = shp, scale = scl)) %>% 
  ggplot(aes(x = sl, y = y, color = time)) +
  geom_line(size = 1) +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  theme_bw()

# We also need to specify our turn angle distribution.

# Let's assume our turn angles are correlated with our step lengths,
# such that short steps are more concentrated around +/- pi and long
# steps are more concentrated around 0.

# Let's say that the concentration parameter of the von Mises is a linear
# function of step length:
# k = -0.1 + 0.001 * sl

# We can write a function for that.
calc_k <- function(sl) {
  return(-0.5 + 0.0015 * sl)
}

calc_k(seq(0, 1000, length.out = 10))

# Recall that the concentration parameter of the von Mises cannot be negative.
# However, a negative kappa with mu = 0 can be written as abs(kappa) with 
# mu = pi (or - pi).

# Plot the distributions
expand.grid(ta = seq(-pi, pi, length.out = 100),
            sl = seq(0.1, 1000, length.out = 10)) %>% 
  mutate(
    k = calc_k(sl),
    mu = case_when(
      k < 0 ~ pi,
      TRUE ~ 0
    ),
    k_abs = abs(k)) %>% 
  rowwise() %>% 
  mutate(y = circular::dvonmises(ta, mu = mu, kappa = k_abs)) %>% 
  ggplot(aes(x = ta, y = y, color = sl, group = sl)) +
  geom_line(size = 1) +
  xlab("Turn Angle (radians)") +
  ylab("Probability Density") +
  scale_color_viridis_c() +
  theme_bw()

# The probability of stepping into any cell depends on the starting point,
# so we can't pre-calculate like we did with the habitat kernel.

# ... times ----

# Now we need to define what time step our movement parameters correspond
# to. Let's say we have a GPS fix every hour and we want one month of data. 
# We'll straddle the spring equinox to have approximately the same number of 
# day and night hours to work with.

dat <- data.frame(id = "A01",
                  time = seq(ymd_hms("2021-03-05 0:00:00", tz = "US/Mountain"),
                             ymd_hms("2021-04-05 00:00:00", tz = "US/Mountain"),
                             by = "1 hour"),
                  x = NA,
                  y = NA)

# We'll start our animal right in the middle of our map
dat$x <- mean(c(extent(hab)@xmin, extent(hab)@xmax))
dat$y <- mean(c(extent(hab)@ymin, extent(hab)@ymax))

# Check
plot(hab[[1]])
points(dat$x, dat$y, pch = 16, col = "red")

# We also want to assign the time of day to our locations. 
# `amt` can help with that.
?time_of_day

dat <- dat %>% 
  make_track(x, y, time, id = id, crs = 32612) %>% 
  time_of_day() %>% 
  # Back to a regular data.frame
  as.data.frame() %>% 
  # Rearrange/rename columns
  dplyr::select(id, x1 = x_, y1 = y_, t1 = t_, tod_start = tod_) %>% 
  # Convert to steps
  mutate(x2 = lead(x1),
         y2 = lead(y1),
         t2 = lead(t1),
         tod_end = lead(tod_start)) %>% 
  filter(!is.na(t2))

# Let's define our first step's endpoint as being 50m directly
# north to get us started moving. This is also the second step's start point.
dat$y2[1] <- dat$y1[2] <- dat$y1[1] + 50

# Lastly, let's add the absolute angle of the first step
dat$abs_angle <- NA
dat$abs_angle[1] <- 0 # directly north

# ... simulate movement ----
# We're ready to simulate!

# One useful thing to have pre-calculated are the xy coordinates of each
# raster cell.
coords <- xyFromCell(hab, 1:ncell(hab))

# We'll also want our jitter function from module 05.

# Function to jitter data
jitter <- function(x, y, min = -25, max = 25) {
  res <- data.frame(x = x + runif(1, min, max),
                    y = y + runif(1, min, max))
  return(res)
}

# This code will take several minutes to run. I have tried to write it
# for maximum ease of understanding, not computational speed. Rather
# than actually running it for all rows of the data.frame, I'll demonstrate
# how it works, but then we'll load in the results that I already generated
# previously.

# We already have the first step. We need to simulate the rest.
set.seed(20220127)

for (i in 2:nrow(dat)) {
  # Report status
  cat("\nStep", i, "of", nrow(dat))
  
  ## Is our step day or night?
  tod_st <- dat$tod_start[i]
  tod_end <- dat$tod_end[i]
  
  ## Calculate selection-free movement kernel
  # Start point
  start <- cbind(dat$x1[i], dat$y1[i])
  # Distances along x and y to every cell
  dx <- coords[, 1] - start[, 1]
  dy <- coords[, 2] - start[, 2]
  # Distance to every cell
  dists <- sqrt(dx^2 + dy^2)
  # Absolute angle to every cell
  abs <- (pi/2 - atan2(dy, dx)) %% (2*pi)
  # Relative angle difference
  rel_diff <- (abs - dat$abs_angle[i-1])
  # Relative angle
  rel_angle <- ifelse(rel_diff > pi, rel_diff - 2*pi, rel_diff)
  # Likelihood of step length
  sl_like <- dgamma(dists, 
                    shape = get(paste0("shp_", tod_st)),
                    scale = get(paste0("scl_", tod_st)))
  # Concentration of von Mises
  k <- calc_k(sl = dists)
  # Mean of von Mises
  mu <- ifelse(k < 0, pi, 0)
  # Likelihood of turn angle (not vectorized over mu -- need loop)
  ta_like <- rep(NA, length(mu))
  for (j in 1:length(mu)) {
    suppressWarnings({
      ta_like[j] <- dvonmises(rel_angle[j],
                              mu = mu[j],
                              kappa = abs(k[j]))
      
    })
  }
  
  # Calculate kernel values
  move_kern_vals <- sl_like * ta_like
  # Normalize (sum to 1)
  move_kern_vals <- move_kern_vals/sum(move_kern_vals)
  
  # # If you want to plot this
  # move_kern <- hab[[1]]
  # values(move_kern) <- move_kern_vals
  # plot(move_kern, main = "Movement Kernel")
  
  ## Movement-free habitat kernel
  # We pre-computed these for day and night
  hab_kern <- get(paste0("hab_kern_", tod_end))
  hab_kern_vals <- values(hab_kern)
  
  # If you want to plot
  # plot(hab_kern, main = "Habitat Kernel")
  
  ## Combine
  step_kern_vals <- move_kern_vals * hab_kern_vals
  # Normalize
  step_kern_vals <- step_kern_vals/sum(step_kern_vals)
  
  # # If you want to visualize
  # step_kern <- hab[[1]]
  # values(step_kern) <- step_kern_vals
  # plot(step_kern, main = "Habitat x Movement Kernel")
  
  # Randomly select cell to move into based on the probabilities
  next_cell <- sample(x = 1:ncell(hab),
                      size = 1,
                      prob = step_kern_vals)
  
  # Get cell coordinates
  next_cell_coords <- xyFromCell(hab, next_cell)
  
  # If you want to plot
  # points(next_cell_coords[,"x"], next_cell_coords[,"y"], pch = 16)
  
  # Jitter
  next_coords <- jitter(next_cell_coords[, 1], next_cell_coords[, 2])
  
  # Insert into data.frame
  if (i != nrow(dat)) {
    dat$x2[i] <- dat$x1[i+1] <- next_coords[, 1]
    dat$y2[i] <- dat$y1[i+1] <- next_coords[, 2]
  } else {
    dat$x2[i] <- next_coords[, 1]
    dat$y2[i] <- next_coords[, 2]
  }
  
  # Calculate absolute angle
  dx <- dat$x2[i] - dat$x1[i]
  dy <- dat$y2[i] - dat$y1[i]
  dat$abs_angle[i] = (pi/2 - atan2(dy, dx)) %% (2*pi)
  
  # If you want to check
  dat[i, ]
  
}

# Save results
# saveRDS(dat, "08 iSSF 2/sim.rds")

# Done!

# Load results
dat <- readRDS("08 iSSF 2/sim.rds")

# Now that we've simulated our data, let's see what our trajectory looks like.
traj <- dat %>% 
  filter(!is.na(t2)) %>% 
  select(id, x = x2, y = y2, time = t2) %>% 
  make_track(x, y, time, id = id, crs = 32612)

plot(hab[[1]])
points(traj)
lines(traj)

# Model fitting ----
# Now we can use our simulated data and fit a model.
head(traj)

# ... format data ----
# Random seed since we're generating random steps
set.seed(20220127 * 2)

issf_dat <- traj %>% 
  steps() %>% 
  # Default is to use gamma and von Mises
  random_steps(n_control = 20) %>% 
  time_of_day(where = "both") %>% 
  extract_covariates(hab, where = "end") %>% 
  # Convert cover to factor
  mutate(cover = factor(cover,
                        levels = 1:3,
                        labels = c("grassland", "forest", "wetland")),
         tod_start_ = factor(tod_start_),
         tod_end_ = factor(tod_end_)) %>% 
  mutate(log_sl_ = log(sl_),
         cos_ta_ = cos(ta_)) %>% 
  filter(!is.na(ta_))

# Take a look at what we've got.
head(as.data.frame(issf_dat))

# Notice that it has attributes, including the CRS, the tentative movement
# parameters, and the number of control steps.
str(issf_dat, 2)

# We have helper functions to extract the tenative movement parameters
sl_distr(issf_dat)
ta_distr(issf_dat)

# For later, we might to be able to access just the observed steps, so we can
# store them separately now.
obs <- issf_dat %>% 
  filter(case_)

# ... fit iSSF ----
# Fit the model

issf <- issf_dat %>% 
  # Make your own dummy variables -- R makes too many levels with
  # interactions between two categorical variables for clogit models.
  # Remember grassland is our reference level
  mutate(forest_day = as.numeric(cover == "forest" & tod_start_ == "day"),
         wetland_day = as.numeric(cover == "wetland" & tod_start_ == "day"),
         forest_night = as.numeric(cover == "forest" & tod_start_ == "night"),
         wetland_night = as.numeric(cover == "wetland" & tod_start_ == "night")) %>%
  # Fit the model
  fit_issf(case_ ~ 
             # Habitat
             forage + temp + I(temp^2) + predator + 
             # All the cover terms
             forest_day + wetland_day + 
             forest_night + wetland_night +
             # Movement
             tod_start_ : log_sl_ + tod_start_ : sl_ +
             # Correlation between sl_ and ta_
             cos_ta_ + cos_ta_ : sl_ +
             # Strata (steps)
             strata(step_id_), model = TRUE)

# Model summary
summary(issf)

# Let's look at the structure of the 'issf' object.
str(issf, 1)

# Notice that it is a list with 4 elements at the top level.
#   - $model: the actual fitted model
#   - $sl_: the tentative step-length distribution
#   - $ta_: the tentative turn-angle distribution
#   - $more: (currently empty) a placeholder for additional information

# Let's make a figure to see how well we did estimating the parameters. Note
# that more generic functions are available to us if we call the model object
# directly rather than the "fit_clogit" object.
b <- coef(issf)
ci <- confint(issf$model)

res <- as.data.frame(cbind(b, ci))
res$truth <- c(beta_forage,
               beta_temp,
               beta_temp2,
               beta_pred,
               beta_forest_day,
               beta_wetland_day,
               beta_forest_night,
               beta_wetland_night,
               rep(NA, 6))
names(res) <- c("est", "lwr", "upr", "truth")

# Plot
(beta_fig <- res %>% 
    # Keep only habitat betas
    filter(!is.na(truth)) %>% 
    mutate(name = row.names(.)) %>% 
    ggplot(aes(x = name, y = est)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 1.25) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  width = 0.2) +
    geom_point(aes(y = truth), color = "red") +
    xlab(expression(beta)) +
    ylab(NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))

# Zoom in to see forage and temp^2
beta_fig +
  coord_cartesian(ylim = c(-0.025, 0.005))

# We did fairly well. The true values for the forest day/night betas are a 
# little more extreme than what we estimated, but we did a good job
# considering the complexity of our model and our sample size.

# ... update movement distributions ----
# We did well with the movement-free habitat kernel. How about with the
# selection-free movement kernel?

# We need to use the movement betas from the fitted model to update the 
# tentative step-length and turn-angle distributions to the estimated 
# selection-free distributions.

# Recall that the formulas are available in Appendix C of Fieberg et al. 2021.
# https://conservancy.umn.edu/bitstream/handle/11299/218272/AppC_iSSA_movement.html

# We also have functions in `amt` to make this a bit easier for you.
?update_gamma

# Note that you cannot simply use 'update_sl_distr()' or 'update_ta_distr()' if
# you have interactions with your movement parameters. You need to pass the
# fully updated betas (the beta-stars from the lecture) to the correct updating
# function.

# In our model, the beta for sl is involved in an interaction with tod and 
# cos(ta). 

# That means we *can* think of the beta-star as a function of tod and cos(ta).
# But in our model formulation, it makes more sense to think of sl and log(sl)
# as fixed, and cos(ta) as a function of sl. Let's try.

# Update step-length distribution

# For tod = "day"

b_sl_day <- b[["tod_start_day:sl_"]]
b_log_sl_day <- b[["tod_start_day:log_sl_"]]

sl_distr_day <- update_gamma(sl_distr(issf),
                             beta_sl = b_sl_day,
                             beta_log_sl = b_log_sl_day)

# And what if tod = "night"?

b_sl_night <- b[["tod_start_night:sl_"]]
b_log_sl_night <- b[["tod_start_night:log_sl_"]]

sl_distr_night <- update_gamma(sl_distr(issf),
                               beta_sl = b_sl_night,
                               beta_log_sl = b_log_sl_night)

# How did we do?
data.frame(tod = c("day", "day", "night", "night"),
           parm = c("shape", "scale", "shape", "scale"),
           est = c(sl_distr_day$params$shape,
                   sl_distr_day$params$scale,
                   sl_distr_night$params$shape,
                   sl_distr_night$params$scale),
           truth = c(shp_day,
                     scl_day,
                     shp_night,
                     scl_night))

# We did very well!

# Now how about the turn-angle distribution? It depends on the step length.
# Let's use an average value for day steps and an average value for night steps.
# Remember, the mean of the gamma is shape * scale
avg_sl_day <- sl_distr_day$params$shape * sl_distr_day$params$scale
avg_sl_night <- sl_distr_night$params$shape * sl_distr_night$params$scale

cos_ta_day <- b[["cos_ta_"]] + b[["cos_ta_:sl_"]] * avg_sl_day

ta_distr_day <- update_vonmises(ta_distr(issf), 
                                beta_cos_ta = cos_ta_day)

cos_ta_night <- b[["cos_ta_"]] + b[["cos_ta_:sl_"]] * avg_sl_night

ta_distr_night <- update_vonmises(ta_distr(issf), 
                                beta_cos_ta = cos_ta_night)

# Recall from our simulations that we defined kappa as a function of sl.
# Let's use that to compare our results to truth.
data.frame(tod = c("day", "night"),
           sl = c(avg_sl_day, avg_sl_night),
           est = c(ta_distr_day$params$kappa, 
                   ta_distr_night$params$kappa)) %>% 
  mutate(truth = calc_k(sl))

# We didn't do such a great job estimating our turn angle distribution. We're
# consistently underestimating the concentration parameter (maybe because
# our landscape was too small! Our animal kept bumping into the boundary!)

# Figures ----
# Given we already know so much about the biology of our system from simulating
# our data, we can hopefully decide on some interesting figures to express it.

# What are some key aspects of our animal's biology that we want to capture?
#   - Different landcover selection for day vs night
#   - Different step lengths for day vs night
#   - Positive selection for forage
#   - Negative selection for predation risk
#   - Selection of an intermediate temperature

# We'll demonstrate the first two.

# ... landcover ----
# Let's use RSS to express selection for landcover types for day vs night
# We can use 'amt::log_rss()' just like we did for HSFs.

# Day
x1_lc_day <- data.frame(forage = mean(obs$forage),
                        temp = mean(obs$temp),
                        predator = mean(obs$predator),
                        cover = factor(c("grassland", "forest", "wetland"),
                                       levels = c("grassland", "forest", "wetland")),
                        tod_start_ = factor("day", levels = c("day", "night")),
                        tod_end_ = factor("day", levels = c("day", "night")),
                        sl_ = 100,
                        log_sl_ = log(100),
                        cos_ta_ = 1) %>% 
  # Make your own dummy variables
  # Remember grassland is our reference level
  mutate(forest_day = as.numeric(cover == "forest" & tod_start_ == "day"),
         wetland_day = as.numeric(cover == "wetland" & tod_start_ == "day"),
         forest_night = as.numeric(cover == "forest" & tod_start_ == "night"),
         wetland_night = as.numeric(cover == "wetland" & tod_start_ == "night"))
x2_lc_day <- data.frame(forage = mean(obs$forage),
                        temp = mean(obs$temp),
                        predator = mean(obs$predator),
                        cover = factor("grassland",
                                       levels = c("grassland", "forest", "wetland")),
                        tod_start_ = factor("day", levels = c("day", "night")),
                        tod_end_ = factor("day", levels = c("day", "night")),
                        sl_ = 100,
                        log_sl_ = log(100),
                        cos_ta_ = 1) %>% 
  # Make your own dummy variables
  # Remember grassland is our reference level
  mutate(forest_day = as.numeric(cover == "forest" & tod_start_ == "day"),
         wetland_day = as.numeric(cover == "wetland" & tod_start_ == "day"),
         forest_night = as.numeric(cover == "forest" & tod_start_ == "night"),
         wetland_night = as.numeric(cover == "wetland" & tod_start_ == "night"))
# Calculate log-RSS
log_rss_lc_day <- log_rss(issf, x1 = x1_lc_day, x2 = x2_lc_day, ci = "se")

# Night
x1_lc_night <- data.frame(forage = mean(obs$forage),
                        temp = mean(obs$temp),
                        predator = mean(obs$predator),
                        cover = factor(c("grassland", "forest", "wetland"),
                                       levels = c("grassland", "forest", "wetland")),
                        tod_start_ = factor("night", levels = c("day", "night")),
                        tod_end_ = factor("night", levels = c("day", "night")),
                        sl_ = 100,
                        log_sl_ = log(100),
                        cos_ta_ = 1) %>% 
  # Make your own dummy variables
  # Remember grassland is our reference level
  mutate(forest_day = as.numeric(cover == "forest" & tod_start_ == "day"),
         wetland_day = as.numeric(cover == "wetland" & tod_start_ == "day"),
         forest_night = as.numeric(cover == "forest" & tod_start_ == "night"),
         wetland_night = as.numeric(cover == "wetland" & tod_start_ == "night"))
x2_lc_night <- data.frame(forage = mean(obs$forage),
                        temp = mean(obs$temp),
                        predator = mean(obs$predator),
                        cover = factor("grassland",
                                       levels = c("grassland", "forest", "wetland")),
                        tod_start_ = factor("night", levels = c("day", "night")),
                        tod_end_ = factor("night", levels = c("day", "night")),
                        sl_ = 100,
                        log_sl_ = log(100),
                        cos_ta_ = 1) %>% 
  # Make your own dummy variables
  # Remember grassland is our reference level
  mutate(forest_day = as.numeric(cover == "forest" & tod_start_ == "day"),
         wetland_day = as.numeric(cover == "wetland" & tod_start_ == "day"),
         forest_night = as.numeric(cover == "forest" & tod_start_ == "night"),
         wetland_night = as.numeric(cover == "wetland" & tod_start_ == "night"))
# Calculate log-RSS
log_rss_lc_night <- log_rss(issf, x1 = x1_lc_night, x2 = x2_lc_night, ci = "se")

# Combine and plot
bind_rows(log_rss_lc_day$df, log_rss_lc_night$df) %>% 
  ggplot(aes(x = cover_x1, y = log_rss, color = tod_end_x1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(0.25)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(0.25), width = 0.1) +
  scale_color_manual(name = "Time of Day",
                     breaks = c("day", "night"),
                     values = c("goldenrod", "navy")) +
  xlab("Land Cover") +
  ylab("log-RSS vs Grassland") +
  theme_bw()

# We can see that, relative to grassland, forest is selected during the day
# and avoided at night. Wetland is avoided during both day and night.

# We would likely want to temper our conclusions based on the confidence 
# intervals in a real analysis.
  
# Step-length distribution ----
# We also want to plot the step-length distributions for day and night.

expand.grid(sl = seq(0.1, 1000, length.out = 100),
            tod = c("day", "night")) %>% 
  mutate(
    shp = case_when(tod == "day" ~ sl_distr_day$params$shape,
                    tod == "night" ~ sl_distr_night$params$shape),
    scl = case_when(tod == "day" ~ sl_distr_day$params$scale,
                    tod == "night" ~ sl_distr_night$params$scale),
    y = dgamma(sl, shape = shp, scale = scl)
  ) %>% 
  ggplot(aes(x = sl, y = y, color = tod, group = tod)) +
  geom_line(size = 1) +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  scale_color_manual(name = "Time of Day",
                     breaks = c("day", "night"),
                     values = c("goldenrod", "navy")) +
  theme_bw()
  
