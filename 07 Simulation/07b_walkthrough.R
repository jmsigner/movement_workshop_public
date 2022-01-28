#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#----------------Module 07 -- Simulations -------------X
#----------------Last updated 2022-01-21---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

library(tidyverse)
library(amt)
library(broom)
library(patchwork)
library(NLMR)
library(conflicted)
library(raster)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")


set.seed(1323)


# Data ----
dat <- amt_fisher %>% filter(name == "Ricky T")
summarize_sampling_rate(dat)
dat <- dat %>% track_resample(rate = minutes(10), tolerance = seconds(60)) %>% 
  filter_min_n_burst()

# Link to legend:  @TODO
landuse <- raster("07 Simulation/data/landuse_study_area.tif")
plot(landuse)

wet <- landuse == 90
wet
plot(wet)
names(wet)
names(wet) <- "wet"
wet

# Do the same with water
water <- landuse == 11
names(water) <- "water"
plot(water)

covar <- stack(wet, water)
covar

plot(wet)
points(dat)
plot(amt::bbox(dat, buffer = 1e4, spatial = TRUE), add = TRUE)


covar.animal <- crop(covar, amt::bbox(dat, buffer = 1e4, spatial = TRUE))
plot(covar.animal)

# Add distance to water
water <- covar.animal$water
water[] <- ifelse(water[], 1, NA)
plot(water)

waterdist <- distance(water) 
plot(waterdist)

names(waterdist) <- "waterdist"
waterdist <- waterdist / 1e3

covar.animal <- stack(covar.animal, waterdist)
covar.animal
plot(covar.animal)

# Selection-free movement kernel from fitted iSSF ----
issf1 <- dat %>% steps_by_burst() %>% random_steps() %>% 
  extract_covariates(covar.animal) %>% 
  mutate(log_sl_ = log(sl_)) %>% 
  fit_issf(case_ ~ wet + waterdist + sl_ + log_sl_ + strata(step_id_))

summary(issf1)

tentative.sl <- sl_distr(issf1)$params
adjusted.sl <- update_sl_distr(issf1)$params

tentative.sl
adjusted.sl

dev.off()
curve(dgamma(x, shape = tentative.sl$shape, scale = tentative.sl$scale), col = "red", from = 0, to = 200)
curve(dgamma(x, shape = adjusted.sl$shape, scale = adjusted.sl$scale), col = "blue", from = 0, to = 200, add = TRUE)


# Simulate UD with constant selection-free movement kernel -----
covar.sim <- crop(covar.animal, amt::bbox(dat[1, ], buffer = 1e3, spatial = TRUE))
plot(covar.sim)

hk <- habitat_kernel(coef(issf1)[1:2], covar.sim[[c("wet", "waterdist")]])
plot(hk)
mk <- movement_kernel(scale = adjusted.sl$scale, shape = adjusted.sl$shape, hk)
plot(mk)

# ... Transient UD ----
tud <- simulate_tud(mk, hk, start = as.numeric(dat[1, c("x_", "y_")]), n = 100, n_rep = 1e3)
plot(tud)

# ... Steady State UD ----
ssud <- simulate_ud(movement_kernel = mk, habitat_kernel = hk, start = as.numeric(dat[1, c("x_", "y_")]), n = 1e7)
plot(ssud)


# Dynamic simulation ----

# In the simulations above, the step length and turn angle were constant in
# time. That means we can not include interactions between step lengths and turn
# angles with habitat covariates.

# ... Simulate data ----

# We will a simple landscape first. 
lscp <- raster(xmn = -75, xmx = 75, ymn = -75, ymx = 75, res = 1)
lscp[] <- 0
lscp[, 1:75] <- 1
lscp <- stack(lscp)
names(lscp) <- "hab"
plot(lscp)

# Movement characteristics
curve(dgamma(x, scale = 10, shape = 2), from = 0, to = 100, 
      ylab = "Density", xlab = "sl_")

scale_to_sl(10)
shape_to_log_sl(2)

dk1 <- dispersal_kernel(
  ~ sl_ + log_sl_, 
  coefficients = c("sl_" = -0.1, "log_sl_" = 1), 
  spatial.covars = lscp, start = c(0, 0), 
  return.raster = FALSE, max.dist = 70)

dk1 <- dispersal_kernel(
  ~ sl_ + log_sl_, 
  coefficients = c("sl_" = -0.1, "log_sl_" = 1), 
  spatial.covars = lscp, start = c(0, 0), 
  return.raster = TRUE, max.dist = 70)
plot(dk1)

kappa_to_cos_ta(10)
dk2 <- dispersal_kernel(
  ~ cos_ta_, 
  coefficients = c("cos_ta_" = 10), 
  spatial.covars = lscp, start = c(0, 0), 
  return.raster = TRUE, max.dist = 70, direction = -pi/2)
plot(dk2)

dk3 <- dispersal_kernel(
  ~ sl_ + log_sl_ + cos_ta_, 
  coefficients = c("sl_" = -0.1, "log_sl_" = 1, "cos_ta_" = 3), 
  spatial.covars = lscp, start = c(0, 0), direction = pi/3, 
  return.raster = TRUE, max.dist = 70)
plot(dk3)

dk4 <- dispersal_kernel(
  ~ hab_end, 
  coefficients = c("hab_end" = 0.4), 
  spatial.covars = lscp, start = c(0, 0), 
  return.raster = TRUE, max.dist = 70)
plot(lscp)
plot(dk4)

dk5 <- dispersal_kernel(
  ~ cos_ta_ + hab_end, 
  coefficients = c("cos_ta_" = 3, "hab_end" = 1), 
  spatial.covars = lscp, start = c(0, 0), 
  return.raster = TRUE, max.dist = 70)
plot(dk5)

dk6 <- dispersal_kernel(
  ~ cos_ta_ + hab_end + cos_ta_:hab_start, 
  coefficients = c("cos_ta_" = 3, "hab_end" = 0, "cos_ta_:hab_start" = -2.9), 
  spatial.covars = lscp, start = c(-3, 0), 
  return.raster = TRUE, max.dist = 70)
plot(dk6)

# Lets move outside the habitat
dk6 <- dispersal_kernel(
  ~ cos_ta_ + hab_end + cos_ta_:hab_start, 
  coefficients = c("cos_ta_" = 3, "hab_end" = 0, "cos_ta_:hab_start" = -2.9), 
  return.raster = TRUE, max.dist = 70)
plot(dk6)


# ... Real data
issf2 <- dat %>% steps_by_burst() %>% random_steps() %>% 
  extract_covariates(covar.animal, where = "both") %>% 
  mutate(log_sl_ = log(sl_)) %>% 
  fit_issf(case_ ~ wet_end + waterdist_end +  # habitat selection
             sl_ + log_sl_ + # movement
             sl_:waterdist_start + log_sl_:waterdist_start +  # interaction between movement and habitat
             strata(step_id_))
summary(issf2)

sl.w.dist0 <- update_gamma(
  sl_distr(issf2), 
  beta_sl = coef(issf2)["sl_"], 
  beta_log_sl = coef(issf2)["log_sl_"])$params

sl.w.dist100 <- update_gamma(
  sl_distr(issf2), 
  beta_sl = coef(issf2)[c("sl_", "sl_:waterdist_start")] * c(1, 2), 
  beta_log_sl = coef(issf2)[c("log_sl_", "log_sl_:waterdist_start")] * c(1, 2))$params

curve(dgamma(x, shape = sl.w.dist100$shape[2], scale = sl.w.dist100$scale[2]), col = "blue", from = 0, to = 200)
curve(dgamma(x, shape = sl.w.dist0$shape, scale = sl.w.dist0$scale), col = "red", from = 0, to = 200, add = TRUE)

# Now lets create a dispersal kernel for these two situation
cfs <- coef(issf2)

# What could the maximum distance be?
qgamma(0.95, shape = sl_distr_params(issf2)$shape, scale = sl_distr_params(issf2)$scale)
covar.sim <- stack(covar.sim)

start <- as.numeric(dat[1, c("x_", "y_")])

sim.formula <- ~ 
  # Selection
  wet_end + waterdist_end + 
  # Movement
  sl_ + log_sl_ + 
  # Interactions
  sl_:waterdist_start + log_sl_:waterdist_start

sim.formula

sim.coefs <- c(
  # Habitat selection
  "wet_end" = unname(cfs["wet_end"]), 
  "waterdist_end" = unname(cfs["waterdist_end"]),
  # Movement
  "sl_" = unname(cfs["sl_"]) + 1/sl_distr_params(issf2)$scale,
  "log_sl_" = unname(cfs["log_sl_"]) - sl_distr_params(issf2)$shape + 1, 
  # Interaction between both
  "sl_:waterdist_start" = unname(cfs["sl_:waterdist_start"]),
  "log_sl_:waterdist_start" = unname(cfs["log_sl_:waterdist_start"])
)

dk5 <- dispersal_kernel(
  sim.formula,
  coefficients = sim.coefs, 
  spatial.covars = stack(covar.sim), start = start, 
  return.raster = TRUE, max.dist = 350)

plot(dk5)

start1 <- c(1779700, 2412500)
dk5a <- dispersal_kernel(
  sim.formula,
  coefficients = sim.coefs, 
  spatial.covars = stack(covar.sim), start = start1, 
  return.raster = TRUE, max.dist = 350)

plot(dk5a)
plot(dk5)

# We could now repeatedly sample from these dispersal kernels and get track. 

