#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#-------------------Module 05 -- HSFs------------------X
#----------------Last updated 2021-01-17---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

# Load packages ----
library(tidyverse)
library(raster)
library(amt)

# Generate example data ----
# We need to generate some example telemetry data to fit our HSF.
# Generating points is fairly easy since the HSF assumes all points are 
# independent. However, we should spend some time thinking about how
# we define the relationships between our study animal and habitat.

# ... habitat axes ----
# First we need habitat variables to define our e-space.

# It is often useful to break our variables into:
#   - resources (more is better)
#   - conditions (some intermediate value is best)
#   - risk (less is better)

# I've already randomly generated 4 variables:
#   - forage (g/m^2) [resource]
#   - mean annual temperature (°C) [condition]
#   - predator density (predators/100 km^2) [risk]
#   - land cover type (grassland = 1, forest = 2, wetland = 3) [condition]

# If you want to see how I did that, check out the script "habitat.R". We can 
# load those here as a RasterStack.

hab <- stack("05 HSF/geo/habitat.tif")
names(hab) <- c("forage", "temp", "predator", "cover")
plot(hab)

# Cover is a factor, let's code it that way
values(hab$cover) <- factor(values(hab$cover), 
                            levels = 1:3,
                            labels = c("grassland", "forest", "wetland"))
hab$cover

# ... habitat relationships ----
# Now we have to define our relationship between density and our habitat axes.
# In general:
#   - resources should have a single, positive slope parameter (more is better)
#   - risks should have a single, negative slope parameter (less is better)
#   - conditions should be defined by a concave-down parabola (intermediate is best)
#       - a parabola requires a linear and quadratic term
#       - the x-coordinate of the vertex will be at -b/2a
#           - a is the beta for the quadratic term
#           - b is the beta for the linear term

# Note: we are ignoring the intercept here, so the y-axis is just 
# *relative* density.

# We will also assume (for now) that our habitats are strictly additive
# (no interactions).

# ... ... resources ----
# Let's begin with resources.
# We know we want our beta to be positive. How big should it be?
hist(hab$forage)

# Our forage variable ranges from 0 to 1000 (g/m^2). Let's say that at 1000 
# g/m^2, the density is 5x the density at 500 g/m^2. I.e., RSS(x1, x2) = 5 if
# x1 has forage == 1000 and x2 has forage == 500. 

# The beta is just the log of RSS for a 1-unit change, and it is also the
# difference between the linear predictors. I.e.,
# log(5) = (beta_forage * 1000) - (beta_forage * 500)
#   ==>
beta_forage = log(5)/500

# Let's check our work. The RSS for 1000 vs 500 g/m^2 should be 5.
exp(beta_forage * 1000) / exp(beta_forage * 500)

# How many times more animals do we have at 500 g/m^2 than 200 g/m^2?
exp(beta_forage * 500) / exp(beta_forage * 200)

# How many times more animals do we have at 1000 g/m^2 than 0 g/m^2?
exp(beta_forage * 1000) / exp(beta_forage * 0)

# ... ... risks ----
# Now onto risks.
# We know we want our beta to be negative, but again, how big?
hist(hab$predator)

# Our predator density ranges from 0 - 12 (predators / 100 km^2).
# Let's say that we have 0.25x as many animals at predator == 10 than
# predator == 5. I.e., RSS(x1, x2) = 0.25.

# As above,
# log(0.25) = (beta_pred * 10) - (beta_pred * 5)
#   ==>
beta_pred = log(0.25)/5

# Let's check our work.
exp(beta_pred * 10) / exp(beta_pred * 5)

# How many MORE animals will we have if predator density is 0 than if predator
# density is 12?
exp(beta_pred * 0) / exp(beta_pred * 12)

# Note that the overall change in density of our critter changes ~ 25x as we
# move through the full range of the risk and resource variables. But the
# magnitude of 'beta_forage' is ~ 100x smaller than the magnitude of 
# 'beta_pred'. That's because the range of forage is about 100x greater than
# the range of predator densities. I.e., be careful not to interpret the
# magnitude of the beta if you haven't standardized your covariates.

# ... ... conditions ----
# Now onto conditions.
# This one is tricky because we are trying to parameterize a parabola, so
# we need 2 parameters. Let's start by deciding where the vertex (the best
# habitat) should be.
hist(hab$temp)

# Say our critter prefers a mean annual temperature of 13 degrees C, i.e.,
# -b/2a = 13, or to give our parameters better names,
# (-1 * beta_temp) / (2 * beta_temp2) = 13
#   ==>
# beta_temp/beta_temp2 = 13 * -2 = -26

# Now we need a second piece of information to make the parameters identifiable.
# So let's say our our animal density is double at 13 degrees vs 7 degrees.
# I.e., RSS(x1, x2) = 2
#   ==>
# (beta_temp2 * 13^2 + beta_temp * 13) - (beta_temp2 * 7^2 + beta_temp * 7) = log(2)

# I won't try to type out the algebra here, but if you solve this system
# of equations, you'll find that:
beta_temp2 = -1 * log(2)/36
beta_temp = beta_temp2 * -26

# Let's check our work.
exp(beta_temp2 * 13^2 + beta_temp * 13) / exp(beta_temp2 * 7^2 + beta_temp * 7)

# How many times more animals will we have at temp == 10 than temp == 5?
exp(beta_temp2 * 10^2 + beta_temp * 10) / exp(beta_temp2 * 5^2 + beta_temp * 5)

# How many times more animals will we have at temp == 13 than temp == 2.5?
exp(beta_temp2 * 13^2 + beta_temp * 13) / exp(beta_temp2 * 2.5^2 + beta_temp * 2.5)

# ... ... land cover ----
# Finally, we tackle land cover.
# One of the categories will be the reference category. I.e., it would be
# captured by the "intercept", but since we don't have one, it takes on the
# values from all our other variables. We'll make grassland our reference.

# The beta for forest is then the log-RSS for forest vs. grassland.
# The beta for wetland is then the log-RSS for wetland vs grassland.

# Let's say density is twice as high in forest as in grassland. I.e.,
beta_forest <- log(2)

# Let's say density is half as high in wetland as in grassland. I.e.,
beta_wetland <- log(1/2)

# Okay, we have all of our betas!

# ... calculate w(x) ----
# Recall that w(x) will be proportional to our number of points. We will decide
# later how many total data points to collect, but for now, let's calculate
# w(x) for all locations in space.

# Get our raster data into a data.frame
dat <- as.data.frame(hab, xy = TRUE) %>% 
  rename(cover = cover_VALUE) %>% 
  # Calculate g(x)
  mutate(g = 
           # forage
           beta_forage * forage +
           # two terms for temperature
           beta_temp * temp +
           beta_temp2 * temp^2 +
           # predator density
           beta_pred * predator +
           # landcover
           beta_forest * (cover == "forest") +
           beta_wetland * (cover == "wetland")) %>% 
  # Calculate w(x)
  mutate(w = exp(g))

# ... calculate expected number of points ----
# Let's say that all of our simulated data will come from a single individual
# wearing a GPS collar.
#   - Assume our GPS fixes are far enough apart that there is no autocorrelation
#   - Assume the entire raster extent is equally available.

# How many points should we have? Let's say we've deployed our collar long
# enough to have ~ 1000 locations.

# What's the expected number of points in each cell?

# Recall, under the IPP, lambda(s) is proportional to w(x(s)). Our expected
# number of points in each raster cell is 1000 * w'(x), where w'(x) is a 
# normalized version of w(x), i.e., it sums to 1.
dat$w_prime <- dat$w/sum(dat$w)

# Expected points in each cell
dat$lambda <- 1000 * dat$w_prime
hist(dat$lambda)

# ... draw points ----
# Finally, we can draw n, the number of locations in each cell. Recall, this
# is a Poisson IPP, so we will draw from the Poisson distribution.
set.seed(20220126)
dat$n <- rpois(n = nrow(dat), lambda = dat$lambda)

# This is all we need if we assume we would use the Poisson GLM to fit our
# HSF.
dat %>% 
  dplyr::select(x, y, forage, temp, predator, cover, n) %>% 
  head()

# However, if we really want GPS coordinates, we can take the coordinates of
# each cell with n > 1, and randomly place n points within that cell. Let's do 
# that.

# Our coordinates are the center of our raster cell, and each cell is 50m x 50m.
# So we can jitter up to 25m in any direction and still be in the same cell.

# Function to jitter data
jitter <- function(x, y, min = -25, max = 25) {
  res <- data.frame(x = x + runif(1, min, max),
                    y = y + runif(1, min, max))
  return(res)
}

# E.g.,
jitter(0, 0)

# Now we split each row with n > 1 into a list element
dat_list <- dat %>% 
  filter(n > 0) %>% 
  split(1:nrow(.))

# And now we can create jittered points for each element of our list.
# We will use 'bind_rows()' (twice) to return a single data.frame
set.seed(20220126)
gps <- lapply(dat_list, function(d) {
  replicate(d$n, jitter(x = d$x, y = d$y), simplify = FALSE) %>% 
    bind_rows()
}) %>% 
  bind_rows()

# This is what the result looks like.
head(gps)

plot(hab$forage)
points(gps$x, gps$y, pch = 16, cex = 0.5)

# Fitting an HSF with 'amt' ----
# Now that we've simulated these data, let's fit an HSF.
mod_dat <- gps %>% 
  make_track(x, y, crs = 32612) %>% 
  # This is technically sampling within a 100% MCP by default,
  # but that is practically the extent of our raster.
  random_points(n = nrow(gps) * 100) %>% 
  extract_covariates(hab) %>% 
  mutate(cover = factor(cover,
                        levels = 1:3,
                        labels = c("grassland", "forest", "wetland"))) %>% 
  # Assign large weights to available points
  mutate(weight = ifelse(case_, 1, 1e5))

# Fit a model
hsf <- glm(case_ ~ forage + temp + I(temp^2) + predator + cover,
           family = binomial(), weights = weight, data = mod_dat)

# Summary
summary(hsf)

# Compare our fitted coefficients to true values
b <- coef(hsf)
ci <- confint(hsf)

res <- as.data.frame(cbind(b, ci))
res$truth <- c(NA,
               beta_forage,
               beta_temp,
               beta_temp2,
               beta_pred,
               beta_forest,
               beta_wetland)
names(res) <- c("est", "lwr", "upr", "truth")

# Plot
(beta_fig <- res %>% 
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
    theme_bw())

# Zoom in to see forage and temp^2
beta_fig +
  coord_cartesian(ylim = c(-0.025, 0.005))

# Looks great!

# Model interpretation ----
# We probably have a pretty good idea of what's going on in our model since
# we simulated the data, but this is a crucial step in understanding what our
# fitted model is telling us.

# Using RSS is the best way to understand what's going on biologically.

# We can calculate log-RSS using a function from `amt`.

?log_rss

# *Importantly*, the data.frame 'x2' should almost always have only 1-row
# to avoid any confusion with R's vector-recycling rules.

# ... scenarios ----
# One useful thing we can do with RSS is to generate some biologically 
# interesting scenarios and ask how many times more points we expect
# in one habitat (real or hypothetical) vs another habitat.

# For example, we might be interested in the tradeoff between foraging
# and predation risk. We could create a scenario with high forage and low
# predator density vs. a habitat with low forage and high predor density.
# How many times more points do we expect if:
#   - x1: 
#       - forage: 750 g/m^2
#       - predator density: 3 predators/100 km^2
#   - x2: 
#       - forage: 250 g/m^2
#       - predator density: 8 predators/100 km^2

# Note that we didn't specify temperature or land cover. We still have to
# pass these variables to 'log_rss()', but *IT DOES NOT MATTER* what values
# we choose for x1 and x2, as long as they are the same (because our model
# has no interactions).

# Define x1
x1 <- data.frame(forage = 750,
                 predator = 3,
                 temp = 15,
                 # Note factors need the same levels as the original data
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

# Define x2
x2 <- data.frame(forage = 250,
                 predator = 8,
                 temp = 15,
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

# Calculate log-RSS and 95% CI
scenario_rss <- log_rss(hsf, x1, x2, ci = "se")

# Examine the structure of the resulting object
str(scenario_rss)

# The first object, "df", is typically what you'll want to work with.
scenario_rss$df

# The data.frame shows us the values at x1 and the resulting log-RSS.
# Exponentiate to answer our question.
scenario_rss$df %>% 
  dplyr::select(log_rss:upr) %>% 
  mutate(across(everything(), exp))

# We expect almost 20x as many points if forage is high and predation is low
# vs. if forage is low and predation is high. We can also see the 95% CI
# for that estimate.

# You can see how generating interesting scenarios with combinations of 
# different habitats can be interesting.

# Another common approach is to use RSS to look across a range of values
# for one habitat axis, with all others held constant. Let's do that for each
# of our habitat dimensions.

# *Note* that I am re-using the names 'x1' and 'x2' in each example. This
# works fine if you always run the script top-to-bottom, but it gets quite
# dangerous if you are skipping around and running the script interactively.
# Be careful if you decide not to give each x1 and x2 unique names!

# ... forage ----
# x1 is a sequence across the whole range of forage
# All other variables held constant
x1 <- data.frame(forage = seq(0, 1000, length.out = 100),
                 temp = mean(mod_dat$temp),
                 predator = mean(mod_dat$predator),
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

# x2 is still a single row
x2 <- data.frame(forage = mean(mod_dat$forage),
                 temp = mean(mod_dat$temp),
                 predator = mean(mod_dat$predator),
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

# Calculate log-RSS
forage_rss <- log_rss(hsf, x1, x2, ci = "se")

# Plot RSS
(forage_plot <- forage_rss$df %>% 
    mutate(rss = exp(log_rss),
           exp_lwr = exp(lwr),
           exp_upr = exp(upr)) %>% 
    ggplot(aes(x = forage_x1, y = rss, ymin = exp_lwr, ymax = exp_upr)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_ribbon(size = 0.5, linetype = "dashed", fill = "gray70", color = "black") +
    geom_line(size = 1) +
    xlab(expression("Forage " * (g/m^2))) +
    ylab("RSS") +
    theme_bw())

# ... temperature ----
x1 <- data.frame(forage = mean(mod_dat$forage),
                 temp = seq(2, 20, length.out = 100),
                 predator = mean(mod_dat$predator),
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

x2 <- data.frame(forage = mean(mod_dat$forage),
                 temp = mean(mod_dat$temp),
                 predator = mean(mod_dat$predator),
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

temp_rss <- log_rss(hsf, x1, x2, ci = "se")

(temp_plot <- temp_rss$df %>% 
    mutate(rss = exp(log_rss),
           exp_lwr = exp(lwr),
           exp_upr = exp(upr)) %>% 
    ggplot(aes(x = temp_x1, y = rss, ymin = exp_lwr, ymax = exp_upr)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_ribbon(size = 0.5, linetype = "dashed", fill = "gray70", color = "black") +
    geom_line(size = 1) +
    xlab("Temperature (°C)") +
    ylab("RSS") +
    theme_bw())

# ... predator density ----
x1 <- data.frame(forage = mean(mod_dat$forage),
                 temp = mean(mod_dat$temp),
                 predator = seq(0, 12, length.out = 100),
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

x2 <- data.frame(forage = mean(mod_dat$forage),
                 temp = mean(mod_dat$temp),
                 predator = mean(mod_dat$predator),
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

pred_rss <- log_rss(hsf, x1, x2, ci = "se")

(pred_plot <- pred_rss$df %>% 
    mutate(rss = exp(log_rss),
           exp_lwr = exp(lwr),
           exp_upr = exp(upr)) %>% 
    ggplot(aes(x = predator_x1, y = rss, ymin = exp_lwr, ymax = exp_upr)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_ribbon(size = 0.5, linetype = "dashed", fill = "gray70", color = "black") +
    geom_line(size = 1) +
    xlab(expression("Predator Density " * (pred/km^2))) +
    ylab("RSS") +
    theme_bw())

# ... land cover ----
x1 <- data.frame(forage = mean(mod_dat$forage),
                 temp = mean(mod_dat$temp),
                 predator = mean(mod_dat$predator),
                 cover = factor(c("grassland", "forest", "wetland"), 
                                levels = c("grassland", "forest", "wetland")))

x2 <- data.frame(forage = mean(mod_dat$forage),
                 temp = mean(mod_dat$temp),
                 predator = mean(mod_dat$predator),
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

cover_rss <- log_rss(hsf, x1, x2, ci = "se")

(cover_plot <- cover_rss$df %>% 
    mutate(rss = exp(log_rss),
           exp_lwr = exp(lwr),
           exp_upr = exp(upr)) %>% 
    ggplot(aes(x = cover_x1, y = rss, ymin = exp_lwr, ymax = exp_upr)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_errorbar(width = 0.25, color = "black") +
    geom_point(size = 2) +
    xlab("Land Cover") +
    ylab("RSS") +
    theme_bw())

# Mapping the HSF ----
# We can use the HSF to predict across space now.
# We may want to plot w(x). We can use R's generic 'predict()' function
# to do this, but recall, we don't want the intercept. Let's call the 
# linear prediction with the intercept y(x).
map_dat <- as.data.frame(hab, xy = TRUE) %>% 
  rename(cover = cover_VALUE) %>% 
  mutate(y_x = predict(hsf, newdata = .),
         # Subtract off the intercept to get g(x)
         g = y_x - coef(hsf)[1],
         # Exponentiate to get w(x)
         w = exp(g))

# Map w(x)
ggplot(map_dat, aes(x = x, y = y, fill = w)) +
  geom_raster() +
  scale_fill_viridis_c(name = expression(w(x))) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw()

# Interpreting w(x) is difficult here, so it might be more intuitive to plot
# RSS vs mean conditions.

x1 <- as.data.frame(hab, xy = TRUE) %>% 
  rename(cover = cover_VALUE) %>% 
  mutate(cover = factor(cover, 
                        levels = c("grassland", "forest", "wetland")))

x2 <- data.frame(forage = mean(mod_dat$forage),
                 temp = mean(mod_dat$temp),
                 predator = mean(mod_dat$predator),
                 cover = factor("grassland", 
                                levels = c("grassland", "forest", "wetland")))

map_rss <- log_rss(hsf, x1, x2)

(map <- map_rss$df %>% 
    mutate(rss = exp(log_rss)) %>% 
    ggplot(aes(x = x_x1, y = y_x1, fill = rss)) +
    geom_raster() +
    scale_fill_viridis_c(name = "RSS") +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw())

# The new map looks the same because the RSS is linearly proportional to w(x),
# but the values are easier to understand intuitively now.

# The problem with this is that we can't really tell what's going on below 1.
# log-RSS would be better for that

(map <- map_rss$df %>% 
    mutate(rss = exp(log_rss)) %>% 
    ggplot(aes(x = x_x1, y = y_x1, fill = log_rss)) +
    geom_raster() +
    scale_fill_gradient2(name = "log-RSS",
                        low = "navyblue",
                        mid = "white",
                        high = "firebrick",
                        midpoint = 0) +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw())

# By using a color palette that distinguishes between negative, positive, and
# 0, we can at a glance see if a pixel has habitats that would be used more than
# expected, less than expected, or directly in proportion to their availability.