#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#------------------Module 06 -- iSSF ------------------X
#----------------Last updated 2022-01-05---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

library(tidyverse)
library(amt)
library(broom)
library(patchwork)
library(NLMR)
library(raster)

library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

set.seed(1323)

source("fun/sim_ssf.R")


# Integrated step-selection analyses are implemented using the following steps:

# 1. Estimate a tentative selection-free movement kernel, using observed step-lengths and turn angles, giving.
# 2. Generate time-dependent available locations by simulating potential movements from the previously observed location.
# 3. Estimate $\beta$ using conditional logistic regression, with strata formed by combining time-dependent used and available locations.
# 4. Re-estimate the movement parameters.

# Simulate data for one animal ----

# ... Landscape ----

# Forest
dim <- 200
forest <- NLMR::nlm_gaussianfield(dim, dim) < 0.5
raster::plot(forest)

# and elevation
ele <- NLMR::nlm_gaussianfield(dim, dim)
ele[] <- scales::rescale(ele[], c(0, 500))
raster::plot(ele)

covars <- raster::stack(forest, ele)
names(covars)
names(covars) <- c("forest", "elevation")
names(covars)

# ... Movement
curve(dexp(x, rate = 0.2), from = 0, to = 20)

dat1 <- simulate_ssf(
  n_steps = 500, n_ch = 10, l = 0.2, xy0 = c(dim/2, dim/2), 
  resc = covars, coef = c(0.01, -0.5)
)

raster::plot(ele)
points(dat1)

raster::plot(forest)
points(dat1)

# Preparing data for SSF ----

# ... Creating steps ----

# From points to steps
dat1 %>% steps() # steps_by_burst()

# Adding random steps. 
tmp <- dat1 %>% steps() %>% random_steps() 

tmp %>% print(n = 15)

# Check the step-length and turn-angle distribution that was fitted to the data.
# We see later how we can plot this.

sl_distr(tmp)
ta_distr(tmp)

# We can vary the number of control steps
dat1 %>% steps() %>% random_steps(n_control = 20) 

# We can now extract the covariates
dat1 %>% steps() %>% 
  random_steps() %>% 
  extract_covariates(covars) 

covars

# Finally, we can some additional covariates
dat1 %>% steps() %>% 
  random_steps() %>% 
  extract_covariates(covars) %>% 
  mutate(log_sl_ = log(sl_), 
         cos_ta_ = cos(ta_))

# Everything at once
ssf.dat <- dat1 %>% steps() %>% 
  random_steps() %>% 
  extract_covariates(covars) %>% 
  mutate(log_sl_ = log(sl_), 
         cos_ta_ = cos(ta_))

ssf.dat

# - `x1_` and `y1_` are the coordinates associated with the starting location of the step.
# - `x2_` and `y2_` are the coordinates associated with the ending location of the step.
# - `sl_` and `ta_` are the step-length and turn angle associated with the step.
# - `t1_` and `t2_` are the timestamps associated with the start and end of the step.
# - `dt_` is the time duration associated with the step.
# - `step_id_` is a unique identifier associated with each step's observed and random locations.
# - `case_` is an indicator variable equal to `TRUE` for observed steps and `FALSE` for random steps.
# - `forest` and `elevation` are the covariate values at the **end** of the step. 


# ... Tentative Movement Parameters   ----
#
#
# Let's have a closer look at what is happening under the hood when applying the
# `random_steps()` function. This function conveniently fits tentative
# step-length and turn-angle distributions and then samples from these
# distributions to generate random steps. The default arguments lead to
# `random_steps()` fitting  a gamma distribution for the step lengths and a von
# Mises distribution for the turn angles.
# 
# It is possible to use other distributions. 
#
# The tentative parameters in these statistical distributions (gamma, von Mises)
# are stored as attributes of the resulting object. We can view the parameters
# of the tentative step-length distribution using:

sl_distr(ssf.dat)

sl_distr(ssf.dat)$params

curve(dgamma(x, shape = sl_distr_params(ssf.dat)$shape, 
             scale = sl_distr_params(ssf.dat)$scale), from = 0, to = 30
)
curve(dexp(x, rate = 0.2), add = TRUE, col = "red")



# Similarly, we can access the parameters of the tentative turn-angle
# distribution using:

ta_distr(ssf.dat)

# In the next module on simulation, we will learn how we can update these
# distributions to adjust for the influence of habitat selection when estimating
# parameters of the movement kernel.
 
# Basic iSSF   ----
# 
# We will fit a model with the same habitat covariates used for the simulations.
# Suppose we hypothesize that habitat selection at the step scale depends on
# forest (`forest`) and elevation (`elevation`). Furthermore, suppose we
# hypothesize that, in the absence of habitat selection, step lengths and turn
# angles follow a constant distribution, i.e., the selection-free movement
# kernel does not depend on covariates. We can fit that model as follows:
 
# Note the use of strata in the model formula. This indicates that each step
# (together with its random steps) forms a strata.
m1 <- ssf.dat %>% 
  fit_issf(case_ ~ forest + elevation +
             strata(step_id_), model = TRUE)


# Interpreting Habitat-Selection Parameters ------
#
# We begin by exploring the coefficients in the fitted model using the
# `summary()` function:

summary(m1)


# ... Calculating Relative Selection Strength (RSS) for Two Locations -----
#
# Calculating the relative use of location $s_1$ versus location $s_2$ is fairly
# straightforward when $s_1$ and $s_2$ share the same values for all but one
# covariate. For more complex scenarios, we have implemented a function in `amt`
# that will calculate the log-relative intensity [referred to as the
# log-Relative Selection Strength, or log-RSS, see Avagar et al 2017]. If you
# prefer to quantify relative-use (i.e., RSS), you can simply exponentiate the
# results.
 
# Here we demonstrate the use  of `log_rss()`: 

#   *  $s_1$: `elevation_end = 100`, `forest = 0`
#   *  $s_2$: `elevation_end = 100`, `forest = 1`

# We want to calculate the log-RSS for a step ending in $s_1$ versus a step
# ending in $s_2$, assuming $s_1$ and $s_2$ are equally accessible to the
# animal.

# The function `log_rss()` expects the two locations to be formatted as separate
# `data.frame` objects. Each `data.frame` must include all covariates used to
# fit the model. Note, Furthermore, `factor` variables should have the same
# `levels` as the original data.
 
# Let us start with creating the two positions.
s1 <- data.frame(
  elevation = 100,
  forest = 1)
 
# data.frame for s2; note the value for forest is different.
s2 <- data.frame(
  elevation = 100,
  forest = 0)

# Now that we have specified each location as a `data.frame`, we can pass them
# along to `log_rss()` for the calculation. The function will return an object
# of class `log_rss`, which is also more generally a `list`. The `list` element
# `"df"` contains a `data.frame` which contains the log-RSS calculation and
# could easily be used to make a plot when considering relative selection
# strength across a range of environmental characteristics.

lr1 <- log_rss(m1, x1 = s1, x2 = s2)

lr1$df

exp(lr1$df$log_rss)

# This is the same as: 

exp(coef(m1)["forest"]) # See Module 5 for this

# `log_rss()` is designed to be able to consider several locations as `x1`,
# relative to a **single** location in `x2`. 

# Lets next consider the scenario where we want to compare locations that are
# outside of forests, but from a range of different elevations

# data.frame for s1
s1 <- data.frame(
  elevation = 10:150,
  forest = 0)
 
# data.frame for s2
s2 <- data.frame(
  elevation = 100,
  forest = 0)

lr1 <- log_rss(m1, x1 = s1, x2 = s2)

head(lr1$df)
exp(lr1$df)

ggplot(lr1$df, aes(elevation_x1, log_rss)) + 
  geom_line() +
  geom_hline(yintercept = 0, col = "red", lty = "dashed")

ggplot(lr1$df, aes(elevation_x1, exp(log_rss))) + 
  geom_line() +
  geom_hline(yintercept = 1, col = "red", lty = "dashed")

# Large-Sample Confidence Intervals
#
# We can use the standard errors from our fitted model to estimate these
# confidence intervals based on a normal approximation to the sampling
# distribution of $\hat{\beta}$.

lr1_ci_se <- log_rss(m1, s1, s2, ci = "se", ci_level = 0.95)
head(lr1_ci_se$df)

ggplot(lr1_ci_se$df, aes(elevation_x1, log_rss)) + 
  geom_hline(yintercept = 0, col = "red", lty = "dashed") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line() 

# If want the RSS instead of log-RSS we have to exponentiate the results again.
# One way to do this is to use `mutate()` in combination with `across`.
lr1_ci_se$df %>% mutate(across(log_rss:upr, exp)) %>% 
  ggplot(aes(elevation_x1, log_rss)) + 
  geom_hline(yintercept = 1, col = "red", lty = "dashed") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line() 


# Fitting an integrated SSF
m2 <- ssf.dat %>% 
  fit_issf(case_ ~ forest + elevation + sl_ + log_sl_ + cos_ta_ +
             strata(step_id_), model = TRUE)

summary(m2)
