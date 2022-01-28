library(tidyverse)
library(lubridate)
library(amt)
library(glmmTMB)

library(conflicted)

conflict_prefer("select", "dplyr")



# Q1: Load the fisher data from the `amt` package and sample all individuals to the same sampling rate. 
# You can use the function `summarize_sampling_rate_many()` to figure out a suitable sampling rage.
# For the actual resampling you need a loop/nesting or alike. 

trk <- amt_fisher


summarize_sampling_rate_many(trk, "id")

dat1 <- trk %>% nest(data = -id)

dat2 <- dat1 %>% 
  mutate(dat.resample = map(data, ~ track_resample(., rate = minutes(30), tolerance = minutes(2))))
dat2

# Q2: Fit four HSFs using only the covariate elevation. Fit the following three models:
# 1. Ignore individuals
# 2. A model with a random intercept only
# 3. A model only with a random slope

# Get data ready

dat.hsf <- dat2 %>% mutate(
  dat.hsf = map(dat.resample, ~ random_points(.x))
)

dat.hsf1 <- dat.hsf %>% unnest(cols = dat.hsf) %>% 
  extract_covariates(amt_fisher_covar$elevation) %>% 
  mutate(weight = ifelse(case_, 1, 1e2)) %>% 
  select(case_, elevation, weight, id)

# 1. Ignore the individuals
m1 <- glm(case_ ~ elevation, data = dat.hsf1, family = binomial(), weights = weight)
summary(m1)

# 2. A model with a random intercept only
m2 <- glmmTMB(case_ ~ elevation + (1|id), data = dat.hsf1, family = binomial(), weights = weight)
summary(m2)

# 3. A model only with a random slope
m3 <- glmmTMB(case_ ~ elevation + (elevation|id), data = dat.hsf1, family = binomial(), weights = weight)
summary(m3)


# Q3: Use the same data set, but now fit an SSF. Fit two different models:
# 1. A conditional logistic regression ignoring different individuals.
# 2. Fit a model where you use the Poisson trick to account for individual variation. 
 
# Prepare the data
dat.ssf <- dat2 %>% 
  mutate(ssf = map(dat.resample, ~ .x %>% steps_by_burst %>% 
                     random_steps %>% extract_covariates(amt_fisher_covar$elevation))) 

dat.ssf2 <- dat.ssf %>% select(ssf, id) %>% unnest(cols = ssf) %>% 
  mutate(step_id1_ = paste0(id, "-", step_id_))

m.ssf0 <- fit_ssf(case_ ~ elevation + strata(step_id1_), data = dat.ssf2)

# Apply the poisson trick
m.ssf1 <- glmmTMB(case_ ~ -1 + elevation + (1 | step_id1_) + (0 + elevation|id),
                  family = poisson(), 
                  data = dat.ssf2, doFit = FALSE)

# Set variance of random intercept to 10^6
m.ssf1$parameters$theta[1] <- log(1e3)
m.ssf1$mapArg <- list(theta = factor(c(NA, 1)))
m.ssf1 <- glmmTMB:::fitTMB(m.ssf1)

# Compare the results
coef(m.ssf0)
fixef(m.ssf1)

confint(m.ssf0$model)
confint(m.ssf1)

vcov(m.ssf0$model)
vcov(m.ssf1)

