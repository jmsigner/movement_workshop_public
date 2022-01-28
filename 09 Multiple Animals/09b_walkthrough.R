#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#----------- Module 09 -- Multiple animals ------------X
#----------------Last updated 2021-12-28---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

library(tidyverse)
library(amt)
library(broom)
library(patchwork)
library(NLMR)
library(here)
library(mvtnorm)
library(glmmTMB)

set.seed(1323)

source("fun/sim_ssf.R")

# HSF for multiple animals ----

# ... Goats ----
library(ResourceSelection) # you may have to install this package, I only added today to the list of required packages.
library(broom.mixed)

data(goats, package = "ResourceSelection")

# Here we ignore individuals
m1 <- glmmTMB(STATUS ~ ELEVATION + SLOPE, 
              data = goats, family = binomial())

# This is the random intercept model
m2 <- glmmTMB(STATUS ~ ELEVATION + SLOPE + (1 | ID), 
              data = goats, family = binomial())

# This is a random slope and intercept model
m3 <- glmmTMB(STATUS ~ ELEVATION + SLOPE + 
                (ELEVATION + SLOPE | ID),
              data = goats, family = binomial())

bind_rows(
  tidy(m1, conf.int = TRUE) %>% mutate(what = "glm"),
  tidy(m2, conf.int = TRUE) %>% mutate(what = "glmm (intercept)"),
  tidy(m3, conf.int = TRUE) %>% mutate(what = "glmm (intercept & slope)")
) %>% filter(effect == "fixed") %>% 
  mutate(term1 = str_remove_all(term, "[\\(\\)]"), 
         term1 = str_to_title(term1) %>% factor() %>% fct_inorder(), 
         what = factor(what) %>% fct_inorder()) %>% 
  ggplot(aes(what, estimate)) + 
           geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  facet_wrap(~ term1, scale = "free") + 
  geom_hline(yintercept = 0, col = "red", lty = 2) + 
  labs(y = "Estimate", x = "Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank())


# Simulations ----

# ... Landscape ----
set.seed(1234301)

# Forest
dim <- 500
forest <- NLMR::nlm_gaussianfield(dim, dim) < 0.5
raster::plot(forest)

# and elevation
ele <- NLMR::nlm_gaussianfield(dim, dim)
raster::plot(ele)

covars <- raster::stack(forest, ele)
names(covars)
names(covars) <- c("forest", "elevation")

# ... Movement
curve(dexp(x, rate = 0.5), from = 0, to = 20)

dat1 <- simulate_ssf(
  n_steps = 500, n_ch = 10, l = 0.5, xy0 = c(dim/2, dim/2), 
  resc = covars, coef = c(0.1, 0)
)

raster::plot(covars)
raster::plot(ele)
points(dat1)

# Now let us simulate 10 animals from the same populations. 
coefs <- rmvnorm(n = 10, mean = c(0.01, -0.5), sigma = diag(c(0.00001, 0.2)))
dat <- map(1:10, ~ {
  simulate_ssf(
    n_steps = 500, n_ch = 10, l = 0.5, xy0 = c(dim/2, dim/2), 
    resc = covars, coef = coefs[.x, ]
  ) %>% mutate(id = .x)
})

dat1 <- dat %>% bind_rows()
dat1 %>% ggplot(aes(x_, y_, col = factor(id))) + geom_point() + coord_equal()

# HSF ----
# ... Prepare data ----

dat.hsf <- dat1 %>% nest(data = -id) %>% 
  mutate(random.points = map(data, random_points)) %>% 
  select(-data) %>% unnest(cols = random.points) %>% 
  extract_covariates(covars) %>% 
  mutate(weight = ifelse(case_, 1, 1e5))

dat.hsf

# ... Global model ----

m1 <- glm(case_ ~ forest + elevation, data = dat.hsf, weights = weight, family = binomial())
summary(m1)

# ... Individual model ----
m2 <- dat.hsf %>% nest(dat = -id) %>% 
  mutate(mod = map(dat, ~ glm(case_ ~ forest + elevation, data = ., weights = weight, family = binomial()) %>% 
                     tidy(conf.int = TRUE))) %>% 
  select(id, mod) %>% unnest(cols = mod)
m2 %>% filter(term != "(Intercept)") %>% 
  ggplot(aes(x = term, y = estimate)) + geom_point(alpha = 0.4) +
  stat_summary(col = "red")

# ... Random intercept ----
m3 <- glmmTMB(case_ ~ forest + elevation + (1 | id), data = dat.hsf, weights = weight, family = binomial())
summary(m3)

# ... Random  slope ----
m4 <- glmmTMB(case_ ~ forest + elevation + (0 + forest | id) + 
                (0 + elevation | id), 
              data = dat.rsf, weights = weight, family = binomial())
summary(m4)


# True parameters under an SSF

# coefs <- rmvnorm(n = 10, mean = c(0.01, -0.5), sigma = diag(c(0.00001, 0.2)))
# Forest: 0.01
# Elevation: -0.5

confint(m1) # ignoring individual animals
confint(m3) # Random intercept
confint(m4) # Random slope


# iSSF ----

# ... Prepare data ----
dat.issf <- dat1 %>% nest(data = -id) %>% 
  mutate(random.steps = map(data, ~ steps(.x) %>% random_steps)) %>% 
  select(-data) %>% unnest(cols = random.steps) %>% 
  extract_covariates(covars) %>% 
  mutate(log_sl_ = log(sl_), sin_ta_ = sin(ta_), 
         step_id1_ = paste(id, step_id_, sep = "."))

# ... Poisson trick ----

# Lets use data for just one individual
dat.issf.1 <- filter(dat.issf, id == 1)

# Use clogit
m1 <- fit_issf(dat.issf.1, case_ ~ forest + elevation + sl_ + log_sl_ + strata(step_id1_))
summary(m1)

# Poisson formulation 
m2 <- glmmTMB(case_ ~ -1 + forest + elevation + sl_ + log_sl_ + 
                ## Make sure we use a unique step id here
                (1|step_id1_), 
               family = poisson(), data = dat.issf.1, doFit = FALSE)
m2$parameters$theta[1] <- log(1e3) # fixing the variance for stratum specific random intercepts
m2$mapArg <- list(theta = factor(NA))
m2 <- glmmTMB:::fitTMB(m2)

# Estimated coefficients are identical
coef(m1)
fixef(m2)

# So are SE
vcov(m1$model)
vcov(m2)

# ... Global model ----

m1 <- fit_issf(case_ ~ forest + elevation + sl_ + log_sl_ + strata(step_id1_), data = dat.issf)
summary(m1)

# ... Individual model ----
m2 <- dat.issf %>% nest(dat = -id) %>% 
  mutate(mod = map(dat, ~ {
    m <- fit_issf(case_ ~ forest + elevation + sl_ + log_sl_ + strata(step_id1_), data = .) 
    tidy(m$model, conf.int = TRUE)
    })) %>% 
  select(id, mod) %>% unnest(cols = mod)
m2 %>% 
  ggplot(aes(x = term, y = estimate)) + geom_point(alpha = 0.4) +
  stat_summary(col = "red")

# ... Random intercept & random slope ----

m3 <- glmmTMB(case_ ~ -1 + forest + elevation + sl_ + log_sl_ + (1 |step_id1_) +
                (0 + elevation | id) + (0 + forest | id), 
               family = poisson(), data = dat.issf, doFit = FALSE)
m3$parameters$theta[1] <- log(1e3)
m3$mapArg <- list(theta = factor(c(NA, 1:2)))
m3 <- glmmTMB:::fitTMB(m3)

summary(m3)

# Comparing again models

# coefs <- rmvnorm(n = 10, mean = c(0.01, -0.5), sigma = diag(c(0.00001, 0.2)))
# Forest: 0.01
# Elevation: -0.5

confint(m1$model) # ignoring individual animals
confint(m3) # Random intercept

summary(m1)
summary(m3)

# Compare Individual models and random-effects model

# Deviation from the FE
ranef(m3)[[1]]$id

# Fixedef
coef(m3)

coef.all <- bind_rows(
  "ind" = m2 %>% select(term, estimate) %>% 
    filter(!term %in% c("sl_", "log_sl_")) %>% mutate(animal = rep(1:10, each = 2)),
  "RE" = coef(m3)[[1]]$id[, c("forest", "elevation")] %>% pivot_longer(
    cols = elevation:forest, 
    names_to = "term", values_to = "estimate") %>% mutate(animal = rep(1:10, each = 2)), 
  .id = "what"
) 

coef.all %>% ggplot(aes(animal, estimate, col = what)) +
  geom_point() +
  facet_wrap(~ term, scale = "free") +
  theme_light()


coef.all %>% pivot_wider(names_from = "term", values_from = "estimate") %>% 
  ggplot(aes(forest, elevation)) + 
  geom_hline(yintercept = -0.5, col = "red") +
  geom_vline(xintercept = 0.01, col = "red") +
  geom_line(aes(group = animal), col = "grey") +
  geom_point(aes(col = what)) +
  theme_light()
