#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#----------------Module 10 -- Validation---------------X
#----------------Last updated 2021-01-27---------------X
#-------------------Exercise Solution------------------X
#######################################################X

# Refitting the cougar iSSF from module 8.

# Load packages ----
library(tidyverse)
library(amt)
library(lubridate)
library(raster)

# 1. Load data ----
# Location data
dat <- read_csv("data/coyote_cougar.csv") %>% 
  # Subset to just cougar F53
  filter(id == "F53")

# Set the timezone
tz(dat$t_) <- "US/Mountain"

# Habitat data as stack
hab <- stack("data/coyote_cougar_habitat.tif")
names(hab) <- c("elevation", "trees", "biomass", "dist_to_road")

# Prep data
set.seed(1)
mod_dat <- dat %>% 
  # Make track
  make_track(x_, y_, t_, crs = 32612) %>% 
  # Find regular bursts
  track_resample(rate = hours(4), tolerance = minutes(15)) %>% 
  filter_min_n_burst(min_n = 3) %>% 
  # Convert points to steps
  steps_by_burst() %>% 
  # Generate available steps
  random_steps(n_control = 20) %>% 
  # Attach habitat variables
  extract_covariates(hab) %>% 
  # Add additional movement covariates
  mutate(log_sl_ = log(sl_),
         cos_ta_ = cos(ta_))

head(mod_dat)

# 2. Separate into training and testing ----
set.seed(2)
step_df <- data.frame(stratum = sort(unique(mod_dat$step_id_))) %>% 
  mutate(train = as.logical(rbinom(n = nrow(.), size = 1, prob = 0.8)))

train <- mod_dat %>% 
  filter(step_id_ %in% step_df$stratum[step_df$train])

test <- mod_dat %>% 
  filter(step_id_ %in% step_df$stratum[!step_df$train])

nrow(train)/nrow(mod_dat)
nrow(test)/nrow(mod_dat)

# 3. Fit models ----
m1 <- train %>% 
  fit_issf(case_ ~ 
             # Habitat selection main effects
             elevation + I(elevation^2) + trees + dist_to_road +
             # Movement main effects
             sl_ + log_sl_ + cos_ta_ +
             # Don't forget the strata
             strata(step_id_),
           # And include model = TRUE so we can use 'log_rss()' later
           model = TRUE)

summary(m1)

# Drop trees
m2 <- train %>% 
  fit_issf(case_ ~ 
             # Habitat selection main effects
             elevation + I(elevation^2) + dist_to_road +
             # Movement main effects
             sl_ + log_sl_ + cos_ta_ +
             # Don't forget the strata
             strata(step_id_),
           # And include model = TRUE so we can use 'log_rss()' later
           model = TRUE)

summary(m2)

# Calculate metrics ----

# ... concordance ----
# Concordance is already in the summary
# Analogous to AUC
summary(m1)$concordance
summary(m2)$concordance

# ... UHC plots ----

# ... Step 1. Summarize distribution of covariates ----

# Plot them all at once with ggplot
test %>% 
  pivot_longer(elevation:dist_to_road) %>% 
  ggplot(aes(x = value, color = case_)) +
  facet_wrap(~ name, scales = "free") +
  geom_density() +
  theme_bw()

# ... Step 2. Fit a model to the training data ----

# We already did this, but now we also need the betas from our model and 
# the variance-covariance matrix
m1_betas <- coef(m1)
m1_vcov <- vcov(m1$model)

m2_betas <- coef(m2)
m2_vcov <- vcov(m2$model)

# ... Step 3. Create predicted distribution using fitted model ----

# We are going to create the predicted distribution by repeating the
# substeps many times. Let's try it with 500 iterations (you probably
# want more for your own data).

M <- 500

# Initialize blank lists
dens_elevation_m1 <- list()
dens_trees_m1 <- list()
dens_dist_to_road_m1 <- list()

for (i in 1:M) {
  # Report status
  cat("\rIteration", i, "of", M, "       ")
  # ... ... a. draw random values for the betas ----
  
  # We are going to draw new betas from a multivariate normal distribution.
  # Base R doesn't have a random number generator for the multivariate normal,
  # but the MASS package does.
  samp_beta <- MASS::mvrnorm(n = 1, mu = m1_betas, Sigma = m1_vcov)
  
  # ... ... b. select points from test data ----
  
  # Now we sample from the test points with probability proportional
  # to the exponential habitat selection function.
  
  # Difference for iSSF from HSF is that we stratify random sampling
  # by stratum and select one observation per stratum.
  samp_test <- test %>% 
    # Calculate w(x), i.e., the exponential HSF
    mutate(w = exp(samp_beta[1] * elevation +
                     samp_beta[2] * elevation^2 +
                     samp_beta[3] * trees +
                     samp_beta[4] * dist_to_road)) %>% 
    # Grouping seems to be failing on this object; clearing classes
    as.data.frame() %>% 
    group_by(step_id_) %>% 
    slice_sample(n = 1, weight_by = w)
  
  # ... ... c. summarize covariates at these locations ----
  dens_elevation_m1[[i]] <- density(samp_test$elevation)
  dens_trees_m1[[i]] <- density(samp_test$trees)
  dens_dist_to_road_m1[[i]] <- density(samp_test$dist_to_road)
}

# Same for m2
# Initialize blank lists
dens_trees_m2 <- list()
dens_elevation_m2 <- list()
dens_dist_to_road_m2 <- list()

for (i in 1:M) {
  # Report status
  cat("\rIteration", i, "of", M, "       ")
  # ... ... a. draw random values for the betas ----
  
  # We are going to draw new betas from a multivariate normal distribution.
  # Base R doesn't have a random number generator for the multivariate normal,
  # but the MASS package does.
  samp_beta <- MASS::mvrnorm(n = 1, mu = m2_betas, Sigma = m2_vcov)
  
  # ... ... b. select points from test data ----
  
  # Now we sample from the test points with probability proportional
  # to the exponential habitat selection function.
  
  # Difference for iSSF from HSF is that we stratify random sampling
  # by stratum and select one observation per stratum.
  samp_test <- test %>% 
    # Calculate w(x), i.e., the exponential HSF
    mutate(w = exp(samp_beta[1] * elevation +
                     samp_beta[2] * elevation^2 +
                     samp_beta[3] * dist_to_road)) %>% 
    # Grouping seems to be failing on this object; clearing classes
    as.data.frame() %>% 
    group_by(step_id_) %>% 
    slice_sample(n = 1, weight_by = w)
  
  # ... ... c. summarize covariates at these locations ----
  dens_trees_m2[[i]] <- density(samp_test$trees)
  dens_elevation_m2[[i]] <- density(samp_test$elevation)
  dens_dist_to_road_m2[[i]] <- density(samp_test$dist_to_road)
}

# ... Step 4. Compare observed and predicted distributions ----

# Now we compare our observed and predicted distributions. 

# Elevation -- same in both models
plot(dens_elevation_m1[[1]], col = "gray70", 
     ylim = c(0, 0.0025), main = "Elevation -- Correct Model")
for (i in 2:length(dens_elevation_m1)) {
  lines(dens_elevation_m1[[i]], col = "gray70") 
}

# Now add the actual used distribution we estimated in step 1
lines(density(test$elevation[test$case_]), col = "black")

# Now add the available
lines(density(test$elevation[!test$case_]), col = "red", lty = 2)

plot(dens_elevation_m2[[1]], col = "gray70", 
     ylim = c(0, 0.0025), main = "Elevation -- Incorrect Model")
for (i in 2:length(dens_elevation_m2)) {
  lines(dens_elevation_m2[[i]], col = "gray70") 
}

# Now add the actual used distribution we estimated in step 1
lines(density(test$elevation[test$case_]), col = "black")

# Now add the available
lines(density(test$elevation[!test$case_]), col = "red", lty = 2)

# Now trees -- the missing variable in second model
plot(dens_trees_m1[[1]], col = "gray70", 
     ylim = c(0, 0.03), main = "Trees -- Correct Model")
for (i in 2:length(dens_trees_m1)) {
  lines(dens_trees_m1[[i]], col = "gray70") 
}

# Now add the actual used distribution we estimated in step 1
lines(density(test$trees[test$case_]), col = "black")

# Now add the available
lines(density(test$trees[!test$case_]), col = "red", lty = 2)

plot(dens_trees_m2[[1]], col = "gray70", 
     ylim = c(0, 0.03), main = "Trees -- Incorrect Model")
for (i in 2:length(dens_trees_m2)) {
  lines(dens_trees_m2[[i]], col = "gray70") 
}

# Now add the actual used distribution we estimated in step 1
lines(density(test$trees[test$case_]), col = "black")

# Now add the available
lines(density(test$trees[!test$case_]), col = "red", lty = 2)

# Not extremely obvious, but the left-hand ascending arm is right
# on the edge of the available values.