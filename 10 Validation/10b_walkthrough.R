#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#----------------Module 10 -- Validation---------------X
#----------------Last updated 2021-01-28---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

# Load packages ----
library(raster)
library(tidyverse)
library(amt)
library(sf)
library(rcompanion)
library(pROC)

# Simulate data ----
# We'll use a simplified version of the simulation we did in module 05.
# This time, we'll only use forage + temp + temp^2.

beta_forage = log(100)/1000
beta_temp2 = -1 * log(10)/36
beta_temp = beta_temp2 * -26

# Load the habitat layers
hab <- stack("05 HSF/geo/habitat.tif")
names(hab) <- c("forage", "temp", "predator", "cover")
plot(hab)

# Simulate locations
set.seed(20220128)
dat <- as.data.frame(hab, xy = TRUE) %>% 
  dplyr::select(x, y, forage, temp) %>% 
  # Calculate g(x)
  mutate(g = 
           # forage
           beta_forage * forage +
           # two terms for temperature
           beta_temp * temp +
           beta_temp2 * temp^2) %>% 
  # Calculate w(x)
  mutate(w = exp(g),
         w_prime = w/sum(w)) %>% 
  # Draw points
  mutate(lambda = 2000 * w_prime,
         n = rpois(n = nrow(.), lambda = lambda))

# Function to jitter data
jitter <- function(x, y, min = -25, max = 25) {
  res <- data.frame(x = x + runif(1, min, max),
                    y = y + runif(1, min, max))
  return(res)
}

# Now we split each row with n > 1 into a list element
dat_list <- dat %>% 
  filter(n > 0) %>% 
  split(1:nrow(.))

# And now we can create jittered points for each element of our list.
# We will use 'bind_rows()' (twice) to return a single data.frame
set.seed(654321)
gps <- lapply(dat_list, function(d) {
  replicate(d$n, jitter(x = d$x, y = d$y), simplify = FALSE) %>% 
    bind_rows()
}) %>% 
  bind_rows()

# This is what the result looks like.
head(gps)

# Split into training and testing ----
# We want to do some out-of-sample validation, so let's withhold ~20% of our
# data for testing our model, using the remaining 80% to fit the model.
set.seed(20220128)
gps$train <- rbinom(n = nrow(gps), size = 1, prob = 0.8)

train <- gps %>% filter(train == 1)
test <- gps %>% filter(train == 0)

nrow(train)/nrow(gps)
nrow(test)/nrow(gps)

# Let's format both our training and testing data for an HSF.
# We want to generate available points across our entire raster.
# We'll use a polygon for that.

r_poly <- st_bbox(hab) %>% 
  st_as_sfc() %>% 
  st_sf()

set.seed(123456)
# Random points for training data
r_train <- random_points(r_poly, n = nrow(train) * 100)
# Random points for testing data
r_test <- random_points(r_poly, n = nrow(test) * 100)

# Format training
train <- train %>% 
  make_track(x, y, crs = 32612) %>% 
  mutate(case_ = TRUE) %>% 
  bind_rows(r_train) %>% 
  extract_covariates(hab) %>% 
  dplyr::select(-predator, -cover) %>% 
  # Assign large weights to available points
  mutate(weight = ifelse(case_, 1, 1e5))

# Format testing
test <- test %>% 
  make_track(x, y, crs = 32612) %>% 
  mutate(case_ = TRUE) %>% 
  bind_rows(r_test) %>% 
  extract_covariates(hab) %>% 
  dplyr::select(-predator, -cover) %>% 
  # Assign large weights to available points
  mutate(weight = ifelse(case_, 1, 1e5))

# Fit model ----
# Let's fit two models. One that is correctly specified, and one that is
# missing the quadratic term for temperature.

# Correct
m1 <- glm(case_ ~ forage + temp + I(temp^2),
          family = binomial(), weights = weight, data = train)
# Incorrect
m2 <- glm(case_ ~ forage + temp,
          family = binomial(), weights = weight, data = train)


# We're ready to evaluate our model

# Model evaluation ----
# ... Pseudo R-squared ----

# Before we jump into measures of external validity, let's take a 
# quick look at a measure of internal validity, i.e., how well does
# our model fit to the data we actually used to fit it.

?rcompanion::nagelkerke

(m1r2 <- nagelkerke(m1)) # Quite low
(m2r2 <- nagelkerke(m2)) # Slightly lower

cbind(m1r2$Pseudo.R.squared.for.model.vs.null,
      m2r2$Pseudo.R.squared.for.model.vs.null)

# We barely see a difference.

# ... AUC ----

# One method that's commonly used to assess models for binary responses 
# is the Area Under the Curve metric, or AUC. It measures the area under the 
# receiver-operator curve (ROC), which is a good measure of the trade-off 
# between sensitivity and specificity in a binary predictor.

# Note that you can calculate AUC as a measure of internal validity (without
# a testing set) or as a measure of external validity (on a testing set).

# We'll use the package `pROC` here, but there are *many* R packages that
# can calculate this statistic for you.

# We will use the function 'roc()' to produce the ROC. We can then plot it
# or calculate the area under it.

# We need to provide 'roc()' two arguments:
#   1. 'response' the binary data used to fit the model
#   2. 'predictor' the predicted probability under the model
#       In our case this is the linear prediction from our model,
#       backtransformed using the inverse-logit.

m1_roc <- roc(response = as.numeric(test$case_),
                  predictor = predict(m1, 
                                      newdata = test,
                                      type = "response"))

# Now we can plot the curve
plot(m1_roc)

# We can also calculate AUC
auc(m1_roc)

# Not bad!

m2_roc <- roc(response = as.numeric(test$case_),
              predictor = predict(m2, 
                                  newdata = test,
                                  type = "response"))
auc(m2_roc) # Only slightly lower

# ... Boyce method ----
# The general approach with binning is to aggregate the testing data into
# bins of:
#   1. equal area (geographic space)
#   2. equal count (geographic space)
#   3. equal interval (environmental space)

# Here, we'll use equal-interval bins.

# We also need to know which cell of the landscape each point falls in
# so that we can adjust our expected number of points by the area covered.

# Let's also break our habitat variables into 5 bins. Note that cover is
# already in discrete classes, so we don't need to bin.

# The base R function 'cut()' can bin our data for us, but the tidyverse
# function 'cut_interval()' will make it easier to control the details
?cut_interval

# Let's bin and summarize:
test_bins <- test %>%
  # Keep only the observed points
  filter(case_) %>% 
  # Assign cell number
  mutate(cell = cellFromXY(hab, as.matrix(.[, c("x_", "y_")]))) %>% 
  # Bin variables
  mutate(forage_bin = cut_interval(forage, n = 5),
         temp_bin = cut_interval(temp, n = 5)) %>% 
  # Now group by bins
  dplyr::group_by(forage_bin, temp_bin) %>% 
  # Summarize observations
  summarize(cells = n_distinct(cell),
            obs = sum(case_),
            obs_per_cell = obs/cells) %>% 
  # Ungroup
  ungroup()

# What did we get?
test_bins

# You can see we have the number of observed points by habitat 
# bin, as well as the density of points. Now how can we compare this to 
# our model?

# The basic idea is that we want to predict the value of each habitat and
# elevation value using our HSF, and then see how strong the correlation
# is between the HSF and the density of observed points.

# Let's convert our text label for each bin into the mean value
# for that bin. Since we don't have nice, round numbers from 'cut_interval()', 
# we need some string manipulation.

# Here's a function that will do it for us.
get_mean <- function(x) {
  # Get rid of parentheses
  x <- gsub(pattern = "(", replacement = "", x, fixed = TRUE)
  x <- gsub(pattern = ")", replacement = "", x, fixed = TRUE)
  # Get rid of square brackets
  x <- gsub(pattern = "[", replacement = "", x, fixed = TRUE)
  x <- gsub(pattern = "]", replacement = "", x, fixed = TRUE)
  # Split by comma
  y <- strsplit(x, ",")
  # Average
  z <- unlist(lapply(y, function(yy) {
    mean(as.numeric(yy))
  }))
  # Return
  return(z)
}

# Example
levels(test_bins$forage_bin)
get_mean(levels(test_bins$forage_bin))

test_bins <- test_bins %>% 
  mutate(forage = get_mean(forage_bin),
         temp = get_mean(temp_bin))

# Now that we have our habitat variables, we can use 'predict()' to
# calculate w(x). Remember, we can get the linear prediction (g(x)), 
# subtract the intercept, and exponentiate to get w(x).

test_bins <- test_bins %>% 
  # Linear predictor
  mutate(g = predict(m1, newdata = ., type = "link")) %>% 
  # Subtract intercept and exp()
  mutate(w = exp(g - coef(m1)[1]))

# Done! Now we can evaluate our model by:
#   1. plotting
#   2. calculating the correlation

# Plot
ggplot(test_bins, aes(x = w, y = obs_per_cell)) +
  geom_point() +
  geom_line() +
  theme_bw()

# Correlation
cor(test_bins$w, test_bins$obs_per_cell, method = "spearman")

# The Boyce method took us significantly more coding than the pR2 or
# AUC approaches, but hopefully you can see much clearer ties to our
# inhomogeneous Poisson point process model here.

# ... UHC plots ----

# Package `uhcplots` cannot be (easily) installed due to a missing
# dependency (`SDMTools`). You can install an archived version of this
# package to resolve the issue, but we can simply recreate the UHC plots
# ourselves.

# Fieberg et al. broke the approach down into 4 steps.
#   1. Summarize used and available distributions in test data
#   2. Fit candidate model (we already did)
#   3. Bootstrap
#   4. Compare observed and predicted

# ... Step 1. Summarize distribution of covariates ----


# But what we really want is one distribution for the used points, and 
# one distribution for the available points (for the *test* dataset)
# for each variable.
test_forage_used <- density(test$forage[test$case_])
test_forage_avail <- density(test$forage[!test$case_])

plot(test_forage_used,
     main = "Forage Density")
lines(test_forage_avail,
      col = "red", lty = 2)

# We can also do this with ggplot
ggplot(test, aes(x = forage, color = case_)) +
  geom_density() +
  theme_bw()

# We can do them all at once with ggplot, too.
test %>% 
  pivot_longer(forage:temp) %>% 
  ggplot(aes(x = value, color = case_)) +
  facet_wrap(~ name, scales = "free") +
  geom_density() +
  theme_bw()

# ... Step 2. Fit a model to the training data ----

# We already did this, but now we also need the betas from our model and 
# the variance-covariance matrix
m1_betas <- coef(m1)
m1_vcov <- vcov(m1)

m2_betas <- coef(m2)
m2_vcov <- vcov(m2)

# ... Step 3. Create predicted distribution using fitted model ----

# We are going to create the predicted distribution by repeating the
# substeps many times. Let's try it with 500 iterations (you probably
# want more for your own data).

M <- 500

# Initialize blank lists
dens_forage_m1 <- list()
dens_temp_m1 <- list()

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
  samp_test <- test %>% 
    # Calculate w(x), i.e., the exponential HSF
    mutate(w = exp(samp_beta[2] * forage +
                     samp_beta[3] * temp +
                     samp_beta[4] * temp^2)) %>% 
    # Select rows in proportion to w(x)
    # How many? One for each used point.
    slice_sample(n = sum(test$case_), weight_by = w)
  
  # ... ... c. summarize covariates at these locations ----
  dens_forage_m1[[i]] <- density(samp_test$forage)
  dens_temp_m1[[i]] <- density(samp_test$temp)
}

# Same for m2
# Initialize blank lists
dens_forage_m2 <- list()
dens_temp_m2 <- list()

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
  samp_test <- test %>% 
    # Calculate w(x), i.e., the exponential HSF
    mutate(w = exp(samp_beta[2] * forage +
                     samp_beta[3] * temp)) %>% 
    # Select rows in proportion to w(x)
    # How many? One for each used point.
    slice_sample(n = sum(test$case_), weight_by = w)
  
  # ... ... c. summarize covariates at these locations ----
  dens_forage_m2[[i]] <- density(samp_test$forage)
  dens_temp_m2[[i]] <- density(samp_test$temp)
}

# ... Step 4. Compare observed and predicted distributions ----

# Now we compare our observed and predicted distributions. 

# The easiest way to plot all of the available densities is with a for() loop.
# We'll plot the first one, and then use a loop to plot numbers 2 -- 500.

plot(dens_forage_m1[[1]], col = "gray70", 
     ylim = c(0, 0.003), main = "Forage -- Correct Model")
for (i in 2:length(dens_forage_m1)) {
  lines(dens_forage_m1[[i]], col = "gray70") 
}

# Now add the actual used distribution we estimated in step 1
lines(density(test$forage[test$case_]), col = "black")

# Now add the available
lines(density(test$forage[!test$case_]), col = "red", lty = 2)

plot(dens_forage_m2[[1]], col = "gray70", 
     ylim = c(0, 0.003), main = "Forage -- Incorrect Model")
for (i in 2:length(dens_forage_m2)) {
  lines(dens_forage_m2[[i]], col = "gray70") 
}

# Now add the actual used distribution we estimated in step 1
lines(density(test$forage[test$case_]), col = "black")

# Now add the available
lines(density(test$forage[!test$case_]), col = "red", lty = 2)

# Now temperature -- the variable that is misspecified
plot(dens_temp_m1[[1]], col = "gray70", 
     ylim = c(0, 0.25), main = "Temperature -- Correct Model")
for (i in 2:length(dens_temp_m1)) {
  lines(dens_temp_m1[[i]], col = "gray70") 
}

# Now add the actual used distribution we estimated in step 1
lines(density(test$temp[test$case_]), col = "black")

# Now add the available
lines(density(test$temp[!test$case_]), col = "red", lty = 2)

plot(dens_temp_m2[[1]], col = "gray70", 
     ylim = c(0, 0.25), main = "Temperature -- Incorrect Model")
for (i in 2:length(dens_temp_m2)) {
  lines(dens_temp_m2[[i]], col = "gray70") 
}

# Now add the actual used distribution we estimated in step 1
lines(density(test$temp[test$case_]), col = "black")

# Now add the available
lines(density(test$temp[!test$case_]), col = "red", lty = 2)


