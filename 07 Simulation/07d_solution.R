library(NLMR)
library(raster)

source("fun/dispersal_kernel.R")

# 1. Simulate a binary landscape encoded 0 and 1. This could be for example a matrix of forest and non-forest. 
# Landscape

set.seed(123)
lscp <- nlm_gaussianfield(200, 200)
raster::plot(lscp)

lscp <- lscp > 0.5
plot(lscp)

lscp <- stack(lscp)

# 2. Create three different dispersal kernels: 

# i. The animal moves with constant step-lengths that follows a gamma distribution with shape 2 and scale 5. 
curve(dgamma(x, scale = 5, shape = 2), from = 0, to = 50, 
      ylab = "Density", xlab = "sl_")

# Get the coefficients for the movement kernel
scale_to_sl(5)
shape_to_log_sl(2)

# Create dispersal kernel
# Note we have no habitat selection here
dk1 <- dispersal_kernel(
  ~ sl_ + log_sl_, 
  coefficients = c("sl_" = -0.2, "log_sl_" = 1), 
  spatial.covars = stack(lscp), start = c(50, 50), 
  max.dist = 50, return.raster = TRUE)

# Plot the results
plot(dk1)

# ii. The animal moves with the same step-length distribution as in i) but now shows a preference for habitat 1.

# Now we add a selection for layer (note the use of `layer_end` here)
lscp
dk2 <- dispersal_kernel(
  ~ sl_ + log_sl_ + layer_end, 
  coefficients = c("sl_" = -0.2, "log_sl_" = 1, "layer_end" = 1), 
  spatial.covars = stack(lscp), start = c(100, 100), 
  max.dist = 50, return.raster = TRUE)
plot(dk2)

# iii. The animal moves with the same step-length distribution as in i) but now
# shows a preference for habitat 1 and high directional persistence (kappa = 5).

# Now we have to include the cos_ta in the issf as it controls for directional
# persistence.

kappa_to_cos_ta(5)
dk3 <- dispersal_kernel(
  ~ sl_ + log_sl_ + layer_end + cos_ta_, 
  coefficients = c("sl_" = -0.2, "log_sl_" = 1, "layer_end" = 1, "cos_ta_" = 5), 
  spatial.covars = stack(lscp), start = c(100, 100), 
  max.dist = 50, return.raster = TRUE)
plot(dk3)

# Note, we could also set the direction
dk3 <- dispersal_kernel(
  ~ sl_ + log_sl_ + layer_end + cos_ta_, 
  coefficients = c("sl_" = -0.2, "log_sl_" = 1, "layer_end" = 1, "cos_ta_" = 5), 
  spatial.covars = stack(lscp), start = c(100, 100), direction = pi,
  max.dist = 50, return.raster = TRUE)
plot(dk3)


dk3 <- dispersal_kernel(
  ~ sl_ + log_sl_ + layer_end + cos_ta_, 
  coefficients = c("sl_" = -0.2, "log_sl_" = 1, "layer_end" = 3, "cos_ta_" = 5), 
  spatial.covars = stack(lscp), start = c(100, 100), direction = pi/2,
  max.dist = 50, return.raster = TRUE)
plot(dk3)
