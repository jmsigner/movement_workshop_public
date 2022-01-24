# Habitat rasters for walkthrough 05 -- HSFs

# Load packages----
library(tidyverse)
library(raster)
library(mvtnorm)

# Coordinates and distance matrix ----
coords <- expand.grid("x" = seq(-1000, 950, by = 50),
                      "y" = seq(-1000, 950, by = 50)) %>% 
  arrange(x, y) %>% 
  as.matrix()

dm <- as.matrix(dist(coords))
dm <- dm/max(dm)

# Covariance matrix ----
# Exponential covariance function
expcov <- function(dists, rho, sigma) {
  n <- dim(dists)[1]
  result <- matrix(nrow = n, ncol = n)
  sigma2 <- sigma^2
  result <- sigma2 * exp(-1 * dists/rho)
  return(result)
}


# Generate habitats ----

# ... forage ----
# units are g/m2
# reasonable values are 0 - 1000
# should probably be 0 inflated (rocky outcrops?)

set.seed(20220105 + 1)
forage <- rmvnorm(n = 1,
                  mean = rep(4, nrow(coords)),
                  sigma = expcov(dm, 3, 2))
# Stretch to range (create some negatives)
forage <- (forage - (min(forage) * 1.08)) / max(forage - (min(forage) * 1.08)) * 1000

# Cut any negative values to 0
forage[forage < 0] <- 0

# Check
hist(forage, breaks = 30)

# ... mean annual temp ----
# units are degrees C
# reasonable values are 0 - 20 deg C
set.seed(20220105 + 2)
temp <- rmvnorm(n = 1,
                mean = rep(13, nrow(coords)),
                sigma = expcov(dm, 1, 5))

# Check
hist(temp, breaks = 30)

# ... predation risk ----
# units are predators/100 km2
# reasonable values are 0 - 15 (think large predator like wolf)
# want mean to vary through space

set.seed(20220105 + 3)
pred <- rmvnorm(n = 1,
                mean = ((sin(1:nrow(coords)/100) + 1) * 5) + 10,
                sigma = expcov(dm, 5, 6))
pred[pred < 0] <- 0

# Check
hist(pred, breaks = 30)

# ... landcover ----
# categories are:
#   - grassland (50%)
#   - forest (30%)
#   - wetland (20%)
# generate with MV norm and then take ceiling to get integer and categorize
# mean should vary to give desire ratios

# Smooth mean with KD smooth
lc_kd <- density(c((rep(1, round(0.5 * nrow(coords)))),
                  rep(2, round(0.3 * nrow(coords))),
                  rep(3, round(0.2 * nrow(coords)))), 
                bw = 0.3)

set.seed(20220105 + 4)
lc_m <- sort(sample(lc_kd$x, size = 1600, prob = lc_kd$y, replace = TRUE))

lc_c <- rmvnorm(n = 1,
                mean = lc_m,
                sigma = expcov(dm, 2, 0.4))
lc <- round(lc_c)
lc[lc > 3] <- 3
lc[lc < 1] <- 1

# Check
table(lc)/nrow(coords)

# Note: add a patch of wetland in SW and NW corner after compiling 
# into data.frame

# Compile in data.frame
g_dat <- data.frame(x = coords[, 1] + 447000,
                    y = coords[, 2] + 4626000,
                    forage = forage[1, ],
                    temp = temp[1, ],
                    pred = pred[1, ],
                    lc = lc[1, ])

g_dat <- g_dat %>% 
  # SW corner
  mutate(lc = case_when(
    x <= 446250 & y <= 4625250 ~ 3,
    TRUE ~ lc)) %>% 
  mutate(lc = case_when(
    x <= 446150 & y > 4626825 ~ 3,
    TRUE ~ lc))

# Convert to rasters and save ----
rast <- rasterFromXYZ(g_dat, res = 50, crs = 32612)

# Check
plot(rast)

writeRaster(rast, "05 HSF/geo/habitat.tif", overwrite = TRUE)
