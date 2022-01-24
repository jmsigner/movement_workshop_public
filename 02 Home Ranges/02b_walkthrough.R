#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#---------------Module 02 -- Home Ranges---------------X
#----------------Last updated 2021-01-06---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

# Load packages ----
library(tidyverse)
library(raster)
library(amt)

# Generate example data ----
# We need to generate some location data with home-ranging behavior for
# analysis here. We will use a biased, correlated random walk (in discrete
# time) to accomplish this. The bias is toward a home range centroid, which
# we will define.

# Source function for BCRW
source("fun/bcrw.R")

# We'll generate data for 4 different individuals in the Bear River mountains 
# just east of Logan, UT. We'll assume our spatial coordinates are in 
# UTMs (zone 12).

centroids <- list("A01" = c("x" = 446589, "y" = 4625899),
                  "A02" = c("x" = 440284, "y" = 4618197),
                  "A03" = c("x" = 448796, "y" = 4613795),
                  "A04" = c("x" = 442367, "y" = 4627042)
)

# Set random seed
set.seed(20220124)

# Generate random walk
dat <- lapply(centroids, function(cent) {
  # The basic BCRW
  x <- bcrw(start_loc = cent + rnorm(n = 2, mean = 0, sd = 500),
            centroid = cent,
            n_steps = 1000,
            sl_distr = c('shape' = 1, 'scale' = 300),
            rho = 0.25, 
            beta = 0.1) 
  
  # Assign actual dates and times (start on June 20, 2021)
  x$date <- as.POSIXct("2021-06-20 00:00:00") +
    x$t * 60 * 60
  
  # Return
  return(x)
}) %>% 
  bind_rows(.id = "ID")

# See what we got
head(dat)
tail(dat)

# Format as 'track_xyt' for `amt`
dat <- dat %>% 
  make_track(x, y, date, ID = ID, crs = 32612)

# Plot locations
ggplot(dat, aes(x = x_, y = y_, color = ID, group = ID)) +
  geom_point() +
  coord_equal() +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw()

# Subset to just one individual
a01 <- dat %>% 
  filter(ID == "A01")

# Home range functions ----
# The home range constructor functions in `amt` start with 'hr_*()'.
# E.g., to create an MCP, use 'hr_mcp()'.
# You can find them all here:
?hr_mcp

# We'll see how the functions work by starting with one individual.

# ... Geometric estimators ----
# We'll begin with the geometric estimators.

# ... ... MCPs ----
# MCPs are the simplest home range estimator.
# They only require one user decision: the level.

# We can construct MCPs with multiple levels with a single call to the 
# function 'hr_mcp()'.

mcps <- hr_mcp(a01, levels = c(0.5, 0.75, 0.95, 1))

# Examine the object we created
class(mcps)
str(mcps, max.level = 1)

# The object 'mcps' is, at its most basic, a list. The list has elements:
#   * 'mcp' -- simple features (sf) data.frame with each MCP
#   * 'levels' -- vector of the levels we requested
#   * 'estimator' -- character string giving the type of estimator ("mcp")
#   * 'crs' -- object of class "crs" giving the coordinate reference system
#   * 'data' -- data.frame of class "track_xyt" with the original input data

# Other 'hr' objects created by `amt` have a similar structure.

# Each 'hr' object created by `amt` also has default "methods" (functions)
#   * 'hr_area()' -- calculates the home range area
#   * 'hr_isopleth()' -- returns the home range isopleths at specified levels
#   * 'hr_overlap()' -- calculates the overlap between two home ranges
#   * 'plot()' -- plots the home range

# Let's see how the functions work with our MCPs.

# Area
hr_area(mcps) # returns a tibble
hr_area(mcps, units = TRUE) # returns area as a 'units' object for easy conversion
# Get area in km^2
hr_area(mcps, units = TRUE) %>% 
  mutate(area = units::set_units(area, "km^2"))

# Isopleths
hr_isopleths(mcps) 
# Returns sf object -- useful, e.g., for export to shapefile
?st_write

# Overlap (obviously overlap will be 1 -- more interesting example later)
hr_overlap(mcps, mcps)

# Plot
plot(mcps) # Plots all isopleths and adds data on top

# Custom plot with ggplot2
# Since we have access to the sf objects for each polygon, we don't have
# to rely on the default plotting method if we don't want to.

hr_isopleths(mcps) %>% 
  # Make level a factor for discrete color scales
  # Can control order and labels here
  mutate(level = factor(level, 
                        levels = c("1", "0.95", "0.75", "0.5"),
                        labels = c("100%", "95%", "75%", "50%"))) %>% 
  ggplot() +
  geom_sf(aes(color = level), 
          fill = NA, size = 2) +
  geom_point(data = mcps$data, aes(x = x_, y = y_),
             color = "gray30", alpha = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  scale_color_brewer(name = "MCP Level",
                     palette = "Set1") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(colour = "black", size = 0.8))

# ... ... LoCoHs ----
# LoCoHs are created by fitting MCPs to smaller subsets of points, defined
# by a neighborhood rule, which are then merged into larger polygons.
# Unlike MCPs, the resulting polygons are not necessarily convex and
# can also have holes.

# The LoCoH algorithms require more user decisions. Like with all estimators,
# the user must choose the isopleth levels. Users also must choose the
# neighborhood rule ("k", "r", or "a") and the parameter for the neighborhood
# algorithm.

# LoCoHs can approximate a utilization distribution if many isopleths
# are drawn.

# Here, we will use the "k" nearest neighbor algorithm, the simplest choice.
# We will create hulls with the 20 nearest neighbors (k = 21). We'll ask the
# algorithm to fit an isopleth for every 10th percentile.

locohs <- hr_locoh(a01, type = "k", n = 21, levels = seq(0.1, 1, by = 0.1))

# Remember, the structure is similar between hr objects.
class(locohs)
str(locohs, max.level = 1)

# And remember that we have the four methods available
hr_area(locohs)
hr_isopleths(locohs)
# Skipping hr_overlap() for now
plot(locohs)

# It is hard to tell what's going on with the LoCoH without some shading, since
# there are multiple polygons and even holes. Let's make a custom plot.
(locoh_plot <- hr_isopleths(locohs) %>% 
    mutate(level = fct_rev(factor(level))) %>% 
    ggplot() +
    geom_sf(aes(fill = level,
                color = level), 
            size = 1) +
    xlab(NULL) +
    ylab(NULL) +
    scale_color_manual(name = "Isopleth",
                       values = rev(gray.colors(10))) +
    scale_fill_manual(name = "Isopleth",
                      values = rev(gray.colors(10))) +
    theme_bw() +
    theme(legend.position = "bottom", 
          legend.box.background = element_rect(colour = "black", size = 0.8)))

# We can add the points, although they tend to obscure the smaller isopleths
locoh_plot +
  geom_point(data = locohs$data, aes(x = x_, y = y_),
             color = "red", size = 0.5, alpha = 0.25)

# ... Probabilistic Estimators ----

# ... ... kernel density estimation in general ----
# The simplest probabilistic estimator is the KDE. It assumes all locations
# are independent, which is probably a bad assumption for frequent GPS fixes.

# Kernel density estimation is a solution to a general problem in statistics:
# how do we estimate a smooth probability density function from a sample of
# data? It is a problem for continuous distributions and is analogous to 
# choosing the bin size for a histogram.

hist(a01$x_, breaks = 5) # too smooth
hist(a01$x_, breaks = 20) # intermediate
hist(a01$x_, breaks = 50) # too choppy

#As a general statistical problem, R has a built-in function for this.
plot(density(a01$x_)) # automatically selects bandwidth

# Compare to the histogram with 20 bins
hist(a01$x_, breaks = 20, freq = FALSE)
lines(density(a01$x_), col = "red", lwd = 2)

# If we choose a smaller or larger bandwidth, our results change.
plot(density(a01$x_, bw = 10), main = "Small Bandwidth -- Too Choppy")
plot(density(a01$x_, bw = 1000), main = "Large Bandwidth -- Too Smooth")
plot(density(a01$x_), main = "Automatic Bandwidth")

# ... ... the KDE home range ----

# Fitting a KDE to telemetry data is the same problem, but in 2 dimensions.
# Choosing a bandwidth is important to the final result, and there are multiple
# possible methods for doing so.

# Methods currently implemented in `amt` are: 
#   * 'hr_kde_ref()' -- the reference bandwidth (appropriate for unimodal UDs).
#   * 'hr_kde_ref_scaled()' -- a scaled version of the reference bandwidth that
#                             finds the smallest bandwidth that produces a
#                             home range isopleth consisting of n (usually 1)
#                             polygons. Suggested by Kie 2013.
#   * 'hr_kde_pi()' -- the plug-in bandwidth. Suggested by Gitzen et al. 2006.
#   * 'hr_kde_lscv()' -- the least squares cross validation bandwidth. Suggested
#                       by Seaman & Powell 1996.

# We'll demonstrate the reference bandwidth here, which is the default.

# Fit KDEs
kdes <- hr_kde(a01, levels = c(0.5, 0.95))

# Examine our object
class(kdes)

# The structure for probabilistic estimators is similar to the geometric
# estimators, but with some important additions.
str(kdes, max.level = 1)

# Two important additions are two 'RasterLayer' objects:
#   * 'ud' -- the rasterized utilization distribution
#   * 'trast' -- the template raster used to rasterize the UD

# We can manipulate these using the `raster` package, including exporting them
# as, e.g., GeoTiff files.
?raster::writeRaster

# All 4 methods are available for our probabilistic estimators, too.
hr_area(kdes)
hr_isopleths(kdes)
# Skipping hr_overlap()
plot(kdes)

# Here is a useful trick for plotting rasters with ggplot2:
# There is a method for 'as.data.frame()' for a `Raster*` object, which
# has an argument called 'xy'. Passing 'xy = TRUE' appends the coordinates
# to the data, allowing you to use the data.frame with ggplot.
ud_df <- as.data.frame(kdes$ud, xy = TRUE)
head(ud_df)

ggplot(ud_df, aes(x = x, y = y, fill = layer)) +
  geom_raster()

# Let's make a plot with our 50% and 95% isopleths drawn over the UD raster.

ggplot() +
  geom_raster(data = ud_df, aes(x = x, y = y, fill = layer)) +
  geom_sf(data = hr_isopleths(kdes), aes(color = fct_rev(factor(level))),
          fill = NA, size = 1) +
  scale_fill_viridis_c(name = "Probability\nDensity (UD)") +
  scale_color_brewer(name = "Isopleth", palette = "Set2") +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw()

# ... ... aKDE ----

# The `amt` function 'hr_akde()' provides a very convenient interface
# to the `ctmm` package. aKDEs are fit in `amt` via a wrapper to the `ctmm` 
# functions. Many more details are available here:
?ctmm::akde

# The aKDE relies on a continuous time movement model (more on those in the 
# next module). `amt` also provides a convenient wrapper for fitting those.
?fit_ctmm

# Without going into detail on the different CTMMs, we'll demonstrate fitting
# an aKDE with an Ornstein-Uhlenbeck model.
akdes <- hr_akde(a01, model = fit_ctmm(a01, "ou"), levels = 0.95)

# Examine the object
class(akdes)

# Structure
str(akdes, max.level = 1)

# Notice that the 'ud' and 'trast' elements are still available, as with the KDE
# Notice that we also have the 'akde' and 'model' objects from `ctmm`

# Our methods still apply
hr_area(akdes) # notice we now have confidence intervals around our estimates
hr_isopleths(akdes)
# Skipping hr_overlap()
plot(akdes) # The 3 lines represent the estimate and the confidence bounds

# We have more plotting options with the individual parts
plot(akdes$akde) # plotting method from `ctmm`

# And we can make a custom plot with ggplot by grabbing the pieces
hr_isopleths(akdes) %>% 
  mutate(type = case_when(
    what == "estimate" ~ "Estimate",
    TRUE ~ "Conf. Int.")) %>% 
  ggplot() +
  geom_raster(data = as.data.frame(akdes$ud, xy = TRUE),
              aes(x = x, y = y, fill = layer)) +
  geom_point(data = akdes$data, aes(x = x_, y = y_),
             color = "gray30", alpha = 0.5) +
  geom_sf(aes(linetype = type, size = type), color = "blue", fill = NA) +
  xlab(NULL) +
  ylab(NULL) +
  scale_linetype_manual(name = "aKDE",
                        breaks = c("Estimate", "Conf. Int."),
                        values = c("solid", "dashed")) +
  scale_size_manual(name = "aKDE",
                    breaks = c("Estimate", "Conf. Int."),
                    values = c(1.75, 0.75)) +
  scale_fill_viridis_c(name = "aKDE UD", option = "A") +
  theme_bw()

# Multiple animals ----
# All of the home range functions in `amt` are designed to work on one
# individual. But what if we have multiple individuals?

# One option is to break the 'track_xyt' data.frame into a list by individual.
dat_list <- split(dat, dat$ID)

# Now each element of 'dat_list' is an individual, and we can use lapply()
# to fit home ranges, e.g.:
mcp_list <- lapply(dat_list, hr_mcp, levels = c(0.5, 0.95))
lapply(mcp_list, plot)

# Lists are very valuable data structures in R, but we're often using tools
# designed to work on data.frames. We can combine the intuitive structure of
# data.frames with the power of lists in a "nested data.frame".

dat_nest <- dat %>% 
  nest(track = x_:t_)
dat_nest

# The column "track" is now a list of 'track_xyt' objects, while the unique
# attribute columns for each individual are still in the data.frame format.

# We can lapply across the track column, e.g.:
lapply(dat_nest$track, hr_mcp, levels = c(0.5, 0.95))

# Or, we can use the tidyverse equivalent, map()
dat_nest <- dat_nest %>% 
  mutate(mcp = map(track, hr_mcp, levels = c(0.5, 0.95)))
dat_nest

# The column "mcp" is a list of the fitted MCPs.
plot(dat_nest$mcp[[1]])
map(dat_nest$mcp, plot)

# And finally, a demonstration of hr_overlap()
ggplot() +
  geom_sf(data = dat_nest$mcp[[1]]$mcp, size = 1, 
          color = "orange", fill = NA) +
  geom_sf(data = dat_nest$mcp[[4]]$mcp, size = 1, 
          color = "blue", fill = NA) +
  theme_bw()

hr_overlap(dat_nest$mcp[[1]], dat_nest$mcp[[4]])
hr_overlap(dat_nest$mcp[[4]], dat_nest$mcp[[1]])

# 'hr_overlap()' also works with a list
ov <- hr_overlap(dat_nest$mcp, which = "all")
head(ov, 12)
tail(ov, 12)

# We will see more of the power of nested data.frames in future modules.
