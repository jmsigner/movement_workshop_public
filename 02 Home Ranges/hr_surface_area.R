# Calculating surface area within a home range from a DEM.

# If you have a digital elevation model (DEM -- a raster of elevation), 
# you can calculate the surface area within a home range. Here is a short
# demo using example data from a coyote in southern Utah, USA.

# Load packages ----
library(tidyverse)
library(amt)
library(sp)
library(sf)
library(raster)

# Load example data ----
# One coyote and two cougars
dat <- read_csv("data/coyote_cougar.csv") %>% 
  # Filter to just the coyote
  filter(id == "C028")

# Environmental data are in a multi-layered GeoTiff.
# Band 1 is elevation
dem <- raster("data/coyote_cougar_habitat.tif", band = 1)

# Fit home range ----
mcp <- dat %>% 
  make_track(x_, y_, t_, crs = 32612) %>% 
  hr_mcp()

# Plot
plot(dem)
plot(hr_isopleths(mcp)$geometry, add = TRUE)

# Calculate surface area ----
# We can use sp::surfaceArea()
?surfaceArea
# This function is expecting the data as a SpatialPixelsDataFrame
# We need to follow a couple of steps to get just the raster values
# from inside the MCP into that format.

# First, crop the raster to the approximate home range
# This isn't strictly necessary, but if you have a raster much 
# larger than your home range, it will save you unneccessary file size
# and (possibly) some computation time.
dem_mcp <- crop(dem, hr_isopleths(mcp))

# Now mask (set values outside of polygon to NA)
dem_mcp <- mask(dem_mcp, hr_isopleths(mcp))

# Plot
plot(dem_mcp)
plot(hr_isopleths(mcp)$geometry, add = TRUE, lwd = 3)

# Now convert to SpatialPixels object
dem_sp <- as(dem_mcp, "SpatialPixelsDataFrame")

plot(dem_sp)

# Now use sp::surfaceArea()
(sa <- surfaceArea(dem_sp))
# Units will always be units of the CRS (squared)
# So in the case of most projected CRS, units are m, and area is m^2

# Compare to area of flat MCP
data.frame(which = c("area", "surface area"),
           value = c(hr_area(mcp)$area, sa))
           