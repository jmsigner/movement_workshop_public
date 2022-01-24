#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#-------------------Module 03 -- CTMM------------------X
#----------------Last updated 2021-01-22---------------X
#-------------------Exercise Solution------------------X
#######################################################X

# Using elephant data

# Load packages ----
library(tidyverse)
library(lubridate)
library(amt)
library(sf)
library(crawl)

# 1. Load and subset data ----
# Load elephant data and only keep individual Salif Keita in 2009
dat <- read_csv("data/elephants.csv") %>% 
  filter(id == "Salif Keita" & year(timestamp) == 2009) %>% 
  # Subset to one point a day
  filter(hour(timestamp) == 0) %>% 
  # Transform the track to a projected CRS. 
  # Using 32630 for UTM zone 30 on the WGS84 spheroid (https://epsg.io/32630).
  make_track(long, lat, timestamp, temp = temperature, crs = 4326) %>% 
  transform_coords(32630)

head(dat)

# Subset by removing entire months of May and June
obs <- dat %>% 
  filter(!month(t_) %in% 6)

plot(obs$t_)

plot(dat, pch = 16)
points(obs, pch = 16, col = "red")

# 2. Calculate statistic for data ----
# Let's say our statistic of interest is the area of a 90% KDE.

(kde_true <- dat %>% 
  hr_kde(levels = 0.9) %>% 
  hr_area(units = TRUE) %>% 
  mutate(area = units::set_units(area, "km^2")))

(kde_obs <- obs %>% 
  hr_kde(levels = 0.9) %>% 
  hr_area(units = TRUE) %>% 
  mutate(area = units::set_units(area, "km^2")))

# 3. Convert to sf ----
sf_true <- dat %>% 
  st_as_sf(coords = c("x_", "y_"), crs = 32630)

sf_obs <- obs %>% 
  st_as_sf(coords = c("x_", "y_"), crs = 32630)

# 4. Fit CRW model ----
# Set seed for good initial values
set.seed(3)
mod <- crwMLE(data = sf_obs, Time.name = "t_", time.scale = "days")

# 5. Multiple imputation ----
# Function to simulate and calculate KDE area
sim_iter <- function(m) {
  # Setup simulation
  sim <- crwSimulator(m, predTime = "1 day")
  # Get a posterior sample
  sim_tracks <- crwPostIS(sim, fullPost = FALSE)
  # Convert to sf points
  pts <- crw_as_sf(sim_tracks, ftype = "POINT")
  # Convert to sf lines
  lns <- crw_as_sf(sim_tracks, ftype = "LINESTRING")
  # Convert to track_xy and fit 90% KDE
  kde <- st_coordinates(pts) %>% 
    make_track(X, Y, crs = 5070) %>% 
    hr_kde(levels = 0.9) %>% 
    hr_area(units = TRUE) %>% 
    mutate(area = units::set_units(area, "km2")) %>% 
    mutate(dataset = "Multiple Imputation")
  
  # List to return
  res <- list(line = lns,
              kde = kde)
  return(res)
}

# Run imputation for 50 iterations 
#   (for speed, really need many more for inference)
# Takes about 80 sec
set.seed(123)
system.time({
  mult_imp <- replicate(50, sim_iter(m = mod), simplify = FALSE)
})

# Pull out KDE areas and lines
l <- do.call(rbind, lapply(mult_imp, getElement, "line"))

kdes <- do.call(rbind, lapply(mult_imp, getElement, "kde"))

# Summary stats
(kde_summ <- kdes %>% 
    group_by(level, what) %>% 
    summarize(mean = mean(area),
              lwr = quantile(area, 0.025),
              upr = quantile(area, 0.975)) %>% 
    rename(area = mean))

# 6. Figures ----
# ... statistic ----
bind_rows("True" = kde_true, 
          "Observed" = kde_obs, 
          "Imputed" = kde_summ,
          .id = "set") %>% 
  mutate(across(area:upr, as.numeric)) %>% 
  ggplot(aes(x = set, y = area, ymin = lwr, ymax = upr)) +
  geom_point(size = 3) +
  geom_errorbar(width = 0.2) +
  xlab("Dataset") +
  ylab(expression("90% KDE Area " * (km^2))) +
  theme_bw()

# ... map ----
ggplot() +
  geom_sf(data = l, aes(color = "Imputed"), size = 0.9, show.legend = TRUE) +
  geom_sf(data = sf_true, aes(color = "True")) +
  geom_sf(data = sf_obs, aes(color = "Observed")) +
  scale_color_manual(name = "Dataset",
                     breaks = c("Imputed", "True", "Observed"),
                     values = c("#00000033",
                                "blue",
                                "orange")) +
  theme_bw() +
  theme(legend.position = "bottom")
