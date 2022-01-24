#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#------------- Module 01 -- Introduction --------------X
#----------------Last updated 2021-01-06---------------X
#-------------------Exercise Solution------------------X
#######################################################X

# Load packages ----
library(tidyverse)
library(amt)
library(lubridate)

# Question 1 ----
# Read the data set `data/elephants.csv` (note file paths always start globally). Create a track using, make sure that you set the CRS to 4326. Ensure that the timezone is "GMT"

dat <- read_csv("data/elephants.csv")
dat
tz(dat$timestamp) 
tz(dat$timestamp) <- "GMT"

# In case you use read.csv, you will have to parse the date first
dat <- read.csv("data/elephants.csv")
dat
head(dat)
str(dat)
dat$timestamp <- ymd_hms(dat$timestamp, tz = "GMT")

# Now lets create a track, set the CRS and keep all columns
trk1 <- make_track(dat, long, lat, timestamp, id = id, temperature = temperature, crs = 4326)
trk1

# or a bit shorter
trk1 <- make_track(dat, long, lat, timestamp, all_cols = TRUE, crs = 4326)
trk1

# Question 2 ----
# Transform the track to a projected UTM CRS. You can use 32630 (https://epsg.io/32630) and filter only for the individual `"Salif Keita"`. Only use data for the year 2009. What is is the sampling rate?

trk2 <- trk1 %>% transform_coords(crs_to = 32630) %>% 
  filter(id == "Salif Keita", year(t_) == 2009)

from_to(trk1)
from_to(trk2)

summarize_sampling_rate(trk2)

# Question 3 ----  
# Calculate the speed (i.e., the distance per unit time) for each relocation, and the time of the day for each relocation. 

trk2$speed <- speed(trk2)
trk2 %>% mutate(speed = speed(trk2))
trk2 %>% mutate(speed = speed(.))

time_of_day(trk2)

# All together
trk3 <- trk2 %>% mutate(speed = speed(.)) %>% 
  time_of_day()


 
# Question 4 ----   
# Does speed differ between time of day and temperature

theme_set(theme_light())

ggplot(trk3, aes(tod_, speed)) +
  geom_boxplot() +
  labs(x = "Time of day", y = "Speed [m/s]")

ggplot(trk3, aes(temperature, speed)) +
  geom_point(alpha = 0.05) +
  geom_smooth() +
  labs(x = "Temperature [°C]", y = "Speed [m/s]")


ggplot(trk3, aes(temperature, speed)) +
  geom_point(alpha = 0.05) +
  geom_smooth() +
  facet_wrap(~ tod_) +
  labs(x = "Temperature [°C]", y = "Speed [m/s]")


# Question 5 ---- 
# Add a new column to the data set that indicates the months of the year when each relocation was taken. Does the association between speed and temperature remains the same for each month? Explore this graphically and calculate monthly means and median for speeds.


trk4 <- mutate(trk3, month = month(t_, label = TRUE, locale = "en_US.UTF-8")) # note you can remove the lcoale
mutate(trk3, month = month(t_, label = TRUE)) %>% pull(month)

ggplot(trk4, aes(temperature, speed)) +
  geom_point(alpha = 0.05) +
  geom_smooth() +
  facet_grid(~ month) +
  labs(x = "Temperature [°C]", y = "Speed [m/s]")

ggplot(trk4, aes(temperature, speed, col = tod_)) +
  geom_point(alpha = 0.05) +
  geom_smooth() +
  facet_grid(~ month) +
  labs(x = "Temperature [°C]", y = "Speed [m/s]") +
  theme(legend.position = "bottom")

# Explore this with descriptive stats
trk4 %>% group_by(month, tod_) %>% 
  summarise(m = mean(speed), 
            sd = sd(speed))

# We could also use a `nest()` approach. This is a bit more code here. 
trk4 %>% nest(data = -c(month, tod_)) %>% 
  mutate(m = map_dbl(data, ~ mean(.x$speed)), 
         sd = map_dbl(data, ~ sd(.x$speed))) %>% 
  select(-data)

# We can display this graphically
trk4 %>% group_by(month, tod_) %>% 
  summarise(mean = mean(speed, na.rm = TRUE), 
            median = median(speed, na.rm = TRUE)) %>% 
  pivot_longer(cols = c(mean, median)) %>% 
  ggplot(aes(month, value, col = name, group = name)) +
  geom_point() +
  geom_line(alpha = 0.2) +
  facet_wrap(~ tod_) +
  labs(x = "Temperature [°C]", y = "Speed [m/s]", col = "Statistic") +
  theme(legend.position = "bottom")
