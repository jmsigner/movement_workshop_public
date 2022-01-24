#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#--------------Module 01 -- Introduction --------------X
#----------------Last updated 2021-12-27---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

# Load packages ----
library(tidyverse)
library(amt)
library(lubridate) # to deal with date and time

# Load data ---- We will use data from fishers that were tracked by Scott
# LaPoint and has been widely used for many methodological comparisons, examples
# et. The data is saved in a `csv`-file in the data directory. Telemetry data
# are stored in many different formats, but often as some kind of text files.

# We can use R's functions to read these data into a `data.frame` or use the
# equivalent from the `readr` package to create a `tibble`.

dat <- read.csv("data/fisher.csv")

head(dat)

# We have six variables (columns):
# 1. and 2. `x_` and `y_` the coordinates for each data point.
# 3. `t_` the time stamp (i.e., when each data point was collected)
# 4. - 6. The `sex`, `id` and `name` of the animal each relocations belongs to. 

# We now have to make sure that the date column is correctly parsed. For `dat`
# we can use the function `ymd_hms` from the `lubridate` package. Yould also use
# the `strptime()` function from base R, but lubirdate's parsing functions are
# often a bit easier to memorize.

str(dat)
dat$ts <- ymd_hms(dat$t_, tz = "UTC")
head(dat$ts)
str(dat)
hour(dat$t_)
hour(dat$ts)
week(dat$ts)

# Important: if your time zone is different from UTC (the default), make sure
# you correctly specify the time zone with the argument `tz`.

# Note, if you `read_csv` from the `readr` package (which loaded with the
# tidyverse), the date is automatically parsed.  (i.e., R already "knows" this
# is a date). This is one advantage of the functions from the `readr` package.
# Note also that `read_csv()` returns a tibble while `read.csv()` returns a data
# frame.

dat1 <- read_csv("data/fisher.csv")
dat1

# We should verify though, that the time zone is set correctly. We can do this
# with the function `tz()`.
tz(dat1$t_)

head(dat1$t_)

# We can continue to work with `dat` or `dat1`. 

# Create tracks ----

# The basic building block of to work with the `amt` package are so called
# tracks. A track consists of a series of relocations. The function
# `make_track()` is used to create a track from a `data.frame` or `tibble`. It
# expects at least the data set, and coordinates (x and y). Optionally time
# stamps, additional columns and a coordinate reference system (crs) can be
# passed to the function using the EPSG code.

make_track(dat, x_, y_, ts)


# Note, that a track is characterized by `x_`, `y_` and `t_`. We could add a crs
# using the crs argument using the the EPSG code.

make_track(dat, x_, y_, ts, crs = 5070)

# If we want to retain the name column. There is also the argument `all_cols`
# which will keep all columns.
tr <- make_track(dat, x_, y_, ts, name = name, crs = 5070)
tr
class(tr)

# Working with tracks ----

# ... `dplyr` verbs ----

# Tracks are by design compatible with `dplyr`s verbs. For example, if we want
# to work with only one animal, we can just use the `filter()` function.

leroy <- tr %>% filter(name == "Leroy")
filter(tr, name == "Leroy")
leroy

# Other dplyr functions such as mutate, select, arrange, group_by, count,
# summarize work as well.

# For example, lets add the month of a relocation, 
leroy %>% mutate(month = month(t_))

# only select coordinates and and the year
leroy %>% mutate(month = month(t_)) %>% select(x_, y_, month)

# count the number of observations per month
leroy %>% mutate(month = month(t_)) %>% select(x_, y_, month) %>% 
  group_by(month) %>% summarize(n = n())

# The same could be achieved with `count()`
leroy %>% mutate(month = month(t_)) %>% select(x_, y_, month) %>% 
  count(month)

# Finally we can arrange the data by the number of relocations
leroy %>% mutate(month = month(t_)) %>% select(x_, y_, month) %>% 
  count(month) %>% arrange(n)

# We will continue to work just with `leroy` for the moment.

# ... Changing the CRS ----

# We can change the the CRS with the function `transform_coords()`. For example
# to change to geographic coordinates, we could just use:

leroy %>% transform_coords(4326)

# ... Visually inspect a track -----

# Function `inspect()` allows us to visually check a track with an interactive
# map.

leroy %>% inspect()


# ... Sampling rate ----

# The rate at which data are sampled for tracks can be different and irregular.
# To get an overview of the sampling rate the function
# `summarize_sampling_rate()` exists.

summarize_sampling_rate(leroy)

# This suggests that the median sampling rate is 15 min. We can now resample the
# track a 15 min interval.

leroy2 <- track_resample(leroy, rate = minutes(15), tolerance = seconds(60))
leroy2

nrow(leroy)
nrow(leroy2)

# Note that we lost 26 observations and gained one column (`burst_`). We lost
# observations that are not within the predefined sampling rate. And each
# relocation that we retained belongs now to a burst. A burst is a sequence of
# relocations with regular sampling rates (i.e., no missing relocations).

# Sometimes we want to restrict burst to a minimum of `n` relocations. This can
# be achieved with the function `filter_min_n_burst()`

leroy2 <- filter_min_n_burst(leroy2, min_n = 3)

# ... Movement attributes ----

# We can now caluculate for example step lengths with the function
# `step_length()`

leroy2 %>% step_lengths()

# If want to add step lengths to the data set we can use mutate()

leroy3 <- leroy2 %>% mutate(sl = step_lengths(.)) 
leroy3
# note the use of the `.` here, indicates, that we want to refer to the dataset
# that is currently under evaluation.

# We are ignoring bursts
leroy3 %>% group_by(burst_) %>% 
  summarize(fs = head(sl, 1), ls = tail(sl, 1)) %>% 
  pivot_longer(-burst_) %>% 
  ggplot(aes(name, value, group = burst_)) + geom_point(alpha = 0.1) + 
  geom_line(alpha = 0.1)


# Let honor bursts now for calculating step lengths (we need to nest our data here)
tmp <-  leroy2 %>% nest(data = -burst_) 
tmp$data[[4]] %>% mutate(sl = step_lengths(.))
tmp$data[[5]] %>% mutate(sl = step_lengths(.))
map(tmp$data, ~ .x %>% mutate(sl = step_lengths(.)))

leroy4 <- leroy2 %>% nest(data = -burst_) %>% 
  mutate(data = map(data, ~ mutate(.x, sl = step_lengths(.x)))) %>% 
  unnest(cols = data)

leroy4

leroy4 %>% group_by(burst_) %>% 
  summarize(fs = head(sl, 1), ls = tail(sl, 2)) %>% 
  pivot_longer(-burst_) %>% 
  ggplot(aes(name, value, group = burst_)) + geom_point(alpha = 0.1) + 
  geom_line(alpha = 0.1)


# ... Steps ---- 
# Next we want to change representations from individual
# locations to steps. A step consists of a start and end coordinate, a step
# length and a turn angle. The time difference between the start and the end
# point is constant.

leroy2 %>% steps() 

# We get a warning, because the function `steps()` be default ignores bursts.
# This is problematic if there is a large time gap between to consecutive
# points. To overcome this, we can use `steps_by_burst()`.

leroy2 %>% steps_by_burst() 

# The resulting tibble has 11 columns by default: 
# - `burst_`: the burst number.
# - `x1_` and `y1_`: the start coordinates of the step.
# - `x2_` and `y2_`: the end coordinates of the step.
# - `sl_`: the step length
# - `direction_p`: the direction of the step (relative to?)
# - `ta_`: the turn angle
# - `t1_` and `t2_`: the start and end time of a step.
# - `dt_`: te duration of a step.


summarize_sampling_rate()
