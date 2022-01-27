#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#----------------Module 08 -- iSSF pt 2----------------X
#----------------Last updated 2021-01-26---------------X
#-------------------Exercise Solution------------------X
#######################################################X

# Using UT cougar data

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

# Format as track_xyt
trk <- dat %>% 
  make_track(x_, y_, t_, crs = 32612)

# Do we have regular steps?
summarize_sampling_rate(trk)

# The median sampling rate is 4h, but the track isn't perfectly regular.
# Let's separate into bursts of 4h steps, with a tolerance of 15 minutes.
# Only keep bursts that have at least 3 relocations. Then we can turn those
# locations into steps.

stp <- track_resample(trk, rate = hours(8), tolerance = minutes(15)) %>% 
  filter_min_n_burst(min_n = 3) %>% 
  steps_by_burst()

# What do our step durations look like now?
hist(as.numeric(stp$dt_))

# 2. Consider habitat and movement variables ----

## Habitats
plot(hab)

# We have:
#   - elevation (m)
#   - tree cover (%)
#   - biomass of annual and perennial grasses & forbs (kg/ha)
#   - distance to road (m)

# I would categorize these as:
#   - elevation: condition during summer, risk during winter (too much
#                 snow at high elevation?)
#   - tree cover: resource during hot summer days, provides shade
#   - biomass: resource -- maybe high biomass of forage attracts prey species
#   - distance to road: risk -- cougars are hunted as a game species in UT
#                       and probably generally perceive humans as a risk.

# Notice the temporal dynamics in some of my a priori predictions. Some of
# my predictions depend on time of day (e.g., trees as shade resource) and some
# depend on season (e.g., only need shade during summer).

# We can account for these dynamics with interactions. For example, a 
# 3-way interaction between trees, time of day, and season could capture
# the dynamics I've hypothesized, but that requires a lot of parameters.
# Additionally, how should I categorize season? Summer by itself? Or are
# spring, and autumn different from each other? (we don't have much winter data
# for this cougar).

# My biological question should guide the model structure, but a valid
# alternative to 3-way interactions would be to split the data into a 
# dataset of summer months and a dataset of spring/fall months. Then I 
# can only use a 2-way interaction between season and time of day to 
# capture all my hypotheses.

# Again, let your understanding of the system and your desired inference guide
# your decision. Here, I'll split into the two seasons and fit separate models.

## Movement variables

# We are going to use a gamma distribution to model step lengths, so we need
# to include 'sl_' and 'log(sl_)' in our model.

# We are going to use a von Mises distribution to model turn angles, so we 
# need to include 'cos(ta_)' in our model.

# How might these interact? Maybe high elevations also have a lot of snow
# during winter. So I expect that at high elevation, step lengths will be
# shorter during winter.

# There may also be daily patterns. Cougars primarily hunt at night, so perhaps
# steps are longer at night than during the day.

# 3. Fit iSSFs ----
issf_dat <- stp %>% 
  # Keep the defaults for generating available steps (gamma and von Mises distrs)
  # Increase the number of random steps
  random_steps(n_control = 20) %>% 
  # Attach habitat variables
  extract_covariates(hab) %>% 
  # Add temporal covariates
  time_of_day(where = "both") %>% 
  # Add additional movement covariates
  mutate(log_sl_ = log(sl_),
         cos_ta_ = cos(ta_))

# Split into summer and other
summ_dat <- issf_dat %>% 
  filter(month(t1_) %in% 6:9)
other_dat <- issf_dat %>% 
  filter(!month(t1_) %in% 6:9)

# With two different models, we could also use different model structures
# (i.e., formulas) for the different seasons. For example, we can have an
# elevation x movement interaction only for the non-summer data.

# Summer model
summ_issf <- summ_dat %>% 
  fit_issf(case_ ~ 
             # Habitat selection main effects
             elevation + I(elevation^2) + trees + biomass + dist_to_road +
             # Habitat interactions
             trees : tod_end_ + 
             # Movement main effects
             sl_ + log_sl_ + cos_ta_ +
             # Movement interactions
             sl_ : tod_start_ + log_sl_ : tod_start_ + cos_ta_ : tod_start_ +
             # Don't forget the strata
             strata(step_id_),
           # And include model = TRUE so we can use 'log_rss()' later
           model = TRUE)

summary(summ_issf) 

# Other seasons
other_issf <- other_dat %>% 
  fit_issf(case_ ~ 
             # Habitat selection main effects
             elevation + I(elevation^2) + trees + biomass + dist_to_road +
             # Habitat interactions
             trees : tod_end_ + 
             # Movement main effects
             sl_ + log_sl_ + cos_ta_ +
             # Movement interactions
             sl_ : tod_start_ + log_sl_ : tod_start_ + cos_ta_ : tod_start_ +
             sl_ : elevation + log_sl_ : elevation +
             # Don't forget the strata
             strata(step_id_),
           # And include model = TRUE so we can use 'log_rss()' later
           model = TRUE)

summary(other_issf) 

# 4. RSS Figures ----
# Let's examine our hypothesis that selection for trees depends on time of day
# during the summer but not during the other seasons.

# Summer - day
x1_summer_day <- data.frame(elevation = 2500,
                            trees = seq(0, 100, length.out = 25),
                            biomass = 2000,
                            dist_to_road = 1500,
                            tod_start_ = factor("day", 
                                                levels = c("day", "night")),
                            tod_end_ = factor("day", 
                                              levels = c("day", "night")), 
                            sl_ = 900,
                            log_sl_ = log(900),
                            cos_ta_ = 1)

x2_summer_day <- data.frame(elevation = 2500,
                            trees = 0,
                            biomass = 2000,
                            dist_to_road = 1500,
                            tod_start_ = factor("day", 
                                                levels = c("day", "night")),
                            tod_end_ = factor("day", 
                                              levels = c("day", "night")), 
                            sl_ = 900,
                            log_sl_ = log(900),
                            cos_ta_ = 1)

rss_summer_day <- log_rss(summ_issf, x1_summer_day, x2_summer_day, ci = "se")

# Summer - night
x1_summer_night <- data.frame(elevation = 2500,
                              trees = seq(0, 100, length.out = 25),
                              biomass = 2000,
                              dist_to_road = 1500,
                              tod_start_ = factor("night", 
                                                  levels = c("day", "night")),
                              tod_end_ = factor("night", 
                                                levels = c("day", "night")), 
                              sl_ = 900,
                              log_sl_ = log(900),
                              cos_ta_ = 1)

x2_summer_night <- data.frame(elevation = 2500,
                              trees = 0,
                              biomass = 2000,
                              dist_to_road = 1500,
                              tod_start_ = factor("night", 
                                                  levels = c("day", "night")),
                              tod_end_ = factor("night", 
                                                levels = c("day", "night")), 
                              sl_ = 900,
                              log_sl_ = log(900),
                              cos_ta_ = 1)

rss_summer_night <- log_rss(summ_issf, x1_summer_night, x2_summer_night, ci = "se")

# Other - day
x1_other_day <- data.frame(elevation = 2500,
                           trees = seq(0, 100, length.out = 25),
                           biomass = 2000,
                           dist_to_road = 1500,
                           tod_start_ = factor("day", 
                                               levels = c("day", "night")),
                           tod_end_ = factor("day", 
                                             levels = c("day", "night")), 
                           sl_ = 900,
                           log_sl_ = log(900),
                           cos_ta_ = 1)

x2_other_day <- data.frame(elevation = 2500,
                           trees = 0,
                           biomass = 2000,
                           dist_to_road = 1500,
                           tod_start_ = factor("day", 
                                               levels = c("day", "night")),
                           tod_end_ = factor("day", 
                                             levels = c("day", "night")), 
                           sl_ = 900,
                           log_sl_ = log(900),
                           cos_ta_ = 1)

rss_other_day <- log_rss(other_issf, x1_other_day, x2_other_day, ci = "se")

# Other - night
x1_other_night <- data.frame(elevation = 2500,
                             trees = seq(0, 100, length.out = 25),
                             biomass = 2000,
                             dist_to_road = 1500,
                             tod_start_ = factor("night", 
                                                 levels = c("day", "night")),
                             tod_end_ = factor("night", 
                                               levels = c("day", "night")), 
                             sl_ = 900,
                             log_sl_ = log(900),
                             cos_ta_ = 1)

x2_other_night <- data.frame(elevation = 2500,
                             trees = 0,
                             biomass = 2000,
                             dist_to_road = 1500,
                             tod_start_ = factor("night", 
                                                 levels = c("day", "night")),
                             tod_end_ = factor("night", 
                                               levels = c("day", "night")), 
                             sl_ = 900,
                             log_sl_ = log(900),
                             cos_ta_ = 1)

rss_other_night <- log_rss(other_issf, x1_other_night, x2_other_night, ci = "se")

# Now grab all the data.frames from the 'log_rss' objects and combine
fig_dat <- bind_rows("Summer Day" = rss_summer_day$df,
                     "Summer Night" = rss_summer_night$df,
                     "Other Day" = rss_other_day$df,
                     "Other Night" = rss_other_night$df,
                     .id = "season_time") %>% 
  # Split "season_time" into two columns
  mutate(season = word(season_time, 1, 1),
         time = word(season_time, 2, 2)) %>% 
  # Convert log-RSS to RSS
  mutate(rss = exp(log_rss),
         rss_lwr = exp(lwr),
         rss_upr = exp(upr))

# Check
head(fig_dat)

# Plot
fig_dat %>% 
  ggplot(aes(x = trees_x1, y = log_rss, ymin = lwr, ymax = upr,
             color = time, fill = time)) +
  facet_wrap(~ season) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_ribbon(linetype = "dashed", alpha = 0.2, size = 1) +
  geom_line(size = 1) +
  xlab("Tree Cover (%)") +
  ylab("log-RSS") +
  theme_bw()

# You can see that the 95% confidence intervals all overlap with log-RSS of 0,
# except for during the day during summer. This seems like good support for
# our hypothesis that our cougar uses trees for shade during summer.

# 5. Movement figures ----
# Now let's make a figure to look at our hypothesis that movement slows down
# at high elevation during other seasons. (Remember, we left this interaction
# out of the summer model).

# Note that none of the interaction terms for the step-length distribution
# are very large, so the answer is most likely no. But let's go through the
# exercise anyway.
summary(other_issf)

# Get betas from our other-season model
b <- coef(other_issf)

# Update step-length distributions for different elevations.
elev <- seq(1500, 3500, length.out = 5)


gamma_elev <- expand.grid(elev = elev,
                          tod_start_ = c("day", "night")) %>% 
  # Recall, that the betas for sl_ and for log_sl_ are now functions of elevation
  # and time of day (day is the reference).
  mutate(b_sl_ = b[["sl_"]] + b[["elevation:sl_"]] * elev + 
           b[["sl_:tod_start_night"]] * (tod_start_ == "night"), 
         b_log_sl_ = b[["log_sl_"]] + b[["elevation:log_sl_"]] * elev + 
           b[["log_sl_:tod_start_night"]] * (tod_start_ == "night"))

# Now that we have the betas for sl and log_sl as a function of elev and tod,
# we can update our gamma distribution.
# (notice that this is vectorized)
gamma_elev_distr <- update_gamma(sl_distr(other_issf), 
                                 gamma_elev$b_sl_,
                                 gamma_elev$b_log_sl_)

# And now add the shape and scale to our data.frame
gamma_elev$shp <- gamma_elev_distr$params$shape
gamma_elev$scl <- gamma_elev_distr$params$scale

# Check
gamma_elev

# Plot mean step length
gamma_elev %>% 
  mutate(mean_sl = shp * scl) %>% 
  ggplot(aes(x = elev, y = mean_sl, color = tod_start_)) +
  geom_line(size = 1) +
  xlab("Elevation (m)") +
  ylab("Mean Step Length (m)") +
  scale_color_discrete(name = "Time") +
  theme_bw()

# Mean step length is longer at night than during the day, but both
# seem to increase with elevation.

# Now we want to calculate the probability density under the gamma with those
# parameters for a range of step lengths

step_fig_dat <- gamma_elev %>% 
  # Convert to tibble for nested data.frame
  as_tibble() %>%
  # Add elevation as list column
  mutate(elev_list = list(tibble(sl = seq(1, 1000, length.out = 100)))) %>% 
  # Unnest
  unnest(cols = elev_list) %>% 
  # Calculate probability density
  mutate(dens = dgamma(sl, shape = shp, scale = scl))


# Plot distribution
ggplot(step_fig_dat, aes(x = sl, y = dens, col = elev, group = elev)) +
  facet_wrap(~ tod_start_) +
  geom_line() +
  xlab("Step Length (m)") +
  ylab("Proability Density") +
  theme_bw()

# We can tell there is barely a difference in the distributions at different
# elevations.

# Maybe that means there is no effect of elevation, or maybe that effect
# isn't very strong because we don't have much data in the core winter
# months.

