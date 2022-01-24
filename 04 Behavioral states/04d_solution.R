library(moveHMM)
library(tidyverse)
library(lubridate)
library(amt)

# 1. Read the data set `data/elephants.csv` (note file paths always start
# globally). Create a track using, make sure that you set the CRS to 4326.
# Ensure that the timezone is "GMT".

dat1 <- read_csv("data/elephants.csv")
dat1 %>% count(id)

# 2. Transform the track to a projected UTM CRS. You can use 32630
# (https://epsg.io/32630) and filter only for the individual `"Salif Keita"`.
# Only use data for the year 2009.

dat2 <- dat1 %>% filter(id == "Salif Keita", year(timestamp) == 2009) %>% 
  make_track(long, lat, timestamp, temperature, crs = 4326) %>% 
  transform_coords(32630) %>% time_of_day() %>% 
  mutate(day = tod_ == "day")

# 3. Create a suitable data set for the `moveHMM` package and prepare data to
# fit a HMM.

dat3 <- data.frame(ID = 1, x = dat2$x_, y = dat2$y_, tmp = dat2$temperature, 
                   t = dat2$t_, day = dat2$day)
dat.hmm <- prepData(dat3, type = "UTM")

hist(dat.hmm$step)
step.start.mean <- c(100, 2000)
step.start.sd <- c(100, 2000)

hist(dat.hmm$angle)
angle.start.mean <- c(0, 0)
angle.start.sd <- c(1, 1)

# 4. Fit two model (`m1` and `m2`). For each model use two states. For the
# second model use temperature as covariate that can effect the transition
# probabilities between states. Which of the two models would you choose?

m1 <- fitHMM(dat.hmm, nbStates = 2, 
       stepPar0 = c(step.start.mean, step.start.sd), 
       anglePar0 = c(angle.start.mean, angle.start.sd))

plot(m1)

# Fit a model with temperature as covariate
head(dat.hmm)

m2 <- fitHMM(dat.hmm, nbStates = 2, formula = ~ tmp,
       stepPar0 = c(step.start.mean, step.start.sd), 
       anglePar0 = c(angle.start.mean, angle.start.sd))

plot(m2)

AIC(m1)
AIC(m2)
