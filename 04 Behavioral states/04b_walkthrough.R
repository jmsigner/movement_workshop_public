#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#-----------Module 04 -- Behavioral States ------------X
#----------------Last updated 2021-12-27---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

# Load packages ----
library(tidyverse)
library(amt)
library(moveHMM)
library(lubridate) # to deal with date and time

# Simulate data ----
# ... One state ----

# We start with simulating movement data for one behavioral state (i.e., the
# animal is not switching between different states).

set.seed(130)

# Simulation set up
n <- 1000 # the number f time steps, that we want to simulate
x <- y <- rep(NA, n)
x[1] <- y[1] <- 0 # We start our simulation at point 0, 0 in space

# Gamma mean and sd The package moveHMM characterizes the step-length
# distribution in terms of mean and sd.
sl.par <- c(4, 2)  # mean, sd.

# The scale of a gamma distribution is given by dividing the sd by the mean.
scales <- sl.par[2] / sl.par[1]
# Shape is given by dividing the mean by the scale
shapes <- sl.par[1] / scales

# The expected value (i.e, the mean is given again by multiplying the scale with
# the shape)
scales * shapes

# We can now simulate the path of the animal by drawing random step lengths from
# a gamma distribution with the shape and scale we defined previously and turn
# angles are uniformly distributed (i.e., no directional persistance).


for (i in 1:(n - 1)) {
# At each time step we draw one value for sl and ta and then add this to the
# previouse position.
  sl <- rgamma(1, shape = shapes[1], scale = scales[1])
  ta <- runif(1, min = -pi, max = pi)
  
  x[i + 1] <- x[i] + sl * cos(ta)
  y[i + 1] <- y[i] + sl * sin(ta)
  
}

plot(x, y, type = "b", asp = 1, pch = 20, col = adjustcolor("black", 0.1))

# Next lets calculate the distances between consecutive relocations
d <- sqrt(diff(x) ^ 2 + diff(y) ^ 2)
hist(d)

# And we can now fit again a gamma distribution to the observed step lengths. 
fitdistrplus::fitdist(d, "gamma")
shapes
1 / scales

# ... Two states ----
# ... ... Data generation -----

# We simulate again, but this time using two states. This means the animal
# switches between different states. Each state has its own set of parameters
# for the step-length distribution. The turn-angle distribution does not differ.

# Simulation set up
n <- 1000
x <- y <- rep(NA, n)
x[1] <- y[1] <- 0

# Gamma mean and sd
sl.par <- c(1, 10, 1, 5) # mean1, mean2, sd1, sd2

scales <- sl.par[3:4] / sl.par[1:2]
shapes <- sl.par[1:2] / scales
scales * shapes

Ts <- n - 1 # The number of time steps
Ks <- 2 # The number of states that we want to simulate
G <- cbind(c(0.9, 0.2), c(0.1, 0.8)) # This is the transition matrix

# Next we need to simulate the true states (that we cant not observe)
Z = rep(0, Ts) # Z holds the true state
Z[1] = 1 # We start with state 1
for (t in 1:(Ts - 1)) {
  # At each time point we sample either if we stay in in the current state or 
  # move to a different state according to the transition matrix G. 
  Z[t + 1] = sample(Ks, size = 1, prob = G[Z[t],]) 
}

# These are the true (but later unknown) states. 
Z

# Now we simulate again the relocations for each animal. But this time the
# parameters of the gamma distribution differ.
for (i in 1:(n - 1)) {
  sl <- rgamma(1, shape = shapes[Z[i]], scale = scales[Z[i]])
  ta <- runif(1, min = -pi, max = pi)
  
  x[i + 1] <- x[i] + sl * cos(ta)
  y[i + 1] <- y[i] + sl * sin(ta)
  
}

plot(x, y, type = "b", asp = 1, pch = 20, col = adjustcolor("black", 0.1))

# ... ... Model fitting -----
# ... ... ... naive -----
d <- sqrt(diff(x) ^ 2 + diff(y) ^ 2)
hist(d)

# Ignoring the different states
(res1 <- fitdistrplus::fitdist(d[Z == 1], "gamma"))
shapes[1]
(1 / scales)[1] # to get the rate


# Lets see how this look like in a plot
hist(d, prob = TRUE, ylim = c(0, 0.9))
curve(dgamma(x, shape = res1$estimate["shape"], rate = res1$estimate["rate"]), 
      add = TRUE, col = "red")

# We can now fit two Gamma distributions (one for each state). We know the true
# state Z, because we are in a simulation setting.
(res1 <- fitdistrplus::fitdist(d[Z == 1], "gamma"))
shapes[1]
(1 / scales)[1] # to get the rate


(res2 <- fitdistrplus::fitdist(d[Z == 2], "gamma"))
shapes[2]
(1 / scales)[2] # to get the rate

# Again lets create a plot
hist(d, prob = TRUE, ylim = c(0, 0.9))
curve(dgamma(x, shape = res1$estimate["shape"], rate = res1$estimate["rate"]), 
      add = TRUE, col = "red")
curve(dgamma(x, shape = res2$estimate["shape"], rate = res2$estimate["rate"]), 
      add = TRUE, col = "blue")

# But of course we do not know `Z`! What we can do instead is to figure out a
# probability for each step to be in either state. That is what a HMM does.

# ... ... ... HMM -----
# The package `moveHMM` provides an implementation of HMMs specially tailored
# for movement data.

# We start of with a data.frame, with regularly sampled data. If data are not
# equally spaced check out the multiple imputation approach introduced in the
# previous module. The package `momentuHMM` also provides multiple imputation
# for HMMs.
dat <- data.frame(x = x, y = y)

# The function `prepData()` provides an easy way to format data correctly for
# the moveHMM package. The function is similar to the `steps()` function in
# `amt`.
dat <- prepData(dat, type = "UTM")

# We then have to provide starting values. Parameter estimates are sensitive to
# the starting values, so it might be good test several possible starting
# values. Starting values are specified in terms of the mean and sd for the step
# length and turning angles.

# Note, since we did not include the turn angle in our simulations, we will not
# estimate them here.
sl.start <- c(1, 10, 1, 5) # mean1, mean2, sd1, sd2

# The function `fitHMM()` fits a HMM to the data. The argument `nbStates`
# declares the number of states that will be fitted.
m1 <- fitHMM(dat, nbStates = 2, stepPar0 = sl.start, 
             stepDist = "gamma", angleDist = "none") 

m1

# The function `stateProbs()` returns a probability for each observation to be
# in each state.

stateProbs(m1)

# The viterbi algorithm is used to decode the states. This means each relocation
# is assigned a number refereing to the state. The function `viterbi()` can be
# used to do this.

Zest <- viterbi(m1)

# Since we know the truth (`Z`), we can calculate the classification accurracy. 
mean(Z == Zest[-1000]) ## Accuracy > 99 %

## Plots
plot(m1)

# Finally we can check the pseudo residuals of the HMM. 
resid1 <- pseudoRes(m1)
qqnorm(resid1$stepRes)
qqline(resid1$stepRes, col = "red")

# As mentioned before, its always good to check different starting values. Thus
# we could use a `lapply()` to iterate over different starting values and check
# of the loglikelihood remains the same.

start <- lapply(1:10, function(x) c(runif(2, min = c(0.05, 10), max = c(10, 20)),
             runif(2, min = c(0.05, 10), max = c(10, 20))))
m1.many <- lapply(start, function(x) {
  fitHMM(dat, 2, stepPar0 = x,  stepDist = "gamma", angleDist = "none")
})

# Comparing the negative loglikelihood among all 10 models and the model we
# fitted before, we see that we are all good.
min(sapply(m1.many, function(x) x$mod$minimum))
m1$mod$minimum

# ... Three states ----

# Below we repeat the previous section, but this time we simulate data from a
# three state HMM. We will explore the effects of fitting a two and three state
# HMM for this data set.

# ... ... Data generation ----

n <- 1000
x <- y <- rep(NA, n)
x[1] <- y[1] <- 0

# Gamma mean and sd
sl.par <- c(1, 10, 50, 1, 5, 10) # mean1, mean2, mean3, sd1, sd2, sd3

scales <- sl.par[4:6] / sl.par[1:3]
shapes <- sl.par[1:3] / scales
scales * shapes

Ts <- n - 1
Ks <- 3 # the number of states
G <- cbind(c(0.8, 0.2, 0.5), c(0.1, 0.5, 0.3), c(0.1, 0.3, 0.2)) 

Z = rep(0, Ts)
Z[1] = 1
for (t in 1:(Ts - 1)) {
  Z[t + 1] = sample(Ks, size = 1, prob = G[Z[t],])
}

# Simulating the relocation
for (i in 1:(n - 1)) {
  sl <- rgamma(1, shape = shapes[Z[i]], scale = scales[Z[i]])
  ta <- runif(1, min = -pi, max = pi)
  
  x[i + 1] <- x[i] + sl * cos(ta)
  y[i + 1] <- y[i] + sl * sin(ta)
  
}

plot(x, y, type = "b", asp = 1, pch = 20, col = adjustcolor("black", 0.1))

# ... ... Model fitting

dat <- data.frame(x = x, y = y)
dat <- prepData(dat, type = "UTM")


# ... ... ... Two state HMM
sl.start <- c(1, 10, 1, 5) # mean1, mean2, sd1, sd2
m1 <-fitHMM(dat, 2, stepPar0 = sl.start, stepDist = "gamma", angleDist = "none")

# Decode states
Zest <- viterbi(m1)
mean(Z == Zest[-1000]) # Accuracy dropped

## Plots
plot(m1)

# Check residuals
resid1 <- pseudoRes(m1)
qqnorm(resid1$stepRes)

# ... ... ... Three state HMM
sl.start <- c(1, 10, 20, 1, 5, 20) # mean1, mean2, mean3, sd1, sd2, sd3
m2 <- fitHMM(dat, 3, stepPar0 = sl.start,
         stepDist = "gamma",
         angleDist = "none")


# Decode states
Zest <- viterbi(m2)
mean(Z == Zest[-1000]) ## Accuracy > 99 %

## Plots
plot(m2)

# Check residuals
resid1 <- pseudoRes(m2)
qqnorm(resid1$stepRes)

# ... Two states with env covariates ----

# In the last scenario we simulate a model where a covariate `x` effects the
# transition probability between states.

# ... ... Data generation

# Gamma mean and sd
sl.par <- c(1, 10, 1, 5) # mean1, mean2, sd1, sd2

scales <- sl.par[3:4] / sl.par[1:2]
shapes <- sl.par[1:2] / scales

# Simulating a covariate x
x <- cumsum(rnorm(n))
plot(x, type = "l")

# coefficients
beta.mat <- matrix(c(-3, 0.1, 0.5, -.2), nrow = 2, byrow = TRUE)
beta.mat

Z <- rep(NA, n)
Z[1] <- 1
for (i in 2:n) {
  gamma <- diag(2)
  gamma[!gamma] <- exp(beta.mat[, 1] + beta.mat[, 2] * x[i])
  gamma <- gamma / apply(gamma, 1, sum)
  Z[i] <- sample(x = c(1, 2), size = 1, prob = gamma[Z[i - 1],])
}

Z

# Now lets simulate the movement path 
for (i in 1:(n - 1)) {
  sl <- rgamma(1, shape = shapes[Z[i]], scale = scales[Z[i]])
  ta <- runif(1, min = -pi, max = pi)
  
  x[i + 1] <- x[i] + sl * cos(ta)
  y[i + 1] <- y[i] + sl * sin(ta)
  
}

plot( x, y, type = "b", asp = 1, pch = 20, col = adjustcolor("black", 0.1))

d <- sqrt(diff(x) ^ 2 + diff(y) ^ 2)
hist(d)

# ... ... Fit HMM
dat <- data.frame(x = x, y = y, covar = x)
dat <- prepData(dat, type = "UTM")

sl.start <- c(1, 10, 1, 5) # mean1, mean2, sd1, sd2 

m1 <- fitHMM(data = dat, nbStates = 2, 
         stepPar0 = sl.start,
         stepDist = "gamma",
         angleDist = "none")

m2 <- fitHMM(data = dat, nbStates = 2, formula = ~ covar, 
         stepPar0 = sl.start,
         stepDist = "gamma",
         angleDist = "none")

# Visualize
plot(m2)


# We can check residuals
r1 <- pseudoRes(m1)$stepRes
r2 <- pseudoRes(m2)$stepRes

qqnorm(r1)
qqline(r1, col = "red")

qqnorm(r2)
qqline(r2, col = "red")

# And see if AIC supports the more complicated model.
AIC(m1)
AIC(m2)

