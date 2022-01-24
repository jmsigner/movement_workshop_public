# Fitting simple example from lecture with logistic regression

# Used points
u <- data.frame(h = c(
  rep("G", 800),
  rep("F", 150),
  rep("aW", 50)),
  y = 1)

# Available points
a <- data.frame(h = c(
  rep("G", 9000),
  rep("F", 500),
  rep("aW", 500)),
  y = 0)

# Data
d <- rbind(u, a)
# Weights
d$w <- ifelse(d$y == 1, 1, 1e5)

# Fit GLM
m <- glm(y ~ h, weights = w, data = d, family = binomial)

# Check summary
summary(m)

# Exponentiate coefficents (= RSS for habitat vs. wetland)
exp(coef(m))
