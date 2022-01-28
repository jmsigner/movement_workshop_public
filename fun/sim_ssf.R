#' Function to simulate a movement track from a SSF
#' 
#' @param n_steps `[numeric]` The number of movement steps.
#' @param n_ch `[numeric]` The number of choices at each time step. 
#' @param l `[numeric>0]` Rate parameter for exponential step length distribution. 
#' @param coef `[numeric]` Selection coefficients.
#' @param xy0 `[numeric]` The starting position.
#' @param resc `[RasterLayer]` The environmental covariates.

simulate_ssf <- function(n_steps, n_ch, l = 1, coef, xy0, resc) {
  
  ###
  # n_steps = 500
  # n_ch = 10
  # l = 0.2
  # xy0 = c(100, 100)
  # resc = covars
  # coef = c(0.01, 0)
  ###
  
  sl <- rexp(n_steps * n_ch, rate = l)
  ta <- runif(n_steps * n_ch, -pi, pi)
  
  steps <- rep(1:n_steps, each = n_ch)
  x_0 <- xy0[1]
  y_0 <- xy0[2]
  x_s <- sl * sin(ta)
  y_s <- sl * cos(ta)
  x <- rep(NA, n_steps)
  y <- rep(NA, n_steps)
  x[1] <- x_0
  y[1] <- y_0
  # multiply resources with selection coef
  for (i in 1:length(coef)) resc[[i]] <- resc[[i]] * coef[i]
  for (i in 2:n_steps) {
    x_pos <- x[i - 1] + sl[steps == i] * sin(ta[steps == i])
    y_pos <- y[i - 1] + sl[steps == i] * cos(ta[steps == i])
    p <- exp(rowSums(raster::extract(resc, cbind(x_pos, y_pos)))) 
    w <- sample(n_ch, 1, prob = p)
    x[i] <- x_pos[w]
    y[i] <- y_pos[w]
  }
  make_track(data.frame(x = x, y = y, t = lubridate::ymd_hms("2020-01-01 00:00:00") + lubridate::hours(1:n_steps)), x, y, t)
}
