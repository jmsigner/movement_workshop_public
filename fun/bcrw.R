# Function to simulate a biased correlated random walk

#' @param start_loc `[numeric]` Vector of length 2 giving x and y coordinates
#' of starting location.
#' @param centroid `[numeric]` Vector of length 2 giving x and y coordinates of
#' home range centroid.
#' @param n_steps `[integer]` Integer giving number of steps to simulate.
#' @param sl_distr `[numeric]` *Named* vector of length 2 giving shape and
#' scale of gamma distribution from which step lengths are drawn.
#' @param rho `[numeric]` Value between 0 and 1 giving correlation of turning
#' angles.
#' @param beta `[numeric]` Value between 0 and 1 giving the strength of the
#' attraction towards the centroid, *i.e.*, strength of the bias.
#' 
bcrw <- function(start_loc,
                 centroid,
                 n_steps,
                 sl_distr = c('shape' = 0.75, 'scale' = 200),
                 rho = 0.25,
                 beta = 0.1) {
  
  # Check for CircStats
  if (!require(CircStats)) {
    stop("This function requires the package 'CircStats' to be installed prior to use.")
  }
  
  # Store start location
  X <- Y <- numeric()
  X[1] <- start_loc[1]
  Y[1] <- start_loc[2]
  
  # Initialize vector of bearings to centroid
  bearing <- numeric()
  
  # Calculate first bearing to centroid
  dx <- centroid[1] - X[1]
  dy <- centroid[2] - Y[1]
  bearing[1] <- pi/2 - atan2(y = dy, x = dx)
  
  # Get x and y components of bias
  xa <- ya <- numeric() # unit delta x and delta y
  ya[1] <- (1 - beta) * sin(0) + beta * sin(bearing[1])
  xa[1] <- (1 - beta) * cos(0) + beta * cos(bearing[1])
  
  # Expected bearing based on bias
  bias_angle <- numeric()
  bias_angle[1] <- 2 * atan(ya[1]/(sqrt(ya[1]^2 + xa[1]^2) + xa[1]))
  
  # Actual turning angle (with correlation)
  ta <- numeric()
  ta[1] <- bias_angle[1]
  
  # Now loop over steps\
  for (i in 2:n_steps){
    
    # Calculate bearing to centroid
    dx <- centroid[1] - X[i - 1]
    dy <- centroid[2] - Y[i - 1]
    bearing[i] <- pi/2 - atan2(y = dy, x = dx)
    
    # Expected angle based on bias
    ya[i] <- (1 - beta) * sin(bias_angle[i - 1]) + beta * sin(bearing[i])
    xa[i] <- (1 - beta) * cos(bias_angle[i - 1]) + beta * cos(bearing[i])
    bias_angle[i] <- 2 * atan(ya[i]/(sqrt(ya[i]^2 + xa[i]^2) + xa[i]))
    
    # Draw turn angle from wrapped Cauchy
    ta[i] <- CircStats::rwrpcauchy(n = 1,
                                   location = bias_angle[i],
                                   rho = rho)
    
    # Draw step length from Gamma
    sl <- stats::rgamma(n = 1,
                        shape = sl_distr["shape"],
                        scale = sl_distr["scale"])
    
    # Calculate next coordinates
    DX <- sl * sin(ta[i])
    DY <- sl * cos(ta[i])
    X[i] <- X[i - 1] + DX
    Y[i] <- Y[i - 1] + DY
  }
  
  # Compile data.frame
  df <- data.frame(t = 1:n_steps,
                   x = X,
                   y = Y)
  
  # Return
  return(df)
}