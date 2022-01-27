#' Get the angle between two points
#' @param xy Positions
get_angle <- function(xy) {
  ta <- base::atan2(xy[, 2], xy[, 1])
  ta <- ta - (pi/2)
  ta
}

#' Functions to transform parameters for the movement model (from the
#' distributions that describe the turn angles and step lengths) to coefficients
#' of the selection function.
#' @param param A parameter of a distribution.
#' @return A named vector of length one with the coefficient for the selection function.
#' @name Transformations
#' @export
kappa_to_cos_ta <- function(param) {
  checkmate::assert_number(param, lower = 0)
  c("cos_ta_" = param)
}

#' @rdname Transformations
#' @export
scale_to_sl <- function(param) {
  checkmate::assert_number(param, lower = .Machine$double.eps)
  c("sl_" = -1 / param )
}

#' @rdname Transformations
#' @export
shape_to_log_sl <- function(param) {
  checkmate::assert_number(param, lower = .Machine$double.eps)
  c("log_sl_" = param - 1)
}

#' @rdname Transformations
#' @export
#' @param coef Name of the coefficient.
sl_to_scale <- function(coef) {
  c("scale" = -1/coef)
}

#' @rdname Transformations
#' @export
log_sl_to_shape <- function(coef) {
  c("shape" = coef + 1)
}

#' @rdname Transformations
#' @export
cos_ta_to_kappa <- function(coef) {
  c("kappa" = coef)
}

#' Internal function to create a movement kernel
#' @param template A raster (stack) that is the template for the simulations.
#' @param max.dist The maximum distance of the dispersal kernel
#' @param position The center of the dispersal kernel.
kernel_setup <- function(template, max.dist = 100, position = c(0, 0)) {

  checkmate::assert_class(template, "RasterStack")
  checkmate::assert_number(max.dist, lower = 0)
  checkmate::assert_numeric(position, len = 2)

  p <- sf::st_sf(geom = sf::st_sfc(sf::st_point(position))) %>%
    sf::st_buffer(dist = max.dist)

  # 2. Rasterize buffer
  r1 <- raster::rasterize(p, raster::crop(template, p))

  # 3. Get xy from buffer
  xy <- raster::rasterToPoints(r1)

  # 4. Substract start, so we are centered around 0,0
  xy[, 1] <- xy[, 1] - position[1]
  xy[, 2] <- xy[, 2] - position[2]
  xy[, 1:2]

  k<- tibble(
    x = xy[, 1],
    y = xy[, 2],
    ta_ = get_angle(xy),
    sl_ = sqrt(x^2 + y^2),
    log_sl_ = log(sl_))
  k
}

#' Internal function to shift a movement kernel
#' @param kernel The movement kernel.
#' @param dxy The amount of x and y to shift ther kernel.
kernel_shift <- function(kernel, dxy) {
  checkmate::assert_tibble(kernel)
  checkmate::assert_numeric(dxy, len = 2)

  kernel$x <- kernel$x + dxy[1]
  kernel$y <- kernel$y + dxy[2]
  kernel
}

#' Internal function to rotate a movement kernel
#' @param kernel The movement kernel.
#' @param direction The direction of the rotation. Clockwise in degrees.
kernel_rotate <- function(kernel, direction) {
  checkmate::assert_tibble(kernel)
  checkmate::assert_number(direction)

  kernel$ta_ <- (kernel$ta_ + direction) %% (2 * pi)
  kernel$cos_ta_ <- cos(kernel$ta_)
  kernel
}

#' Internal function to annotate a dispersal kernel with (spatial) covariates
#' @param kernel Ther dispersal kernel
#' @param spatial.covars RasterStack with covariates
#' @param position The current position.
#' @param temporal.covars Additional covariates (non sptatially varying).
kernel_add_covars <- function(kernel, spatial.covars, position, temporal.covars = NULL) {
  checkmate::assert_tibble(kernel)
  checkmate::assert_class(spatial.covars, "RasterStack")
  checkmate::assert_data_frame(temporal.covars, null.ok = TRUE)

  cells <- raster::cellFromXY(spatial.covars, cbind(kernel$x, kernel$y))
  sc1 <- as.data.frame(spatial.covars[cells])

  cc <- raster::cellFromXY(spatial.covars,
                           cbind(position[1], position[2]))
  ii <- which(cc == cells)[1]
  sc2 <- sc1[ii, , drop = FALSE]
  names(sc1) <- paste0(names(sc1), "_end")
  names(sc2) <- paste0(names(sc2), "_start")


  dplyr::bind_cols(kernel[, c("x", "y", "sl_", "ta_", "log_sl_", "cos_ta_")],
                   sc1, sc2, if (!is.null(temporal.covars)) temporal.covars)

}

#' Internal function to "finish" a dispersal kernel.
#' @param kernel The dispersal kernel.
#' @param formula The definition of the kernel from the covariates.
#' @param coefficients The values of the coefficients (need to match formula).
#' @param return.raster Should a raster be returned or not?
#' @param normalize Should the kernel be normaliszed
kernel_finish <- function(
  kernel, formula, coefficients, return.raster = FALSE, normalize = TRUE,
  correct = TRUE) {

  checkmate::assert_tibble(kernel)
  checkmate::assert_formula(formula)
  checkmate::assert_vector(coefficients, names = "unique")
  checkmate::assert_logical(return.raster)

  if (attr(terms(formula), "intercept")) {
    wx <- as.formula(paste("~ 0 + ", as.character(formula)[2]))
  }
  design_matrix <- model.matrix(wx, kernel)

  if (!all(c(colnames(design_matrix) %in% names(coefficients),
             names(coefficients) %in% colnames(design_matrix)))) {
    stop("Some name of coefficients do not match the name of covariates.  In case you have interactions, please make sure that terms within the interactions are given in alphabetical order. E.g., `layer:sl_` instead of `sl_:layer`.")
  }

  # Multiply with coefficients
  kernel$dk <- exp((design_matrix %*% coefficients[colnames(design_matrix)])[, 1])

  if (correct) {
    kernel$dk <- kernel$dk / kernel$sl_
  }
  if (normalize) {
    kernel$dk <- kernel$dk / sum(kernel$dk, na.rm = TRUE)
  }

  if(return.raster) {
    raster::rasterFromXYZ(kernel[, c("x", "y", "dk")])
  } else {
    kernel
  }
}

#' Function to calculate a dispersal kernel, given covariates and coefficients.
#' @param wx A fromula with the selection function.
#' @param coefficients The coefficient value for each term in `wx`.
#' @param start The coordinates of the center of the disspersal kernel.
#' @param spatial.covars A `rasterStack` with the spatial covariates.
#' @param temporal.covars A `tibble` with non spatial covariates.
#' @param direction The direction of the kernel (0 being north) in clockwise direction.
#' @param max.dist Maximum dispersal distance.
#' @param return.raster Should a raster be returned?
#' @param normalize Should the raste be normalized at the end?
#' @export
#' @name simulate

dispersal_kernel <- function(
  wx, coefficients, start,
  spatial.covars, temporal.covars = NULL,
  direction = 0, max.dist = 100, return.raster = FALSE, normalize = TRUE, correct = FALSE) {

  dk0 <- kernel_setup(spatial.covars, max.dist = max.dist, position = start)
  dk1 <- kernel_shift(dk0, dxy = start)
  dk2 <- kernel_rotate(dk1, direction = direction)
  dk3 <- kernel_add_covars(dk2, spatial.covars = spatial.covars,
                           position = start,
                           temporal.covars = temporal.covars)
  dk4 <- kernel_finish(dk3, formula = wx, coefficients = coefficients,
                       return.raster = return.raster, normalize = normalize, correct = correct)

}

#' @rdname simulate
#' @param n Number of time steps.
#' @param as.track Should a track be returned.
#' @param start.time The time stamp of the first simulated location.
#' @param delta.time Time increment for each subsequent step.
#' @export
#'

simulate_track <- function(
  wx, coefficients, start, spatial.covars,
  direction = 0, temporal.covars = NULL, max.dist = 100, n = 10,
  as.track = TRUE, start.time = lubridate::ymd_hms("2020-01-01 00:00:00"),
  delta.time = hours(2)) {

  ##
  if (FALSE) {

    wx = sim.formula
    coefficients = sim.coefs
    start = start1
    spatial.covars = covar.sim
    max.dist = 500
    n = 10
    as.track = TRUE
    delta.time <- hours(2)
    direction = 0
    temporal.covars = NULL
    start <- start
  }
  ###

  checkmate::assert_number(n, lower = 1)
  checkmate::assert_logical(as.track)
  checkmate::assert_class(delta.time, "Period")

  k <- kernel_setup(spatial.covars, max.dist = max.dist, position = start)

  n.px <- nrow(k)
  dir <- direction
  res <- tibble::tibble(x = NA_real_, y = rep(NA_real_, n))
  start1 <- start

  xmn <- raster::xmin(spatial.covars)
  xmx <- raster::xmax(spatial.covars)
  ymn <- raster::ymin(spatial.covars)
  ymx <- raster::xmax(spatial.covars)

  for (i in 1:n) {
    dk <- kernel_shift(k, dxy = start1)
    dk <- kernel_rotate(dk, direction = dir)
    dk <- kernel_add_covars(dk,
                            spatial.covars = spatial.covars, position = start1,
                            temporal.covars = if(!is.null(temporal.covars))
                              temporal.covars[i, , drop = FALSE] else NULL)

    dk <- kernel_finish(dk, formula = wx, coefficients = coefficients,
                        return.raster = FALSE)
    nxt.cell <- wrswoR::sample_int_expj(n.px, size = 1, dk$dk)
    dir <- dir + dk$ta_[nxt.cell]
    start1 <- c(dk$x[nxt.cell], dk$y[nxt.cell])
    res[i, ] <- as.list(start1)

    if ((start1[1] - max.dist) < xmn | (start1[1] + max.dist) > xmx |
        (start1[2] - max.dist) < ymn | (start1[2] + max.dist) > ymx) {
      break()
    }
  }

  trk <- bind_rows(tibble(x = start[1], y = start[2]), res[1:i, ])

  if (as.track) {
    trk$ts <- start.time + 0:i * delta.time
    make_track(trk, x, y, ts)
  } else {
    trk
  }
}

#' @rdname simulate
#' @param n.sim The number of simulations.
#' @export
simulate_track_many <- function(wx, coefficients, start, spatial.covars,
  direction = 0, temporal.covars = NULL, max.dist = 100, n = 10, n.sim = 1,
  as.track = TRUE, start.time = lubridate::ymd_hms("2020-01-01 00:00:00"),
  delta.time = hours(2)) {

  res1 <- lapply(1:n.sim, function(i) {
    xx <- simulate_track(
      wx = wx, coefficients = coefficients, start = start,
      spatial.covars = spatial.covars, direction = direction,
      temporal.covars = temporal.covars, max.dist = max.dist, n = n,
      as.track = as.track, start.time = start.time, delta.time = delta.time)
    xx$rep_ <- i
    xx
  })
  do.call(rbind, res1)
}
