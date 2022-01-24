#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#---------------Module 02 -- Home Ranges---------------X
#----------------Last updated 2021-01-06---------------X
#-------------------Exercise Solution------------------X
#######################################################X

# Using amt_fisher data

# Load packages ----
library(tidyverse)
library(amt)

# Part 1 ----
# Load in your data or the `amt_fisher` dataset. If you're working with your 
# own data, format it as a `track_xyt` object using `amt::make_track()`.
dat <- amt_fisher

# Part 2 ----
# Create a nested `data.frame` by nesting by individual ID (or burst, if you 
# have separate bursts in your data).
nd <- dat %>% 
  dplyr::select(id, name, sex, everything()) %>% 
  nest(track = x_:t_)

# You don't need to keep all 3 id columns. For example, this is fine:
dat %>% 
  nest(track = -id)

# Part 3 ----  
# Fit a home range of your choice for each individual (feel free to fit more 
# than one home range type). Fit at least two isopleths per home range, *e.g.*, 
# 95% and 50%. 
nd <- nd %>% 
  mutate(mcp = lapply(track, hr_mcp, levels = c(0.5, 0.95)),
         locoh = map(track, hr_locoh, n = 11, levels = seq(0.25, 1, by = 0.25)),
         kde = map(track, hr_kde, levels = c(0.5, 0.95)))
 
# Part 4 ----   
# Compute the area of each home range polygon. Convert the area to reasonable 
# units for your study animal, such as ha or km^2.

# Function to grab area out of hr_area list
get_area <- function(l, level) {
  ll <- lapply(l, function(x) {
    x$area[which(x$level == level)]
  })
  
  res <- do.call(c, ll)
  return(res)
}

(ar <- nd %>% 
  mutate(mcp_area = map(mcp, hr_area, units = TRUE),
         locoh_area = map(locoh, hr_area, units = TRUE),
         kde_area = map(kde, hr_area, units = TRUE)) %>% 
  mutate(mcp95 = get_area(mcp_area, 0.95),
         mcp50 = get_area(mcp_area, 0.5),
         locoh100 = get_area(locoh_area, 1),
         locoh75 = get_area(locoh_area, 0.75),
         locoh50 = get_area(locoh_area, 0.50),
         locoh25 = get_area(locoh_area, 0.25),
         kde95 = get_area(kde_area, 0.95),
         kde50 = get_area(kde_area, 0.5)) %>% 
  select(id:sex, mcp95:kde50) %>% 
  mutate(across(mcp95:kde50, units::set_units, "km^2")))

# Here's a solution from Christina Aiello
a_mcp <- map_df(nd$mcp, hr_area, units = TRUE)

# To add the names back, you can tweak this solution to this:
# Return a list rather than a single data.frame
a_mcp <- map(nd$mcp, hr_area, units = TRUE)
# Give the list elements names
names(a_mcp) <- nd$id
# Add the argument '.id' to 'bind_rows()' to keep the list element names as
# a new column
bind_rows(a_mcp, .id = "ID")

# Here's another solution with unnest.
a <-  nd %>% 
  mutate(mcp_area = map(mcp, hr_area, units = TRUE))

a %>% 
  dplyr::select(id, mcp_area) %>% 
  unnest(cols = mcp_area)


# Part 5 ---- 
# Make a map with each individual's home range. It's up to you whether you 
# want all your individuals on a single map or one map per individual.

# Polygons
poly <- nd %>% 
  mutate(mcp_poly = map(mcp, hr_isopleths),
         locoh_poly = map(locoh, hr_isopleths),
         kde_poly = map(kde, hr_isopleths)) %>% 
  select(id:track, mcp_poly:kde_poly)

# Plot for each individual
plots <- lapply(poly$name, function(nm) {
  # Pull out polygons
  mcp <- poly$mcp_poly[[which(poly$name == nm)]]
  locoh <- poly$locoh_poly[[which(poly$name == nm)]]
  kde <- poly$kde_poly[[which(poly$name == nm)]]
  
  # Combine in 1 data.frame
  polys <- rbind("MCP" = mcp,
                     "LoCoH" = locoh,
                     "KDE" = kde)
  # Get HR type from row.names
  polys$HR <- unlist(lapply(
    strsplit(row.names(polys), ".", fixed = TRUE), 
    getElement, 1))
  
  # Get the points
  trk <- poly$track[[which(poly$name == nm)]]
  
  g <- ggplot() +
    geom_sf(data = polys, fill = NA,
            aes(color = HR, linetype = factor(level))) +
    geom_point(data = trk, aes(x = x_, y = y_), size = 1,
               shape = 21, color = "#00000055", fill = "#55555555") +
    scale_linetype_discrete(name = "Isopleth") +
    ggtitle(nm) +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw()
  
  print(g)
  return(g)
})


# Here's a simpler solution.
nd %>% 
  pull(mcp) %>% 
  map(plot)

# Here's an intermediate solution
nd %>% 
  pull(mcp) %>% 
  map(function(x) {
    hr_isopleths(x) %>% 
      ggplot() +
      geom_sf(aes(color = factor(level)), size = 1, fill = NA) +
      theme_bw()
  })
  