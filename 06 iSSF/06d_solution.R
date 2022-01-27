library(amt)

# 1. Load the `amt_fisher` dataset (it is shipeed with `amt` and you can get it with `data(amt_fisher)` or read it from the csv file `data/fisher.csv`). If you're working with your own data, format it as a `track_xyt` object using `amt::make_track()`.

data("amt_fisher")
data("amt_fisher_covar")

fisher <- make_track(amt_fisher, x_, y_, t_, name = name)

# 2. Filter all data points for  `Leroy` and resample the data to 30 min. 

leroy <- filter(fisher, name == "Leroy")
summarize_sampling_rate(leroy)
leroy <- track_resample(leroy, rate = minutes(30), tolerance = minutes(2))

print(leroy, n = 100)
   
# 3. Use `Leroy` and create 15 random steps for each observed step. 

leroy <- leroy %>% steps_by_burst() %>% random_steps(n_control = 15)


# 4. Extract covariates at the end of each step for elevation. A raster with the elevation data is also shipped with `amt` and you can get the elevation by executing `data(amt_fisher_covar)`.  This returns a list with three covariates, we are only interested in elevation. To access elevation you can use: `amt_fisher$elevation`.
leroy <- leroy %>% extract_covariates(amt_fisher_covar$elevation)

# 5. Fit a SSF an SSF, where you use elevation as the only covariate.
m1 <- leroy %>% fit_ssf(case_ ~ elevation + strata(step_id_), model = TRUE)
summary(m1)

# 6. What is the log-RSS for a step ending at different elevations. Choose a suitable range of elevations.
range(leroy$elevation, na.rm = TRUE)
s1 <- data.frame(
  elevation = 80:120
)
s2 <- data.frame(
  elevation = 100
)
logrss <- log_rss(m1, s1, s2, ci = "se")

logrss$df
ggplot(logrss$df, aes(x = elevation_x1, y = log_rss)) + geom_line()
ggplot(logrss$df, aes(x = elevation_x1, y = log_rss)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line() 
    
ggplot(logrss$df, aes(x = elevation_x1, y = exp(log_rss))) + geom_line()

# 7. Finally fit an `iSSF` where you include in addition to `elevation` also the step length the log of the step length and the cosine of the turn angle.
m2 <- leroy %>% fit_ssf(case_ ~ elevation + sl_ + log(sl_) + cos(ta_) + strata(step_id_))
summary(m2)
