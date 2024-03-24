# 0. Load in dependency Correlogram_Function.R
source("path/to/Correlogram_Function.R")

# 0. Load in Packages

require(sf)
require(terra)
require(gstat)
require(spdep)
require(units)
require(ggplot2)
require(wesanderson)
require(here)

# 1. Load in data & Inspect Dataset

ca_cnty  = st_read(here("data", "CA", "CA_Counties", "CA_Counties_TIGER2016.shp"))
ca_ozone = st_read(here("data","CA", "CA_ozone_2017.GPKG"))
ca_ozone = st_transform(ca_ozone, st_crs(ca_cnty))

# 2. Inspect the residual error of ozone & AQI

# 2.1: Create your linear model with lm()
aqi_lm = lm(ozone ~ AQI, data = ca_ozone)

# 2.2: Inspect the residuals by looking at the fitted residuals & QQPlot
par(mfrow = c(1,2))
plot(aqi_lm, which = 1)
plot(aqi_lm, which = 2)

# 2.3: Extract residuals to ca_ozone from the model fit!
ca_ozone$resids = residuals(aqi_lm)

correlogram_calculator <- function(data, variable_name = "ozone", n_dist_class = 10, nsim = 999) {
  # 1. Moran's Residuals Correlogram
  distmat <- st_distance(data)
  maxdist <- 0.5 * max(distmat)
  mindist <- min(distmat)
  dist_classes <- seq(mindist, maxdist, length.out = n_dist_class)
  
  # Calculate Euclidean distance between neighbors
  neighbors <- dnearneigh(
    x = data,
    d1 = as.numeric(dist_classes[1]),
    d2 = as.numeric(dist_classes[2]),
    longlat = FALSE)
  
  # Make a weights object using our euclidean distances
  wts <- nb2listw(neighbors, style = 'W', zero.policy = TRUE)
  
  # Create a Moran's I test with normal approximation vs Monte Carlo permutation test
  mor_mc <- moran.mc(
    data[[variable_name]],
    listw = wts,
    nsim = nsim,
    zero.policy = TRUE)
  
  # Moran's I Normal Distribution Test
  mor_norm <- moran.test(
    data[[variable_name]],
    listw = wts,
    randomisation = FALSE,
    zero.policy = TRUE)
  
  # 2. Generate a Distance Correlogram
  correlog <- calculate_morans_i(
    data = data,
    variable_name = variable_name,
    mindist = mindist,
    maxdist = maxdist,
    n_dist_class = n_dist_class,
    nsim = nsim)
  
  return(list(mor_mc = mor_mc, mor_norm = mor_norm, correlogram = correlog))
}

ozone_correlogram_results = correlogram_calculator(ca_ozone, variable_name = "ozone")
resid_correlogram_results = correlogram_calculator(ca_ozone, variable_name = "resids")

# 3. Plot Correlograms!

# Extract the correlogram data for "ozone" & plot
ozone_correlogram <- ozone_correlogram_results$correlogram
ggplot(ozone_correlogram, aes(x = dist, y = Morans.i)) +
  geom_smooth(se = FALSE) +
  geom_point() +
  geom_ribbon(aes(ymin = Null.lcl, ymax = Null.ucl), fill = adjustcolor("steelblue", 0.2)) +
  ggtitle("CA Ozone Correlogram") +
  theme_minimal()

# Extract the correlogram data for "resids" & plot
resids_correlogram <- resid_correlogram_results$correlogram
ggplot(resids_correlogram, aes(x = dist, y = Morans.i)) +
  geom_smooth(se = FALSE) +
  geom_point() +
  geom_ribbon(aes(ymin = Null.lcl, ymax = Null.ucl), fill = adjustcolor("steelblue", 0.2)) +
  ggtitle("CA Ozone ~ AQI Residual Correlogram") +
  theme_minimal()
