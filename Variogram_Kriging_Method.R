# 0. Load in Packages

require(sf)
require(terra)
require(gstat)
require(spdep)
require(units)
require(ggplot2)
require(wesanderson)
require(here)

# 1. Load in Data of Interest

ca_cnty  = st_read(here("data", "CA", "CA_Counties", "CA_Counties_TIGER2016.shp"))
ca_ozone = st_read(here("data","CA", "CA_ozone_2017.GPKG"))
ca_ozone = st_transform(ca_ozone, st_crs(ca_cnty))
st_crs(ca_ozone) == st_crs(ca_cnty)

# 2. Create Automated Variogram & Correlogram Generation Function

vgm_krige_analysis <- function(study_area_sf, points_sf, variable_name, nugget = 1e-5, range = 2e5) {
  # Ensure Points and Study Area SF CRS Match
  points_sf <- st_transform(points_sf, st_crs(study_area_sf))
  
  # Create an empirical, exponential, and spherical variogram
  oz_gs <- gstat(formula = as.formula(paste(variable_name, "~ 1")), locations = points_sf)
  vgm_emp <- variogram(oz_gs)
  vgm_mod_exp <- vgm(model = "Exp", nugget = nugget, range = range)
  vgm_mod_sph <- vgm(model = "Sph", nugget = nugget, range = range)
  vgm_fit_exp <- fit.variogram(vgm_emp, vgm_mod_exp)
  vgm_fit_sph <- fit.variogram(vgm_emp, vgm_mod_sph)
  
  # Krige using the Exponential & Spherical Variograms
  temp_rast <- rast(study_area_sf, nrow = 200, ncol = 180)
  temp_rast <- project(temp_rast, st_crs(study_area_sf)$wkt)
  krige_pts <- as.points(temp_rast)
  krige_pts <- st_as_sf(krige_pts)
  krige_pts <- st_intersection(krige_pts, study_area_sf)
  
  AQI_krig_exp <- krige(formula = as.formula(paste(variable_name, "~ 1")), locations = points_sf, newdata = krige_pts, model = vgm_fit_exp)
  AQI_krig_sph <- krige(formula = as.formula(paste(variable_name, "~ 1")), locations = points_sf, newdata = krige_pts, model = vgm_fit_sph)

  return(list(empirical_variogram = vgm_emp, 
              exponential_fit     = vgm_fit_exp,
              spherical_fit       = vgm_fit_sph,
              kriged_maps = list(exp = AQI_krig_exp, sph = AQI_krig_sph),
              krige_elements = list(temp_rast, krige_pts)))
}

vgm_krige_results = vgm_krige_analysis(study_area_sf = ca_cnty, 
                                       points_sf = ca_ozone, 
                                       variable_name = "AQI")

# 3. Retrieve Variogram Krige Outputs

# Extract the empirical variogram
vgm_emp <- vgm_krige_results$empirical_variogram

# Extract the exponential and spherical model fits
vgm_fit_exp <- vgm_krige_results$exponential_fit
vgm_fit_sph <- vgm_krige_results$spherical_fit

# Plot the empirical variogram with the exponential & spherical model fit
plot(vgm_emp, pch = 16, cex = 1.2, main = "AQI Empirical Variogram")
plot(vgm_emp, vgm_fit_exp, main = "AQI Exponential Variogram")
plot(vgm_emp, vgm_fit_sph, main = "AQI Spherical Variogram")

# Extract the kriged maps
AQI_krig_exp <- vgm_krige_results$kriged_maps$exp
AQI_krig_sph <- vgm_krige_results$kriged_maps$sph

# 4. Plot the kriged maps

ggplot() +
  geom_sf(data = ca_cnty, fill = "white", color = "black", size = 0.2, show.legend = FALSE) + 
  geom_sf(data = AQI_krig_exp, aes(color = var1.pred), size = 0.5, show.legend = "point") +  
  scale_color_gradientn(colors = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) +  
  labs(title = "Kriged AQI Map (Exponential Model)", x = "", y = "", color = "AQI") +
  theme_minimal() +  
  theme(legend.position = "right") 

ggplot() +
  geom_sf(data = ca_cnty, fill = "white", color = "black", size = 0.2, show.legend = FALSE) +
  geom_sf(data = AQI_krig_sph, aes(color = var1.pred), size = 0.5, show.legend = "point") +
  scale_color_gradientn(colors = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) +  
  labs(title = "Kriged AQI Map (Spherical Model)", x = "", y = "", color = "AQI") +
  theme_minimal() + 
  theme(legend.position = "right") 
