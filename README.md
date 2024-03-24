# Variograms & Kriging Workflow
## 1. Load in SF Points & SF Study Area
```r
ca_cnty  = st_read(here("data", "CA", "CA_Counties", "CA_Counties_TIGER2016.shp"))
ca_ozone = st_read(here("data","CA", "CA_ozone_2017.GPKG"))
ca_ozone = st_transform(ca_ozone, st_crs(ca_cnty))
st_crs(ca_ozone) == st_crs(ca_cnty)
```
## 2. Utilize the vgm_krige_analysis() function
This custom function wraps together all the essential components to generate a krige. You can specify your study area, points, and nugget/range for your krige here.
```r
vgm_krige_results = vgm_krige_analysis(study_area_sf = ca_cnty, 
                                       points_sf = ca_ozone, 
                                       variable_name = "AQI")
```
## 3-4. Extract Variograms & Kriged Outputs for Plotting 
![KrigeEXP](https://raw.githubusercontent.com/JTSALAH/R-Spatial-Autocorrelation-Analytics/main/Example_Output/Kriged_AQI_EXP.png)
![KrigeSPH](https://raw.githubusercontent.com/JTSALAH/R-Spatial-Autocorrelation-Analytics/main/Example_Output/Kriged_AQI_SPH.png)

# Correlogram Workflow
## 0. Load in dependency custom function from Correlogram_Function.R
Correlogram_Function.R contains the calculate_morans_i() custom function which loops through unique distance classes and calculates Moran's I at those distances. Confidence intervals are additionally calculated in order to generate an envelope and determine the presence/abscence of autocorrelation in our data.
```r
source("path/to/Correlogram_Function.R")
```

## 1. Load in SF Points & SF Study Area
```r
ca_cnty  = st_read(here("data", "CA", "CA_Counties", "CA_Counties_TIGER2016.shp"))
ca_ozone = st_read(here("data","CA", "CA_ozone_2017.GPKG"))
ca_ozone = st_transform(ca_ozone, st_crs(ca_cnty))
st_crs(ca_ozone) == st_crs(ca_cnty)
```

## 2. Create & Inspect a Linear Model for our data 
```r
# 2.1: Create your linear model with lm()
aqi_lm = lm(ozone ~ AQI, data = ca_ozone)

# 2.2: Inspect the residuals by looking at the fitted residuals & QQPlot
par(mfrow = c(1,2))
plot(aqi_lm, which = 1)
plot(aqi_lm, which = 2)

# 2.3: Extract residuals to ca_ozone from the model fit!
ca_ozone$resids = residuals(aqi_lm)
```

## 3. Run correlogram_calculator() custom function 
```r
ozone_correlogram_results = correlogram_calculator(ca_ozone, variable_name = "ozone")
resid_correlogram_results = correlogram_calculator(ca_ozone, variable_name = "resids")
```
## 4. Plot Correlograms for your variable of interest & residuals!
![Correlog](https://raw.githubusercontent.com/JTSALAH/R-Spatial-Autocorrelation-Analytics/main/Example_Output/CA_Ozone_Correlogram.png)
![ResidCorrelog](https://raw.githubusercontent.com/JTSALAH/R-Spatial-Autocorrelation-Analytics/main/Example_Output/CA_Ozone_AQI_Residual_Correlogram.png)
