calculate_morans_i <- function(data, variable_name, mindist, maxdist, n_dist_class = 10, nsim = 999) {
  # Make a sequence of distance classes between the min and max distances
  dist_classes <- seq(mindist, maxdist, length.out = n_dist_class)
  dist_classes <- as.numeric(dist_classes)
  
  # Create an empty data frame for storing output
  correlog_result <- data.frame(
    dist = dist_classes[-n_dist_class],
    Morans.i = NA,
    Null.lcl = NA,
    Null.ucl = NA,
    Pvalue = NA
  )
  
  # Use a for loop to calculate Moran's I for lag distances
  for (i in 1:(length(dist_classes) - 1)) {
    d.start <- dist_classes[i] 
    d.end   <- dist_classes[i + 1]
    
    # dnearneigh() accepts a sp or sf object as its x argument
    neigh_i <- dnearneigh(x = data, d1 = d.start, d2 = d.end)
    
    # Make a weights object
    wts_i <- nb2listw(
      neighbours = neigh_i,
      style = 'W',
      zero.policy = TRUE
    )
    
    mor.i <- moran.mc(data[[variable_name]], 
                      listw = wts_i,
                      nsim = nsim,
                      zero.policy = TRUE)
    
    # Summarize results
    correlog_result[i, "Morans.i"] = mor.i$statistic
    correlog_result[i, "Null.lcl"] = quantile(mor.i$res, probs = 0.025, na.rm = TRUE)
    correlog_result[i, "Null.ucl"] = quantile(mor.i$res, probs = 0.975, na.rm = TRUE)
    correlog_result[i, "Pvalue"]   = mor.i$p.value
  }
  
  return(correlog_result)
}