
############################################################################################################################
############################################################################################################################
###########      
###########                     Regional changes in temperature-related mortality before and during the COVID-19 pandemic: 
###########                                  a continental modelling analysis in 805 European regions
###########
###########                             Blanca Paniello-Castillo (blanca.paniello@isglobal.org)
###########
###########                                                     June 2025
###########
#############################################################################################################################
#############################################################################################################################


COVID_adjustment <- function(df8, dfmax, wave_date1, wave_date2, crossbasis_fun, lag_fun, datatable, threshold = 5) {
  
  # Load required packages using requireNamespace
  pkgs <- c("lubridate", "tidyverse", "dlnm", "splines", "mixmeta", "DistributionUtils", 
            "tsModel", "MASS", "ISOweek", "giscoR", "sf", "grid", "writexl", "scales", "gridExtra")
  sapply(pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop(paste("Package", pkg, "is required but not installed."))
  })
  
  # Error checks
  if (!is.data.frame(datatable_cali)) stop("`datatable_cali` must be a data frame.")
  if (!is.numeric(df8) || !is.numeric(dfmax)) stop("`df8` and `dfmax` must be numeric.")
  if (!is.Date(wave_date1) || !is.Date(wave_date2)) stop("Wave dates must be Date objects.")
  if (df8 > dfmax) stop("`df8` must be less than `dfmax`.")
  
  # Ensure only one region is being processed
  iREG <- unique(datatable_cali$location)
  if (length(iREG) != 1) stop("Function currently supports only one region per call.")
  
  # Prepare data list (indexed by region)
  datalist_cali <- list()
  datalist_cali[[iREG]] <- datatable_cali
  datalist_cali[[iREG]]$wop <- 1:length(datatable_cali$woy)
    
  # Initialize output containers
  results_df_list <- list()
  
  df_change_value <- df8
  perc_tot <- 100
  
  while (df_change_value <= dfmax && perc_tot > threshold) {
    
    # Split data into before, during, and after wave (define wave periods)
    vBEFORE_DATE <- datalist_cali[[iREG]]$date[datalist_cali[[iREG]]$date <= wave_date1]
    vDURING_DATE <- datalist_cali[[iREG]]$date[datalist_cali[[iREG]]$date >= wave_date1 & datalist_cali[[iREG]]$date <= wave_date2]
    vAFTER_DATE <- datalist_cali[[iREG]]$date[datalist_cali[[iREG]]$date >= wave_date2]
    
    # Create spline knots of the seasonality term proportionally based on df 
    vBEFORE_DATE <- vBEFORE_DATE[seq(from = 1, to = length(vBEFORE_DATE), length.out = round(df8 * length(vBEFORE_DATE) * 7 / 365.25))]
    vDURING_DATE <- vDURING_DATE[seq(from = 1, to = length(vDURING_DATE), length.out = round(df_change_value * length(vDURING_DATE) * 7 / 365.25))]
    vAFTER_DATE <- vAFTER_DATE[seq(from = 1, to = length(vAFTER_DATE), length.out = round(df8 * length(vAFTER_DATE) * 7 / 365.25))]
    
    # Create complete list of knots
    vKNOTS <- datalist_cali[[iREG]]$wop[datalist_cali[[iREG]]$date %in% unique(c(vBEFORE_DATE, vDURING_DATE, vAFTER_DATE))]
    vKNOTS <- vKNOTS[2:(length(vKNOTS)-1)]  # remove boundary knots
    
    # Models
    FORMULA_SEA <- mort ~ splines::ns(wop, knots = vKNOTS)
    FORMULA_CRB <- mort ~ splines::ns(wop, knots = vKNOTS) + CROSS_BASIS
    
    # Fit seasonality model
    GLM_MODEL_SEA <- glm(FORMULA_SEA, datalist_cali[[iREG]], family = quasipoisson, na.action = "na.exclude")
    datalist_cali[[iREG]]$mort_pred_seas <- predict(GLM_MODEL_SEA, type = "response")
    
    # Build crossbasis object (DLNM)
    VAR_PRC <- c(0.1, 0.5, 0.9)
    MIN_LAG <- 0
    MAX_LAG <- 3
    
    CROSS_BASIS <- dlnm::crossbasis(
      datalist_cali[[iREG]]$temp, c(MIN_LAG, MAX_LAG),
      argvar = list(fun = crossbasis_fun, knots = quantile(datalist_cali[[iREG]]$temp, VAR_PRC, na.rm = TRUE)),
      arglag = list(fun = lag_fun)
    )
      
    # Fit model with seasonality + temperature crossbasis
    GLM_MODEL_CRB <- glm(FORMULA_CRB, datalist_cali[[iREG]], family = quasipoisson, na.action = "na.exclude")
    datalist_cali[[iREG]]$mort_pred_seas_cb <- predict(GLM_MODEL_CRB, type = "response")
    
    # Residual evaluation
    residuals <- datalist_cali[[iREG]]$mort - datalist_cali[[iREG]]$mort_pred_seas_cb # This can be changed for seasonality only
    ci_bounds <- quantile(residuals, probs = c(0.025, 0.975), na.rm = TRUE)
    residuals_filtered <- residuals[datalist_cali[[iREG]]$date >= wave_date1 & datalist_cali[[iREG]]$date <= wave_date2]
    
    # Calculate % of residuals outside CI
    peaks_outside_ci1 <- sum(residuals_filtered > ci_bounds[2] | residuals_filtered < ci_bounds[1])
    peaks_total <- length(residuals_filtered)
    perc_tot <- (peaks_outside_ci1 / peaks_total) * 100
    
    # Save results
    results_df_list[[as.character(df_change_value)]] <- data.frame(
      Region_Name = iREG,
      df_change = df_change_value,
      peaks_outside_ci1 = peaks_outside_ci1,
      peaks_total = peaks_total,
      perc_tot = perc_tot
    )
    
    # Stop if threshold met
    if (perc_tot <= threshold) break  # stop when ≤5%
    
    # Increment df for next iteration
    df_change_value <- df_change_value + 1
  }
  
  # Return full results
  results_df <- do.call(rbind, results_df_list)
  
  # If threshold never met, pick best df (min perc_tot)
  if (perc_tot > threshold) {
    best_index <- which.min(results_df$perc_tot)
    df_change_value <- results_df$df_change[best_index]
  }
  
  return(list(final_df = df_change_value, results_df = results_df))
}
