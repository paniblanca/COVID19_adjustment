
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

# Load the function
source("COVID19_adjustment_function.R")

# Load libraries
library(lubridate) ; # Required Functions: wday

# Read data
datatable_data <- read.csv("./data_lisbon.csv")
datatable_data$date <- as.Date(datatable_data$date)

# Calibration date range
date1_cali <- as.Date("2020-01-02") # start
date2_cali <- as.Date("2023-09-21") # end

# Validate dates fall on Thursdays (wday == 4 with week starting on Monday)
if(wday(date1_cali, week_start = 1) != 4 | wday(date2_cali, week_start = 1) != 4){
  stop("ERROR: Invalid Calibration Dates !!!")}

# Filter calibration data (Account for maximum lag (3 weeks) and convert lag to days (MAX_LAG * 7) used in the model for weekly data)
MAX_LAG <- 3            # maximum lag in weeks
datatable_cali = datatable_data[which(date1_cali <= datatable_data$date & datatable_data$date <= date2_cali + max(7 * MAX_LAG, 0)),]

# Run function
result <- COVID_adjustment(
  df8 = 8,
  dfmax = 52, 
  wave_date1 = as.Date("2020-10-29"),
  wave_date2 = as.Date("2021-03-25"),
  crossbasis_fun = "ns",
  lag_fun = "integer",
  datatable = datatable_cali,
  threshold = 5
)

print(result$final_df)      # Final df value that satisfied condition
print(result$results_df)    # All tested df with corresponding stats

