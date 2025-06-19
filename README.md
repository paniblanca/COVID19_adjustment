# COVID-19 pandemic adjustment and temperature-related mortality
This function incrementally adjusts the degrees of freedom (df) of a cross-basis model during the COVID-19 pandemic period to ensure the percentage of residuals outside the 95% confidence interval during the wave period is below a defined threshold.

More concretely, we selected specific time windows within identified COVID-19 waves that showed large residuals, i.e., those exceeding the 2.5th or 97.5th centiles of the residuals. 
The residuals were calculated by subtracting the all-cause mortality counts of the full model including the temperature cross-basis.
Departing from the standard choice of 8 d.f. per year, we applied equally-spaced knots within these time windows, increasing the d.f. as necessary to reduce the negative-positive-negative mortality residuals below a defined threshold. We imposed that each COVID-19 wave window had a maximum of 5 % of mortality residual values outside the 2.5–97.5 residual threshold range.

## Getting Started
The function necessitates from the specification of 7 parameters. Moreover, the input of data also requires from a specific structure. 
The function is designed to work with weekly data, as used in the provided example. However, it can be easily adapted for daily data or other temporal resolutions.

### Parameters 
What you need to specify in the function (format and explanation). 
* <ins>df8</ins>: Numeric. Usually the standard choice of 8 df per year
* <ins>dfmax</ins>: Numeric. Maximum df allowed for testing during the wave. The maximum will depend on the data resolution (52 df for weekly data).
* <ins>wave_date1</ins>: Date. Start date of the adjusting window. 
* <ins>wave_date2</ins>: Date. End date of the adjusting window. 
* <ins>crossbasis_fun</ins>: Character. Spline function for the crossbasis. Either natural cubic spline (*ns*) or B-spline (*bs*).
* <ins>lag_fun</ins>: Character. Internal spline function for generating the lad matrix. Here we use *integer* but it could also be *poly*, *strata*, *thr* or *lin*.
* <ins>datatable</ins>: Data frame. Input data for a **single region**, containing variables `location` `mort`, `temp`, and `date`, among others (see Input Data Requirements section).
* <ins>threshold</ins>: Numeric. Acceptable percentage of out-of-bound residuals during the wave. We imposed that each window had a maximum of 5 % of mortality residual values outside the 2.5–97.5 residual threshold range. 

```r
Example usage (Lisbon region):

  df8 = 8
  dfmax = 52
  wave_date1 = as.Date("2020-10-29")
  wave_date2 = as.Date("2021-03-25")
  crossbasis_fun = "ns"
  lag_fun = "integer"
  datatable = datatable_data
  threshold = 5
```
###  Input Data Requirements (datatable)
The input dataframe must include the following columns:
* <ins>location</ins>: region identifier (only one region per function call is allowed)
* <ins>year</ins>: year
* <ins>woy</ins>: week of the year
* <ins>mortality</ins>: daily or weekly all-cause mortality counts
* <ins>temperature</ins>: corresponding temperature values
* <ins>date</ins>: data in Date format
  
Important:
The calibration dataset (datatable) should extend beyond the last calibration date by MAX_LAG * 7 days (e.g., 21 days for a maximum lag of 3 weeks). This is necessary to properly model lagged temperature effects in the cross-basis.

### Output
The function returns a list with:
* final_df: Optimal df used during wave that satisfies threshold.
* results_df: Data frame with percentage of residuals outside CI for each df tested.

# Reference: 
Paniello-Castillo B, Quijal-Zamorano M, Gallo E, Basagaña X, Ballester J. Regional changes in temperature-related mortality before and during the COVID-19 pandemic: a continental modelling analysis in 805 European regions. Environ Res [Internet]. 2025;278(121697):121697. Available from: http://dx.doi.org/10.1016/j.envres.2025.121697

