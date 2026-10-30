# 2026_whg.27.6a_assessment
2026 - Whiting (Merlangius merlangus) in Division 6.a (West of Scotland) - WGCSE (ADGCS)

The 'Data' folder contains everything to run the SAM assessment.

This year (2026) the SSB is estimated to be below Blim. As a result the forecast was run with a low productivity scenario. To do so the median recruitment used in the forecast was resampled from the last 6 years (2020-2025) instead of the last 10 years as this is more representative of the recent low recruitment estimates. 

Run the R scripts in this order:
1- model
2- forecast_low_productivity
3- forecast_summary_low_productivity
4- forecastExploration_productivity
5- prob_above_Blim_low_prod

The script 'forecast_sam_with_natage_v3' is a function called in the scripts above.
