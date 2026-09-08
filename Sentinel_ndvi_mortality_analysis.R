# =============================================================================
# Sentinel NDVI / drought mortality analysis
# -----------------------------------------------------------------------------
#The MODIS analysis revealed the coarse resolution NDVI (250m) was only
#sensitive to changes in canopy greenness as a proxy for 
#canopy mortality above 60% severity. 
#This means that if mortality is more scattered in smaller areas
#it will not be detected.
#The hypothesis is that Sentinel2 data with its much higher 
#resolution (10m) will be more accurate at detecting changing 
#canopy greenness in smaller areas.
#Our research question is:
#"Does the finer spatial resolution of Sentinel2 improve 
#drought induced forest mortality detection?"
#The disadvantage of Sentinel is its shorter temporal resolution
#with a time series from 2017 only (to the present)
#This limits the sample size to around 80 plots - around half of the total

#This script does two things, in order:
#   PART 1: Identifies Sentinel era plots where mortality was observed in 2017 or later
#           These plots are loaded in a csv file and uploaded to GEE          
#   PART 2: Extracts NDVI and EVI Sentinel2 SR Harmonised data for every plot and 
#           combines it with the MODIS database (after monthly aggregation).
#   PART 3: runs the analysis again; i.e. build the "years since mortality" variable, restricts to a
#           -4 to +5 year window around each site's mortality event, fits the
#           mixed model, and summarises the immediate change and recovery
#           time to the control baseline, by forest type and mortality severity.

# Every step is commented. Run the script top to bottom.
# =============================================================================
#lets use GEE
library(rgee)
library(sf)
library(dplyr)

library(rgee)

ee_check() #'tells us which GEE dependencies are missing 
#solution
library(rgee)

ee_install()

ee_clean_pyenv()

ee_install_upgrade()

library(rgee)

ee_Initialize(project='ee1-jennyshep')

ee.Authenticate()

#1. Load sites
library(tidyverse)
library(lubridate)

sites_S2 <- read_csv("C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/Sentinel2/sentinel_sites_v2.csv")

#check
glimpse(sites_s2)

#load required packages
#library(rstac)
#library(terra)
#library(dplyr)
#library(purrr)
#library(tibble)

#function to create annual NDVI composite = median annual 
#Sentinel2 NDVI image for one year
get_ndvi <- function(year){
  
  img <- ee$
    ImageCollection(
      "COPERNICUS/S2_SR_HARMONIZED"
    )$
    filterDate(
      paste0(year,"-01-01"),
      paste0(year,"-12-31")
    )$
    median()
  
  img$
    normalizedDifference(
      c("B8","B4")
    )$
    rename("NDVI")
  
}

#Extract by year
years <- unique(
  sites_s2$mort_year
)

results <- list()

for(y in years){
  
  ndvi <- get_ndvi(y)
  
  pts_year <-
    sites_sf |>
    filter(mort_year == y)
  
  ee_pts <-
    sf_as_ee(pts_year)
  
  out <- ndvi$reduceRegions(
    collection = ee_pts,
    reducer = ee$Reducer$mean(),
    scale = 10
  )
  
  results[[as.character(y)]] <- out
  
}

library(rgee)

ee_check()

library(rgee)

ee_install_upgrade(project='ee1-jennyshep')

#next load the MODIS master NDVI dataset 'master_dataset_cleaned_1'
master <- read_csv(
  "master_dataset_cleaned_1.csv"
)

#merge Sentinel and MODIS dataset
master_final <- master |>
  left_join(
    sentinel_results,
    by = "Plot_ID"
  )

#save combined dataset
write_csv(
  master_final,
  "C:/Users/s422478/OneDrive - Cranfield University/Data/Experimental design/New Data/Sentinel2/master_dataset_with_sentinel.csv"
)
