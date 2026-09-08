# Jen-PhD
All my PhD work involving analysis of time series (MODIS) satellite data searching for signals of resilience and critical slowing down (may be the same thing?) in (UK) forests of relevance to Europe and beyond.
Coding work is being done in Python using Colab and Jupyter Lab (plus some R for mixed models and some visualisations) 
The raw NDVI data is extracted from Google Earth Engine, specifically MODIS (MOD13Q1) Version 6.1 and Harmonised Landsat and Sentinel-2 (HLS). 
Noise reduction has been achieved using Savitsky-Golay filter, Seasonal Trend decomposition using Loess (STL). And during the original DLM work, a smoothing/declouding algorithm called SWEEP and subsequently SWEEP2.
Further analysis will be conducted using autocorrelation method(s) 
Subsequent to the DLM work there has been analysis of the relationship between NDVI values and forest type and mortality severity (other independent variables are explored e.g. climate). This includes threshold analysis to determine the minimal mortality necessary for detection using MODIS.  
Final experimental work was a rigorous comparison of MODIS and HLS NDVI including processes such as auditing the dataset for data quality, computing annual means, plotting comparative trajectories, statistical evaluation and coding some complex static and animated visualisations.
