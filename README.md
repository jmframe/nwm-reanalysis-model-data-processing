# Here are scripts to process the NOAA National Water Model Reanalysis Model Data for daily averages at the camels catchments.
# The scripts here require data described here: https://docs.opendata.aws/nwm-archive/readme.html
# Either NWM version 1.2 or 2.0 will work.
# There are three components to the NWM that are processed here:
## 1) RTOUT: Geospatial, 250m Gridded NetCDF
## 2) CHRTOUT: Point Type (including Reach ID)
## 3) LDASOUT: Geospatial, 1Km Gridded NetCDF
## These data should go in the appropriate sub-directory (v1 or v2)
# There is a specific order that these scripts should be run:
## 1) collate CHRT
## 2) ...
