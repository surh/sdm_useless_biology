library(tidyverse)

#' Read data downloaded from:
#' https://raw.githubusercontent.com/jcoliver/learn-r/gh-pages/data/Carnegiea-gigantea-GBIF.csv
Obs <- read_csv("data/teosintle_maxent_input.csv") %>%
  rename(latitude = lat, longitude = lon)
Obs

#' Dowload bioclimatic data. Can take some time (627.9 MB)
bioclim_data <- geodata::worldclim_global(var = "bio",
                                          res = 2.5,
                                          path = "data/")

#' Download world map. Relatively fast.
world_map <- geodata::world(resolution = 3,
                            path = "data/")

#' Find quadrant that eoncompasses species distribution, increase each
#' size 1.5 fold. Save data to raster object and save objects for future use.
max_lat <- ceiling(max(Obs$latitude))
min_lat <- floor(min(Obs$latitude))
max_lon <- ceiling(max(Obs$longitude))
min_lon <- floor(min(Obs$longitude))
quadrant <- terra::ext(x = c(min_lon, max_lon, min_lat, max_lat))
quadrant <- quadrant * 1.5

quadrant_map <- terra::crop(x = world_map, y = quadrant)
bioclim_data <- terra::crop(x = bioclim_data, y = quadrant)

terra::writeRaster(bioclim_data, filename = "data/teosintle_bioclim_raster.tif")
terra::writeVector(quadrant_map, filename = "data/teosintle_map")


#' # Scenario data
#' Download proyected bioclimatic variables under different climate scenarios
#` from the Coupled Model Intercomparison Project (CMIP) 6
forecast_data <- geodata::cmip6_world(model = "MPI-ESM1-2-HR",
                                      ssp = "245",
                                      time = "2081-2100",
                                      var = "bioc",
                                      res = 2.5,
                                      path = "data")
forecast_data <- terra::crop(forecast_data, quadrant)
names(forecast_data) <- names(bioclim_data)
terra::writeRaster(forecast_data,
  filename = "data/teosintle_forecast_2081-2100_raster.tif")


forecast_data <- geodata::cmip6_world(model = "MPI-ESM1-2-HR",
                                      ssp = "245",
                                      time = "2061-2080",
                                      var = "bioc",
                                      res = 2.5,
                                      path = "data")
forecast_data <- terra::crop(forecast_data, quadrant)
names(forecast_data) <- names(bioclim_data)
terra::writeRaster(forecast_data,
  filename = "data/teosintle_forecast_2061-2080_raster.tif")


sessionInfo()