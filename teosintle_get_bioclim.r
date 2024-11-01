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
#' side of the quadrant 1.5 fold. 
max_lat <- ceiling(max(Obs$latitude))
min_lat <- floor(min(Obs$latitude))
max_lon <- ceiling(max(Obs$longitude))
min_lon <- floor(min(Obs$longitude))
quadrant <- terra::ext(x = c(min_lon, max_lon, min_lat, max_lat))
quadrant <- quadrant * 1.5

#' Crop bioclimatic and map data to the quadrant
quadrant_map <- terra::crop(x = world_map, y = quadrant)
bioclim_data <- terra::crop(x = bioclim_data, y = quadrant)

#' Saved cropped data for future use. These files are the input
#' for analysis in the class.
terra::writeRaster(bioclim_data, filename = "data/teosintle_bioclim_raster2.tif")
terra::writeVector(quadrant_map, filename = "data/teosintle_map")


#' # Scenario data
#' Download proyected bioclimatic variables frin different models under
#' different Shared Socioeconomic Pathways (SSP) from the Coupled Model
#' Intercomparison Project (CMIP) 6.
#' Takes a long time to download (100s of MB per model and pathway combination)
forecast_data <- geodata::cmip6_world(model = "MPI-ESM1-2-HR",
                                      ssp = "245",
                                      time = "2081-2100",
                                      var = "bioc",
                                      res = 2.5,
                                      path = "data")

#' Crop forecast data to the quadrant and save for future use                                    
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


#' List all  of models, pathways, and time periods to
#' download forecast data
model <- c("ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", "BCC-CSM2-MR",
           "CanESM5", "CanESM5-CanOE", "CMCC-ESM2", "CNRM-CM6-1",
           "CNRM-CM6-1-HR", "CNRM-ESM2-1", "EC-Earth3-Veg", "EC-Earth3-Veg-LR",
           "FIO-ESM-2-0", "GFDL-ESM4", "GISS-E2-1-G", "GISS-E2-1-H",
           "HadGEM3-GC31-LL", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR",
           "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR",
           "MRI-ESM2-0", "UKESM1-0-LL")
ssp <- c("126", "245", "370", "585")
time <- c("2041-2060", "2061-2080", "2081-2100")



#' Download all forecast data, for all combinations of models, pathways and
#' timeframes. Crop each forecast the species quadrant, and save for future use.
#' This will take a long time to run (100s of MB per combination)
expand.grid(model = model, ssp = ssp, time = time) %>%
  mutate(path = paste0("data/teosintle_forecast_",
                       time, "_", ssp, "_", model, ".tif")) %>%
  rowwise() %>%
  mutate(forecast = list(geodata::cmip6_world(model = model,
                                              ssp = ssp,
                                              time = time,
                                              var = "bioc",
                                              res = 2.5,
                                              path = "data")))  %>%
  mutate(forecast = list(terra::crop(forecast, quadrant))) %>%
  mutate(forecast = list(terra::writeRaster(forecast, filename = path))) %>%
  ungroup()








forecasts %>% 
  mutate(path = paste0("data/teosintle_forecast_",
         time, "_", ssp, "_", model, ".tif")) %>%
  mutate(forecast = map2(model, ssp, ~geodata::cmip6_world(model = .x,
                                                            ssp = .y,
                                                            time = time,
                                                            var = "bioc",
                                                            res = 2.5,
                                                            path = paste0("data/teosintle_forecast_", time, "_", ssp, "_", .x, ".tif")))) %>%
  mutate(forecast = map2(forecast, path, ~{
    terra::writeRaster(.x, filename = .y)
    .x
  })) %>%
  select(-path) %>%
  write_rds("data/teosintle_forecasts.rds")




forecast_data <- geodata::cmip6_world(model = NULL,
                                      ssp = NULL,
                                      time = "2081-2100",
                                      var = "bioc",
                                      res = 2.5,
                                      path = "data")

#' Crop forecast data to the quadrant and save for future use                                    
forecast_data <- terra::crop(forecast_data, quadrant)
names(forecast_data) <- names(bioclim_data)

