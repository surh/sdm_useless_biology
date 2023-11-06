# setwd("/Users/sur/lab/exp/2023/today3/")

#' Based on: https://jcoliver.github.io/learn-r/011-species-distribution-models.html
# library(terra)
# library(geodata)
# library(predicts)
library(tidyverse)


#' First we download bioclimate data, and world map. This are ignored by git
#' in my repo
bioclim_data <- geodata::worldclim_global(var = "bio",
                                          res = 2.5,
                                          path = "data/")
world_map <- geodata::world(resolution = 3,
                            path = "data/")

#' Read data downloaded from: https://raw.githubusercontent.com/jcoliver/learn-r/gh-pages/data/Carnegiea-gigantea-GBIF.csv
Obs <- read_csv("data/Carnegiea-gigantea-GBIF.csv") %>%
    filter(!is.na(latitude))
summary(Obs)

#' Find quadrant where species is located and plot, and plot species data 
#' with climate. The index in the plot call indicates which variable to
#' plot. For variable definitions see: https://www.worldclim.org/data/bioclim.html
max_lat <- ceiling(max(Obs$latitude))
min_lat <- floor(min(Obs$latitude))
max_lon <- ceiling(max(Obs$longitude))
min_lon <- floor(min(Obs$longitude))
quadrant <- terra::ext(x = c(min_lon, max_lon, min_lat, max_lat))
quadrant <- quadrant * 1.25

quadrant_map <- terra::crop(x = world_map, y = quadrant)
bioclim_data <- terra::crop(x = bioclim_data, y = quadrant)

terra::plot(bioclim_data[[1]])
# terra::plot(quadrant_map,
#             axes = TRUE, 
#             col = "grey95")
points(x = Obs$longitude, 
       y = Obs$latitude, 
       col = "red", 
       pch = 20, 
       cex = 0.75)









# Set the seed for the random-number generator to ensure results are similar
set.seed(20231106)

# Randomly sample points (same number as our observed points)
background <- terra::spatSample(x = bioclim_data,
                                size = 1000,    # generate 1,000 pseudo-absence points
                                values = FALSE, # don't need values
                                na.rm = TRUE,   # don't sample from ocean
                                xy = TRUE)      # just need coordinates

# Look at first few rows of background
head(background)

# Plot the base map
plot(my_map,
     axes = TRUE, 
     col = "grey95")

# Add the background points
points(background,
       col = "grey30",
       pch = 1,
       cex = 0.75)

# Add the points for individual observations
points(x = Obs$longitude, 
       y = Obs$latitude, 
       col = "olivedrab", 
       pch = 20, 
       cex = 0.75)





Dat <- Obs %>%
  mutate(presence = 1) %>%
  select(latitude, longitude, presence) %>%
  bind_rows(background %>%
              as_tibble() %>%
              mutate(presence = 0) %>%
              rename(longitude = x, latitude = y))
Dat


Dat <- Dat %>% 
  bind_cols(terra::extract(x = bioclim_data,
                       y = Dat %>% select(longitude, latitude),
                       ID = FALSE))
Dat


Dat_folds <- predicts::folds(x = Dat,
                             k = 5,
                             by = Dat$presence)
table(Dat_folds)
str(Dat_folds)





