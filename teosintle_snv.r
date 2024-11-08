library(tidyverse)
library(neuralnet)

#' This script fits a simple feed-forward neural network (FFN) to predict the
#' SNV frequency required to match forecasted bioclimatic data. 
#' These predictions are then used to calculate the genetic offsets
#' 
#' The FFN is trained with the neuralnet package. This package is **obsolete**.
#' It is recommended to use tensorflow through the keras package. However.
#' installing keras and tensorflow might be tricky.

#' Reading input data
Dat <- read_tsv("data/teosintle_snv_freqs.tsv")
Pops <- read_csv("data/teosintle_maxent_input.csv")
bioclim <- terra::rast("data/teosintle_bioclim_raster.tif")
bioclim_fut <- terra::rast("forecasts/teosintle_forecast_2061-2080_370_CMCC-ESM2.tif")
names(bioclim_fut) <- names(bioclim) # Make sure the names are the same!!


#' Add group information to the SNV data
Dat <- Dat %>%
  rename(lon = X, lat = Y) %>%
  left_join(Pops %>%
              filter(ID != "all"), by = c("lon", "lat")) %>%
  select(pop, lon, lat, ID, starts_with("refSNPS"), starts_with("candSNPs"))
Dat

#' Use terra to match bioclimatic data to the pops and create separate
#' tibble with the bioclimatic observation
pred_vars <- terra::extract(
    bioclim,
    terra::vect(Dat[, c("lon", "lat")],
      geom = c("lon", "lat"),
      crs = terra::crs(bioclim)
    )
  ) %>%
    select(-ID)

#' Same as above but for the forecasted data
fut_vars <- terra::vect(Dat[, c("lon", "lat")],
  geom = c("lon", "lat"),
  crs = terra::crs(bioclim_fut)
) %>%
  terra::extract(bioclim_fut, .) %>%
  select(-ID) 


#' Train neural network model to predict each SNV frequency
#' from the bioclimatic data
#' Do this for SNPs under selection
snps_to_model <- names(Dat)[str_detect(names(Dat), "^candSNPs")]
Res <- NULL
for (snp in snps_to_model) {
  # snp <- snps_to_model[1]
  cat(snp, "\n")

  # Match biclimatic data to current SNV
  d <- pred_vars %>%
    mutate(snv = Dat[[snp]])
  # d

  # Train neural network. Ideally, some form of cross-validation or
  # model evaluation should be performed. The keras package is recommended
  ffn1 <- neuralnet(snv ~ .,
    data = d,
    hidden = c(4, 2),
    linear.output = FALSE
  )

  # Make snv frequency predictions based on current and forecasted bioclimatic 
  # variables
  curr_est <- predict(ffn1, d)
  fut_pred <- predict(ffn1, fut_vars %>% mutate(snv = Dat[[snp]]))

  # Combine results
  r <- tibble(
    snp = snp, af = Dat[[snp]],
    curr_est = as.vector(curr_est),
    fut_pred = as.vector(fut_pred),
    pop = Dat$ID
  )

  # Clean to avoid issues, ideally this should be written in function form
  # to automate garbage collection
  rm(ffn1, curr_est, fut_pred)

  # Store results
  Res <- bind_rows(Res, r)
}

#' Calculate genetic offsets and plot them
Res %>%
  mutate(g_offset = abs(fut_pred - af)) %>%
  ggplot(aes(x = pop, y = g_offset)) +
  geom_boxplot(outlier.color = NA) +
  geom_point(position =  position_jitter(width = 0.2)) +
  theme_classic()


#' Should be done for all scenarios