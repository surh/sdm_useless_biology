library(tidyverse)
library(neuralnet)

Dat <- read_tsv("data/teosintle_snv_freqs.tsv")
Pops <- read_csv("data/teosintle_maxent_input.csv")
bioclim <- terra::rast("data/teosintle_bioclim_raster.tif")
bioclim_fut <- terra::rast("forecasts/teosintle_forecast_2061-2080_370_CMCC-ESM2.tif")
names(bioclim_fut) <- names(bioclim) # Make sure the names are the same!!


Dat <- Dat %>%
  rename(lon = X, lat = Y) %>%
  left_join(Pops %>%
              filter(ID != "all"), by = c("lon", "lat")) %>%
  select(pop, lon, lat, ID, starts_with("refSNPS"), starts_with("candSNPs"))
Dat

Dat


pred_vars <- terra::extract(
    bioclim,
    terra::vect(Dat[, c("lon", "lat")],
      geom = c("lon", "lat"),
      crs = terra::crs(bioclim)
    )
  ) %>%
    select(-ID)

fut_vars <- terra::vect(Dat[, c("lon", "lat")],
  geom = c("lon", "lat"),
  crs = terra::crs(bioclim_fut)
) %>%
  terra::extract(bioclim_fut, .) %>%
  select(-ID) 



snps_to_model <- names(Dat)[str_detect(names(Dat), "^candSNPs")]
Res <- NULL
for (snp in snps_to_model) {
  # snp <- snps_to_model[1]
  cat(snp, "\n")

  d <- pred_vars %>%
    mutate(snv = Dat[[snp]])
  d

  ffn1 <- neuralnet(snv ~ .,
    data = d,
    hidden = c(4, 2),
    linear.output = FALSE
  )

  curr_est <- predict(ffn1, d)
  fut_pred <- predict(ffn1, fut_vars %>% mutate(snv = Dat[[snp]]))


  r <- tibble(
    snp = snp, af = Dat[[snp]],
    curr_est = as.vector(curr_est),
    fut_pred = as.vector(fut_pred),
    pop = Dat$ID
  )

  rm(ffn1, curr_est, fut_pred)

  Res <- bind_rows(Res, r)
}

Res %>%
  mutate(g_offset = abs(fut_pred - af)) %>%
  ggplot(aes(x = pop, y = g_offset)) +
  geom_boxplot(outlier.color = NA) +
  geom_point(position =  position_jitter(width = 0.2)) +
  theme_classic()

