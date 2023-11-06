
#' Now we run our own analysis
library(tidyverse)

#' Load the data
Dat <- read_tsv("data/carnegieae_gigantea_clean.tsv")
Dat


#' Create folds to validate predictions
ii_folds <- predicts::folds(x = Dat,
                             k = 5,
                             by = Dat$presence)
str(ii_folds)
table(ii_folds, Dat$presence)

# Let's call the first 4 folds *training* and the last fold *testing*
# folds

#' **Excercise: Build a GLM model with all bioclimatic variables as predictors,
#' and train it in the first 4 folds** 
m1 <- glm(presence ~ .,
          family = binomial(link = "logit"),
          data = Dat %>%
            filter(ii_folds != 5) %>%
            select(-gbifid, -latitude, -longitude))
summary(m1)

# **Excercise: Evaluate the model with testing data**
m1_eval <- predicts::pa_evaluate(p = Dat %>%
                                   filter(ii_folds == 5 & presence == 1),
                                 a = Dat %>%
                                   filter(ii_folds == 5 & presence == 0),
                                 model = m1,
                                 type = "response")
m1_eval@thresholds
predicts::plot(m1_eval, "ROC")



#' We compare prediction probs for presences and pseudo absences 
Dat %>%
  filter(ii_folds == 5) %>%
  bind_cols(prediction = predict(m1,
                                 Dat %>%
                                   filter(ii_folds == 5),
                                 type = "response")) %>%
  ggplot(aes(x = factor(presence), y = prediction)) +
  geom_violin() +
  geom_point(position = position_jitter(width = 0.2)) +
  theme_classic()
  

#' Plot in the map
bioclim_data <- terra::rast("data/carnegieae_bioclim_raster.tif")
quadrant_map <- terra::vect("data/carnegieae_map/")
map_predictions <- terra::predict(bioclim_data, m1, type = "response")
terra::plot(map_predictions)

#' Plot only regions with higher confidence
terra::plot(quadrant_map, 
            axes = TRUE, 
            col = "grey95")
terra::plot(map_predictions > m1_eval@thresholds$max_spec_sens, 
            add = TRUE, 
            legend = FALSE, 
            col = c(NA, "olivedrab"))
points(x = Dat$longitude[ Dat$presence == 1], 
       y = Dat$latitude[ Dat$presence == 1], 
       col = "red",
       pch = 20, 
       cex = 0.75)

#' **Excercise: Forecast species range in the future according to model**
forecast_bioclim <- terra::rast("data/carnegieae_forecast_raster.tif")
map_forecast <- terra::predict(forecast_bioclim, m1, type = "response")

terra::plot(quadrant_map, 
            axes = TRUE, 
            col = "grey95")
terra::plot(map_forecast > m1_eval@thresholds$max_spec_sens, 
            add = TRUE, 
            legend = FALSE, 
            col = c(NA, "olivedrab"))
points(x = Dat$longitude[ Dat$presence == 1], 
       y = Dat$latitude[ Dat$presence == 1], 
       col = "red",
       pch = 20, 
       cex = 0.75)

#' #' # Make a better model
#' Dat %>%
#'   select(-gbifid, -latitude, -longitude, -presence) %>%
#'   cor %>% replace(list = . < 0.7, values = 0) %>%
#'   heatmap()
#' 
#' 
#' m2 <- glm(presence ~ wc2.1_2.5m_bio_12 +
#'             I(wc2.1_2.5m_bio_12^2) +
#'             wc2.1_2.5m_bio_18 +
#'             I(wc2.1_2.5m_bio_18^2),
#' 
#'           family = binomial(link = "logit"),
#'           data = Dat %>%
#'             filter(ii_folds != 5) %>%
#'             select(-gbifid, -latitude, -longitude))
#' summary(m2)
#' AIC(m1,m2)
#' m2_eval <- predicts::pa_evaluate(p = Dat %>%
#'                                    filter(ii_folds == 5 & presence == 1),
#'                                  a = Dat %>%
#'                                    filter(ii_folds == 5 & presence == 0),
#'                                  model = m2,
#'                                  type = "response")
#' m2_eval

#' # Random forest
rf1 <- randomForest::randomForest(x = Dat %>%
                                    select(-gbifid, -latitude, -longitude, -presence),
                                  y = Dat$presence %>% as.factor,
                                  # y = Dat$presence, 
                                  ntree = 1e3,
                                  nodesize = 10,
                                  importance = TRUE)
rf1
randomForest::importance(rf1)
randomForest::varImpPlot(rf1)
randomForest::partialPlot(x = rf1,
                          pred.data = Dat %>%
                            select(-gbifid, -latitude, -longitude, -presence) %>%
                            as.data.frame(), 
                          x.var = wc2.1_2.5m_bio_12,
                          which.class = "1")

#' # Maxent
me1 <- dismo::maxent(x = terra::extract(x = bioclim_data,
                                        y = Dat %>%
                                          select(longitude, latitude),
                                        ID = FALSE),
                     p = Dat$presence)
dismo::response(me1)
me1_predictions <- dismo::predict(me1, forecast_bioclim)

terra::plot(quadrant_map, 
            axes = TRUE, 
            col = "grey95")
terra::plot(me1_predictions, 
            add = TRUE, 
            legend = TRUE)
points(x = Dat$longitude[ Dat$presence == 1], 
       y = Dat$latitude[ Dat$presence == 1], 
       col = "red",
       pch = 20, 
       cex = 0.75)

