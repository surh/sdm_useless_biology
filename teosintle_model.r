library(tidyverse)

#' Read data and combine bioclimatic variables with observations
Dat <- read_csv("data/teosintle_maxent_input.csv")
bioclim <- terra::rast("data/teosintle_bioclim_raster.tif")


#' Run maxent
#'
#' Performs variable selection and then runs maxent on presence only data
#'
#' @param Dat A data.frame or tibble with one row per presence
#' observation
#' @param bioclim A "SpatRaster" object. Should be already cropped to quadran
#' of interest
#' @param lon_lat A character vector with the names of the columns in Dat
#' corresponding to longitude and latitude. I that order. 
#'
#' @return A MaEnt model
#' @export
#'
#' @examples
model_maxent <- function(Dat, bioclim,
                         lon_lat = c("lon", "lat"),
                         cor_thres = 0.7) {
  # Remove colinear variables, use Variance Inflation Factor.
  selected_vars <- fuzzySim::corSelect(
    data = terra::extract(
      x = bioclim,
      y = Dat %>%
        select(lon, lat),
      ID = FALSE
    ),
    var.cols = names(bioclim),
    coeff = TRUE,
    cor.thresh = cor_thres,
    select = "VIF",
    method = "pearson"
  )
  bioclim <- bioclim[[selected_vars$selected.vars]]


  me <- dismo::maxent(
    x = as(bioclim, "Raster"),
    p = Dat %>%
      select(all_of(lon_lat)) %>%
      as.data.frame(),
    nbg = 1e4
  )

  return(me)
}

#' Run models on each group and combined
me_g1 <- model_maxent(Dat = Dat %>%
                        filter(ID == "g1"),
                      bioclim = bioclim,
                      lon_lat = c("lon", "lat"),
                      cor_thres = 0.7)
me_g2 <-  model_maxent(Dat = Dat %>%
                         filter(ID == "g2"),
                       bioclim = bioclim,
                       lon_lat = c("lon", "lat"),
                       cor_thres = 0.7)
me_all <-  model_maxent(Dat = Dat %>%
                          filter(ID == "all"),
                        bioclim = bioclim,
                        lon_lat = c("lon", "lat"),
                        cor_thres = 0.7)
op <- par(mfrow = c(1, 3))
dismo::plot(me_g1)
dismo::plot(me_g2)
dismo::plot(me_all)
par(op)

g1_pred <- terra::predict(as(bioclim, "Raster"), me_g1)
g2_pred <- terra::predict(as(bioclim, "Raster"), me_g2)

#' Forecast into the future. Instead of the CCSM4 2050 model under the rcp45
#' scenario we use the following
bioclim6180 <- terra::rast("data/teosintle_forecast_2061-2080_245_CMCC-ESM2.tif")
names(bioclim6180) <- names(bioclim)
bioclim81100 <- terra::rast("data/teosintle_forecast_2081-2100_245_CMCC-ESM2.tif")
names(bioclim81100) <- names(bioclim)

g1_forecast_6180 <- terra::predict(as(bioclim6180, "Raster"), me_g1)
g1_forecast_81100 <- terra::predict(as(bioclim81100, "Raster"), me_g1)


g2_forecast_6180 <- terra::predict(as(bioclim6180, "Raster"), me_g2)
g2_forecast_81100 <- terra::predict(as(bioclim81100, "Raster"), me_g2)

#' Read map and plot
quadrant_map <- terra::vect("data/teosintle_map/")


op <- par(mfrow=c(3,2))
terra::plot(quadrant_map, col = "white", main = "Now")
terra::plot(g1_pred,
  add = TRUE,
  col = colorRampPalette(c("white", "red"))(10)
)
points(Dat %>% filter(ID == "g1") %>%
  select(lon, lat), pch = 21, bg = "darkred", col = "black")
terra::plot(quadrant_map, add=TRUE, border='dark grey')

terra::plot(quadrant_map, col = "white", main = "now")
terra::plot(g2_pred,
  add = TRUE,
  col = colorRampPalette(c("white", "blue"))(10)
)
points(Dat %>% filter(ID == "g2") %>%
  select(lon, lat), pch = 21, bg = "darkblue", col = "black")
terra::plot(quadrant_map, add=TRUE, border='dark grey')



terra::plot(quadrant_map, col = "white", main = "2061-2080")
terra::plot(g1_forecast_6180,
  add = TRUE,
  col = colorRampPalette(c("white", "red"))(10)
)
points(Dat %>% filter(ID == "g1") %>%
  select(lon, lat), pch = 21, bg = "darkred", col = "black")
terra::plot(quadrant_map, add=TRUE, border='dark grey')

terra::plot(quadrant_map, col = "white", main = "2061-2081")
terra::plot(g2_forecast_6180,
  add = TRUE,
  col = colorRampPalette(c("white", "blue"))(10)
)
points(Dat %>% filter(ID == "g2") %>%
  select(lon, lat), pch = 21, bg = "darkblue", col = "black")
terra::plot(quadrant_map, add=TRUE, border='dark grey')



terra::plot(quadrant_map, col = "white", main = "2081-2100")
terra::plot(g1_forecast_81100,
  add = TRUE, 
  col = colorRampPalette(c("white", "red"))(10)
)
points(Dat %>% filter(ID == "g1") %>%
  select(lon, lat), pch = 21, bg = "darkred", col = "black")
terra::plot(quadrant_map, add=TRUE, border='dark grey')

terra::plot(quadrant_map, col = "white", main = "2081-2100")
terra::plot(g2_forecast_81100,
  add = TRUE,
  col = colorRampPalette(c("white", "blue"))(10)
)
points(Dat %>% filter(ID == "g2") %>%
  select(lon, lat), pch = 21, bg = "darkblue", col = "black")
terra::plot(quadrant_map, add=TRUE, border='dark grey')
par(op)



#' For more systematic evaluation of the models
#' Takes a while to run
# eval <- ENMeval::ENMevaluate(occs = Dat %>% filter(ID == "g1") %>% select(lon, lat) %>% as.data.frame(),
#                      envs = as(bioclim , "Raster"),
#                      tune.args = list(rm = c(0.5, 1 , 1.5),
#                                       fc = unlist(sapply(1:5, function(x) apply(combn(c("L","Q","H","P","T"), x), 2, function(y) paste(y, collapse = ""))))),
#                      bg.coords = NULL,
#                      clamp=FALSE,
#             algorithm = "maxnet",
#             method = 'jackknife',
#             parallel = TRUE,
#             numCores = 8)

#' Looad library
library(tidyverse)
library(biomod2)

#' Read observation data, current and forecasted bioclimatic data
Dat <- read_csv("data/teosintle_maxent_input.csv")
bioclim <- terra::rast("data/teosintle_bioclim_raster.tif")
bioclim6180 <- terra::rast("data/teosintle_forecast_2061-2080_245_CMCC-ESM2.tif")
names(bioclim6180) <- names(bioclim) # Make sure the names are the same!!
quadrant_map <- terra::vect("data/teosintle_map/")

#' # Using the biomod2 package

#' Prepare data
g1.bmd <- BIOMOD_FormatingData(
  resp.var = ifelse(Dat$ID == "g1", 1, NA),
  expl.var = bioclim,
  resp.xy = Dat[, c("lon", "lat")],
  resp.name = "g1", 
  PA.nb.rep = 4,
  PA.nb.absences = 1000,
  PA.strategy = "random",
  filter.raster = TRUE
)

#' Run models
g1.bm <- BIOMOD_Modeling(
  bm.format = g1.bmd,
  modeling.id = "g1.bm",
  models = c("MAXNET"),
  CV.strategy = "random",
  CV.nb.rep = 2,
  CV.perc = 0.8,
  metric.eval = c("ROC", "TSS"),
#  OPT.strategy = "tuned", # Normally you should tune the models
  var.import = 3
)

#' Project models
g1.bmp <- BIOMOD_Projection(
  bm.mod = g1.bm,
  proj.name = "Current",
  new.env = bioclim6180,
  models.chosen = "all",
  metric.binary = "all",
  metric.filter = "all",
  build.clamping.mask = TRUE
)

#' # Explore results with biomod2 functions

# # Get evaluation scores & variables importance
# get_evaluations(g1.bm)
# get_variables_importance(g1.bm) %>% 
#   as_tibble() %>%
#   filter(run != "allRun") %>%
#   ggplot(aes(x = expl.var, y = var.imp, col = run)) +
#   facet_wrap(~ expl.var, scales = "free_x") +
#   geom_boxplot(outlier.colour = NA) +
#   geom_point(position = position_jitterdodge()) +
#   theme_classic() +
#   theme(axis.text.x = element_blank())

# # Represent evaluation scores & variables importance
# get_calib_lines(g1.bm) %>% plot(g1.bmd, calib.lines = .)
# bm_PlotEvalMean(bm.out = g1.bm)
# bm_PlotEvalBoxplot(bm.out = g1.bm, group.by = c('algo', 'algo'))
# bm_PlotEvalBoxplot(bm.out = g1.bm, group.by = c('algo', 'run'))
# bm_PlotVarImpBoxplot(bm.out = g1.bm, group.by = c('expl.var', 'algo', 'algo'))
# bm_PlotVarImpBoxplot(bm.out = g1.bm, group.by = c('expl.var', 'algo', 'run'))
# bm_PlotVarImpBoxplot(bm.out = g1.bm, group.by = c('algo', 'expl.var', 'run'))

# # Represent response curves
# bm_PlotResponseCurves(bm.out = g1.bm, 
#                       models.chosen = get_built_models(g1.bm)[c(1:3)],
#                       fixed.var = 'median')

# # Plot projection
# plot(g1.bmp)




