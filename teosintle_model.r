#' Looad libraries
library(tidyverse)
library(biomod2)

#' Read observation data, current and forecasted bioclimatic data
#' Fore the forecast model. Instead of the CCSM4 2050 model under the rcp45
#' scenario we use the CMCC-ESM2 model under de SSP370 scenario
Dat <- read_csv("data/teosintle_maxent_input.csv")
bioclim <- terra::rast("data/teosintle_bioclim_raster.tif")
bioclim_fut <- terra::rast("data/teosintle_forecast_2061-2080_370_CMCC-ESM2.tif")
names(bioclim_fut) <- names(bioclim) # Make sure the names are the same!!
quadrant_map <- terra::vect("data/teosintle_map/")

#' # Using the biomod2 package

pop <- "g2"

#' Prepare data
g1.bmd <- BIOMOD_FormatingData(
    resp.var = rep(1, sum(Dat$ID == pop)),
    expl.var = bioclim,
    resp.xy = Dat %>% filter(ID == pop) %>% select(lon, lat),
    resp.name = pop,
    PA.nb.rep = 4,
    PA.nb.absences = 1000,
    PA.strategy = "random",
    filter.raster = TRUE
  )

#' Run models
g1.bm <- BIOMOD_Modeling(
  bm.format = g1.bmd,
  modeling.id = paste0(pop, ".bm"),
  models = c("MAXNET"),
  CV.strategy = "random",
  CV.nb.rep = 2,
  CV.perc = 0.8,
  metric.eval = c("ROC", "TSS"),
  # OPT.strategy = "tuned", # Normally you should tune the models
  var.import = 3
)

#' Project models
g1_curr.bmp <- BIOMOD_Projection(
  bm.mod = g1.bm,
  proj.name = "Now",
  new.env = bioclim,
  models.chosen = "all",
  metric.binary = "all",
  metric.filter = "all",
  build.clamping.mask = TRUE
)

g1_fut.bmp <- BIOMOD_Projection(
  bm.mod = g1.bm,
  proj.name = "Future",
  new.env = bioclim_fut,
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
# plot(g1_fut.bmp)

#' # Plot results with terra
op <- par(mfrow=c(1, 2))

terra::plot(quadrant_map, col = "white", main = "Now")
terra::plot(terra::unwrap(g1_curr.bmp@proj.out@val)[[3]],
  add = TRUE,
  col = colorRampPalette(c("white", "red"))(10)
)
points(Dat %>% filter(ID == pop) %>%
  select(lon, lat), pch = 21, bg = "darkred", col = "black")
terra::plot(quadrant_map, add = TRUE, border = "dark grey")

terra::plot(quadrant_map, col = "white", main = "Future")
terra::plot(terra::unwrap(g1_fut.bmp@proj.out@val)[[3]],
  add = TRUE,
  col = colorRampPalette(c("white", "red"))(10)
)
points(Dat %>% filter(ID == pop) %>%
  select(lon, lat), pch = 21, bg = "darkred", col = "black")
terra::plot(quadrant_map, add = TRUE, border = "dark grey")


op <- par(op)




