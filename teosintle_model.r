#' Looad libraries
library(tidyverse)
library(biomod2)

#' Read observation data, current and forecasted bioclimatic data
#' Fore the forecast model. Instead of the CCSM4 2050 model under the rcp45
#' scenario we use the CMCC-ESM2 model under de SSP370 scenario
Dat <- read_csv("data/teosintle_maxent_input.csv")
bioclim <- terra::rast("data/teosintle_bioclim_raster.tif")
bioclim_fut <- terra::rast("forecasts/teosintle_forecast_2061-2080_370_CMCC-ESM2.tif")
names(bioclim_fut) <- names(bioclim) # Make sure the names are the same!!
quadrant_map <- terra::vect("data/teosintle_map/")

#' # Using the biomod2 package

pop <- "g1"

#' Prepare data
bmd <- BIOMOD_FormatingData(
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
bm1 <- BIOMOD_Modeling(
  bm.format = bmd,
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
bmp1_curr <- BIOMOD_Projection(
  bm.mod = bm1,
  proj.name = "Now",
  new.env = bioclim,
  models.chosen = "all",
  metric.binary = "all",
  metric.filter = "all",
  build.clamping.mask = TRUE
)

bmp1_fut <- BIOMOD_Projection(
  bm.mod = bm1,
  proj.name = "Future",
  new.env = bioclim_fut,
  models.chosen = "all",
  metric.binary = "all",
  metric.filter = "all",
  build.clamping.mask = TRUE
)

# #' # Explore results with biomod2 functions
# # Represent evaluation scores & variables importance
# get_calib_lines(bm1) %>% plot(bmd, calib.lines = .)
# bm_PlotEvalMean(bm.out = bm1)
# bm_PlotEvalBoxplot(bm.out = bm1, group.by = c('algo', 'algo'))
# bm_PlotEvalBoxplot(bm.out = bm1, group.by = c('algo', 'run'))
# bm_PlotVarImpBoxplot(bm.out = bm1, group.by = c('expl.var', 'algo', 'algo'))
# bm_PlotVarImpBoxplot(bm.out = bm1, group.by = c('expl.var', 'algo', 'run'))
# bm_PlotVarImpBoxplot(bm.out = bm1, group.by = c('algo', 'expl.var', 'run'))

# # Represent response curves
# bm_PlotResponseCurves(bm.out = bm1,
#                       models.chosen = get_built_models(bm1)[c(1:3)],
#                       fixed.var = 'median')

# # Plot projection
# plot(g1_fut.bmp)

# # Plot range size
# bm_PlotRangeSize(BIOMOD_RangeSize(
#   proj.current = get_predictions(bmp1_curr,
#     metric.binary = "TSS",
#     model.as.col = TRUE
#   ),
#   proj.future = get_predictions(bmp1_fut,
#     metric.binary = "TSS",
#     model.as.col = TRUE
#   )
# ))

#' Plot evaluation scores
get_evaluations(bm1) %>%
  filter(run != "allRun") %>%
  ggplot(aes(x = run, y = validation)) +
  facet_wrap(~metric.eval) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitter(), size = 3) +
  ggtitle(label = paste0("Population: ", pop)) +
  theme_classic()

#' Plot variables importance
get_variables_importance(bm1) %>%
  filter(run != "allRun") %>%
  ggplot(aes(x = expl.var, y = var.imp, col = run)) +
  facet_wrap(~expl.var, scales = "free_x") +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge()) +
  ggtitle(label = paste0("Population: ", pop)) +
  theme_classic() +
  theme(axis.text.x = element_blank())

#' # Plot results with terra
op <- par(mfrow=c(1, 2))

terra::plot(quadrant_map, col = "white", main = "Now")
terra::plot(terra::unwrap(bmp1_curr@proj.out@val)[[3]],
  add = TRUE,
  col = colorRampPalette(c("white", "red"))(10)
)
points(Dat %>% filter(ID == pop) %>%
  select(lon, lat), pch = 21, bg = "darkred", col = "black")
terra::plot(quadrant_map, add = TRUE, border = "dark grey")

terra::plot(quadrant_map, col = "white", main = "Future")
terra::plot(terra::unwrap(bmp1_fut@proj.out@val)[[3]],
  add = TRUE,
  col = colorRampPalette(c("white", "red"))(10)
)
points(Dat %>% filter(ID == pop) %>%
  select(lon, lat), pch = 21, bg = "darkred", col = "black")
terra::plot(quadrant_map, add = TRUE, border = "dark grey")

op <- par(op)

#' # Plot response curves
bm_PlotResponseCurves(bm1,
  models.chosen = get_built_models(bm1)[1:3],
  fixed.var = "median",
  do.plot = FALSE
)$tab %>%
  ggplot(aes(x = expl.val, y = pred.val, col = pred.name)) +
  facet_wrap(~expl.name, scales = "free_x") +
  geom_line() +
  ggtitle(label = paste0("Population: ", pop)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

#' # Plot range size

BIOMOD_RangeSize(
  proj.current = get_predictions(bmp1_curr,
    metric.binary = "TSS",
    model.as.col = TRUE
  ),
  proj.future = get_predictions(bmp1_fut,
    metric.binary = "TSS",
    model.as.col = TRUE
  )
)$Compt.By.Models %>%
  as.data.frame() %>%
  rownames_to_column(var = "meta") %>%
  as_tibble() %>%
  separate(meta,
    into = c("spec", "pa_group", "run", "model"),
    sep = "_"
  ) %>%
  filter(run != "allRun") %>%
  mutate(id = paste0(pa_group, "_", run)) %>%
  select(id, run, Stable1, Gain, Loss) %>%
  pivot_longer(
    cols = -c("run", "id"),
    names_to = "type",
    values_to = "n_pixels"
  ) %>%
  mutate(type = factor(type2080, levels = c("Loss", "Stable1", "Gain"))) %>%
  ggplot(aes(x = run, y = n_pixels, fill = type)) +
  geom_bar(stat = "identity", position = "fill", col = NA) +
  ggtitle(label = paste0("Population: ", pop)) +
  ylab("Proportion of pixels") +
  theme_classic()


# Y la población g2?
# Qué tan robustos son los resultado ante lo diferentes escenarios
# climático futuros?
