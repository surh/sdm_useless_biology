#' Looad libraries
library(tidyverse)
library(biomod2)

#' Model scenario and population
#' 
#' This function models the distribution of a population under current and future
#' bioclimatic conditions
#' 
#' @param Dat Data frame with the following columns: ID, lon, lat
#' @param bioclim SpatRaster with bioclimatic variables
#' @param bioclim_fut SpatRaster with forecasted bioclimatic variables under
#' a specific scenario
#' @param pop ID of population to model. Should be in column ID
#' @return List with the results of the model
#' @export
model_scenario <- function(Dat, bioclim, bioclim_fut, pop){ 
    
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

    return(list(
        bm1 = bm1, bmd = bmd,
        bmp_curr = bmp1_curr, bmp_fut = bmp1_fut
    ))
}


#' Calculate area change
#'
#' This function calculates forecasted area change (in pixels) under a
#' specific future scenario
#'
#' @param res List with the results of the model produced by model_scenario()
#' @return Data frame with the number of pixels that are forecasted to be
#' lost, gained, or remain stable
calc_change <- function(res) {
    BIOMOD_RangeSize(
        proj.current = get_predictions(res$bmp_curr,
            metric.binary = "TSS",
            model.as.col = TRUE
        ),
        proj.future = get_predictions(res$bmp_fut,
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
        mutate(type = factor(type, levels = c("Loss", "Stable1", "Gain")))
}

#' Read input data for current pops
Dat <- read_csv("data/teosintle_maxent_input.csv")
bioclim <- terra::rast("data/teosintle_bioclim_raster.tif")
quadrant_map <- terra::vect("data/teosintle_map/")

#' Lit all forecast files, and select those for the 2061-2080 period
forecast_files <- list.files("forecasts/") 
forecast_files <- forecast_files[str_detect(forecast_files, "_2061-2080")]

#' Loop over forecast files, fit the maxent sdm model, and calculate the area
#' change for each population and scenario.
Res <- NULL
# for (forecast_file in forecast_files[sample(1:length(forecast_files), 3)]) {
for (forecast_file in forecast_files) {
  # forecast_file <- forecast_files[1]
  bioclim_fut <- terra::rast(file.path("forecasts", forecast_file))
  names(bioclim_fut) <- names(bioclim) # Make sure the names are the same!!

  # Model each population
  res_g1 <- model_scenario(Dat, bioclim, bioclim_fut, "g1")
  res_g2 <- model_scenario(Dat, bioclim, bioclim_fut, "g2")

  # Calculate area change for each population
  loss_g1 <- calc_change(res_g1)
  loss_g2 <- calc_change(res_g2)

  # Combine results
  res <- loss_g1 %>%
      filter(type == "Loss") %>%
      mutate(pop = "g1") %>%
      bind_rows(
          loss_g2 %>%
              filter(type == "Loss") %>%
              mutate(pop = "g2")
      ) %>%
      mutate(model = forecast_file)
  
  # Store results
  Res <- bind_rows(Res, res)
}

# Plot results
Res
p1 <- Res %>%
    separate(model,
        into = c("species", "forecast", "years", "ssp", "model"),
        sep = "_"
    ) %>%
    select(-species, -forecast) %>%
    mutate(model = str_remove(model, ".tif$")) %>%
    ggplot(aes(x = pop, y = n_pixels)) +
    facet_wrap(~ssp) +
    geom_boxplot(outlier.color = NA) +
    geom_point(aes(color = model),
        position = position_jitter(width = 0.2),
        size = 0.5
    ) +
    geom_hline(yintercept = 500) +
    labs(
        title = "Area loss proyected in 2061-2080",
        y = "Number of pixels lost", x = "Population"
    ) +
    theme_classic()
p1
ggsave("teosintle_area_loss.png", p1, width = 10, height = 5)





    












