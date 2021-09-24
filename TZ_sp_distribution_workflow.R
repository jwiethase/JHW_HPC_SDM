args = commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages(library(INLA, quietly=TRUE))

library(data.table)
library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)

source("/users/jhw538/scratch/TZ_INLA/source/JHW_FitModel.R")
source("/users/jhw538/scratch/TZ_INLA/source/MakeSpatialRegion.R")
source("/users/jhw538/scratch/TZ_INLA/source/MakeIntegrationStack.R")
source("/users/jhw538/scratch/TZ_INLA/source/GetNearestCovariate.R")
source("/users/jhw538/scratch/TZ_INLA/source/JHW_MakeProjectionGrid.R")
source("/users/jhw538/scratch/TZ_INLA/source/JHW_MakeBinomStack.R")
source("/users/jhw538/scratch/TZ_INLA/source/MakePointsStack.R")
source("/users/jhw538/scratch/TZ_INLA/source/misc_functions.R")

# The species for the distribution model
species_list = c('Cisticola juncidis')
#species_list = c('Passer domesticus', 'Cisticola juncidis', 'Estrilda astrild', 'Histurgops ruficauda', 'Ploceus nigricollis', 
#                 'Cisticola brunnescens', 'Chrysococcyx cupreus', 'Tauraco hartlaubi', 'Ploceus castaneiceps', 'Nigrita canicapilla', 
#                 'Nectarinia kilimensis', 'Lanius collaris', 'Terpsiphone viridis', 'Oriolus auratus', 'Bubo capensis', 'Bubo africanus')
i <- as.numeric(args[1])

time_list = c('20s')

# Max.edge based on estimated range from pilot model
estimated_range = 2
max.edge = estimated_range/8

load(paste0("/users/jhw538/scratch/TZ_INLA/data/TZ_INLA_model_file_E", round(max.edge, digits = 3), ".RData"))

# The prior ranges to use for the model
# Prior range could be 10*max.edge. Should not be smaller than mesh resolution.
range0_ls = c(5) 
Prange_ls = c(0.5)
range_combs <- crossing(range0_ls, Prange_ls)

sigma0_ls = c(3)
Psigma_ls =  c(0.01)
sigma_combs <- crossing(sigma0_ls, Psigma_ls)

all_combs <- crossing(range_combs, sigma_combs)

# Process data -------------------------------------------------------
#for (i in 1:length(species_list)){
  species = species_list[i]
  print(species)
  for (j in 1:length(time_list)){
   time = time_list[j]

if(time == '60s'){
   ebird_full <- ebird_full %>% 
      filter(date > '1960-01-01', date < '2000-01-01')
}

if(time == '20s'){
   ebird_full <- ebird_full %>% 
      filter(date > '2000-01-01')
}

# Filter eBird data
ebird_filtered <- ebird_full %>% 
      filter(APPROVED == 1,  # Only keep reviewed and approved records
             `ALL SPECIES REPORTED` == 1,           # Only keep complete checklists
             `EFFORT DISTANCE KM` < 15,
             `DURATION MINUTES` >= 5,
             `DURATION MINUTES` <= 240)   
print(summary(ebird_filtered))
ebird_sp <- make_ebird_sp(species, TZ_outline)

atlas_filtered <- atlas_full %>% 
   mutate(Scientific = trimws(Scientific, which = 'both')) %>% 
   filter(time_period == time, Scientific == species) %>% 
   mutate(presence = ifelse(occurrence == 1, TRUE, FALSE)) %>% 
   dplyr::select(-V1); if(is_empty(atlas_filtered$presence)){print("ERROR: No Atlas data available")}

atlas_sp <- SpatialPointsDataFrame(
   coords = atlas_filtered[, c("Long", "Lat")],
   data = data.frame(presence = atlas_filtered$presence, effort = atlas_filtered$effort),
   proj4string = crs(proj)
)

# Standardize effort variables
range01 <- function(x){(x - min(x))/(max(x) - min(x))}
ebird_sp$duration_minutes <- range01(ebird_sp$duration_minutes)
atlas_sp$effort <- range01(atlas_sp$effort)

# Create data stacks for each of the different data types
stk.eBird <- MakeBinomStack(
   observs = ebird_sp,
   data = variables_no_BG,
   mesh = Mesh$mesh,
   presname = "presence",
   tag = "eBird",
   InclCoords = TRUE,
   add_effort = TRUE, 
   effort_names = names(ebird_sp@data[ , -which(names(ebird_sp) %in% c("presence"))]))

stk.atlas <- MakeBinomStack(
   observs = atlas_sp,
   data = variables_no_BG,
   mesh = Mesh$mesh,
   presname = "presence",
   tag = "atlas",
   InclCoords = TRUE,
   add_effort = TRUE, 
   effort_names = "effort"
)

C.F. <- list(
   mean = list(int.eBird = 0, int.atlas = 0),
   mean.intercept = 0,
   prec = list(int.eBird = 1, int.atlas = 1),
   prec.intercept = 0.001
)

# pcmatern stands for Penalised Complexity priors for parameters of the Matern field. 

for (i in 1:NROW(all_combs)){
   
   comb_values <- all_combs[i, ]
   spde <- inla.spde2.pcmatern(
         mesh = Mesh$mesh,
         alpha = 2,
         prior.range = c(comb_values$range0_ls, comb_values$Prange_ls),   # e.g. 5, 0.01 = 1% probability that range is smaller than 5 degrees
         prior.sigma = c(comb_values$sigma0_ls, comb_values$Psigma_ls)
      )

   
   form <- formula(
      resp ~ 0 + z.rain1.s + z.rain2.s +    # 0 suppresses default intercept
         z.temp_min1.s + z.temp_min2.s + z.temp_max1.s + z.temp_max2.s + z.pop1.s + z.pop2.s + z.canopy1.s + z.canopy2.s + z.BG1.s + z.BG2.s + indicator +
         Intercept + x + y + int.eBird + duration_minutes + 
         int.atlas + effort +
         f(i, model = spde)
      )

   
   model <- FitModel(
      stk.ip, stk.eBird, stk.atlas, stk.pred$stk,
      formula = form,
      CovNames = NULL,
      mesh = Mesh$mesh,
      predictions = TRUE,
      control.fixed = C.F.,
      waic = FALSE,
      cpo = FALSE, 
      dic = FALSE,
      nthreads = 16,
      #lincomb = all_lc,
      verbose = TRUE
      )
   
   
   filename = paste0("/users/jhw538/scratch/TZ_INLA/INLA_TZ_toMove/Model_", gsub(" ", "_", species, fixed = TRUE), "_", time, 
                     "_Output_r", comb_values$range0_ls, "_", comb_values$Prange_ls, "_s", comb_values$sigma0_ls, "_", comb_values$Psigma_ls, "_E", 
                     gsub("\\.", "_", as.character(round(max.edge, digits = 3))),  ".RData")

   print(filename)
   
   save(
        C.F., spde, form, model, comb_values, stk.pred, Mesh,
         file = filename
         )

}
}

#}

