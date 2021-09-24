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


# Max.edge based on estimated range from pilot model
estimated_range = 2
max.edge = estimated_range/8

# Import data --------------------------------------------------------
# Bird data
ebird_full <- fread("/users/jhw538/scratch/TZ_INLA/data/ebd_TZ_relMay-2021/ebd_TZ_relMay-2021.txt") %>% 
      mutate(date = ymd(get("OBSERVATION DATE"))) %>%
      filter(!is.na(`DURATION MINUTES`), !is.na(`EFFORT DISTANCE KM`))

atlas_full <- fread("/users/jhw538/scratch/TZ_INLA/data/TZ_bird_atlas_data.csv") %>%
      filter(!is.na(effort))

# Map files and covariates
proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
TZ_outline <- readOGR('/users/jhw538/scratch/TZ_INLA/data/TZ_simpler.shp')

# Check how mesh will use the border. Shows potential problems in shapefile
# region.bdry <- inla.sp2segment(TZ_outline); plot(region.bdry) 

TZ_canopy_height <- raster('/users/jhw538/scratch/TZ_INLA/data/TZ_GEDI_500m.tif') %>% mask(., TZ_outline) 
TZ_ann_rain_2000s <- raster('/users/jhw538/scratch/TZ_INLA/data/TZ_annual_median_rain_00_20.tif') %>% mask(., TZ_outline) %>% projectRaster(., TZ_canopy_height)
TZ_min_temp_2000s <- raster('data/TZ_MODIS_coldest_temperature_00_20.tif') %>% mask(., TZ_outline) %>% projectRaster(., TZ_canopy_height)
TZ_max_temp_2000s <- raster('data/TZ_MODIS_hottest_temperature_00_20.tif') %>% mask(., TZ_outline) %>% projectRaster(., TZ_canopy_height)TZ_population <- raster('/users/jhw538/scratch/TZ_INLA/data/TZ_worldpop_2020_500m.tif') %>% mask(., TZ_outline) %>% projectRaster(., TZ_canopy_height)
degind_2010_14 <- raster('/users/jhw538/scratch/TZ_INLA/data/BG_2010_14_500m.tif') %>% mask(., TZ_outline) %>% projectRaster(., TZ_canopy_height)

# Prepare BG layer ---------------------------------------------------
# BG layer has large gaps in data, this needs to be accounted for. Manually create indicator layer and BG interaction layer,
# to make sure model ignores areas with NA BG

# Make sure all NaN are NA
degind_2010_14[is.nan(degind_2010_14)] <- NA

# Make indicator layer where 0 is NA in BG, and 1 is value in BG
indicator <- degind_2010_14 
values(indicator)[!is.na(values(indicator))] <- 1
values(indicator)[is.na(values(indicator))] <- 0

# Stack and convert to spdf
variables <- stack(degind_2010_14, indicator, TZ_ann_rain_2000s, TZ_min_temp_2000s, TZ_max_temp_2000s, TZ_population, TZ_canopy_height)
variables <- as(variables, "SpatialPointsDataFrame")
names(variables) <- c("BG", "indicator", "rain", "temp_min", "temp_max", "pop",  "canopy")

# BG_spdf <-  as(combined, "SpatialPointsDataFrame")
# names(BG_spdf) <- c()

# Make spdf without NAs, do transformation and GAM prep
variables_no_NA <- variables[!is.na(variables$BG), ] 
variables_no_NA <- prepare_GAM(variables_no_NA, "BG")

# Put NAs back, replace with any other value, here 0 (will be ignored)
variables$z.BG1.s <- variables$BG
variables$z.BG1.s[!is.na(variables$z.BG1.s)] <- variables_no_NA$z.BG1.s
variables$z.BG1.s[is.na(variables$z.BG1.s)] <- 0

variables$z.BG2.s <- variables$BG
variables$z.BG2.s[!is.na(variables$z.BG2.s)] <- variables_no_NA$z.BG2.s
variables$z.BG2.s[is.na(variables$z.BG2.s)] <- 0

# Remove BG layer, no longer needed
variables_no_BG <- variables[ , !names(variables) %in% "BG"]

variables_no_BG <- variables_no_BG[!is.na(rowSums(variables_no_BG@data)), ]
# 
# # Make map of all covariates
# coords <- coordinates(variables)
# 
# #adds ID column for each cell
# variables$ID <- 1:NROW(coords)

variables_no_BG <- prepare_GAM(variables_no_BG, "rain")
variables_no_BG <- prepare_GAM(variables_no_BG, "temp_min")
variables_no_BG <- prepare_GAM(variables_no_BG, "temp_max")
variables_no_BG <- prepare_GAM(variables_no_BG, "pop")
variables_no_BG <- prepare_GAM(variables_no_BG, "canopy")

## Make some linear combinations:
lincombs_wrapper <- function(file1, file2, new_var){
   lc <- inla.make.lincombs(int.eBird = rep(1, 100),  # Change for all
                            int.atlas = rep(1, 100),
                            file1 = file1,
                            file2 = file2)
   names(lc) <- paste0(new_var, 1:100)
   assign(new_var, lc, envir = .GlobalEnv)
}
# lincombs_wrapper(z.temp1.s, z.temp2.s, "z.temp_lc")
lincombs_wrapper(z.rain1.s, z.rain2.s, "z.rain_lc")
lincombs_wrapper(z.pop1.s, z.pop2.s, "z.pop_lc")
lincombs_wrapper(z.canopy1.s, z.canopy2.s, "z.canopy_lc")
all_lc <- c(z.temp_lc, z.rain_lc, z.pop_lc, z.canopy_lc)

# Create the mesh to approximate the area and the spatial field
# Meshpars <- list(max.edge = c(0.05, 0.4), offset = c(0.1, 0.4), cutoff = 0.1)    # Fine mesh
Meshpars <- list(max.edge = c(max.edge, max.edge*4), 
                 offset = c(max.edge, max.edge*5), 
                 cutoff = max.edge/2)

Mesh <- MakeSpatialRegion(
   data = NULL,
   bdry = TZ_outline,
   meshpars = Meshpars,
   proj = proj
)

pdf(paste0("/users/jhw538/scratch/TZ_INLA/out_figures/TZ_mesh_E", max.edge, ".pdf"))
plot(Mesh$mesh)
dev.off()
# Make stack for background mesh, the prediction stack with NA as response.
# Projection stack = Integration stack
stk.ip <- MakeIntegrationStack(
   mesh = Mesh$mesh,
   data = variables_no_BG,
   area = Mesh$w,
   tag = "ip",
   InclCoords = TRUE
)

# Make data for projections

# This sets the resolution of the predictions. For the final version we use 
# Nxy.scale <- 0.01, but have changed it here so the code will run more quickly.
if(!exists("Nxy.scale")) Nxy.scale <- 0.1  # about 10km resolution

Boundary <- Mesh$mesh$loc[Mesh$mesh$segm$int$idx[, 2], ]
Nxy.size <- c(diff(range(Boundary[, 1])), diff(range(Boundary[, 2])))
Nxy <- round(Nxy.size / Nxy.scale)

# Make stack for projections
stk.pred <- MakeProjectionGrid(
   nxy = Nxy,
   mesh = Mesh$mesh,
   data = variables_no_BG,
   tag = "pred",
   boundary = Boundary
)

save(proj, TZ_outline, ebird_full, atlas_full, Mesh, stk.ip, stk.pred, all_lc, variables_no_BG, file = paste0("/users/jhw538/scratch/TZ_INLA/data/TZ_INLA_model_file_E", round(max.edge, digits = 3), ".RData"))

