sortBase <- function(vec, n.knots = 2) {
   ## Function to calculate bases for regression splines. Modified from code
   ## provided in Crainiceanu, C., Ruppert, D. & Wand, M.P. Bayesian analysis for
   ## penalized spline regression using WinBUGS. J. Stat. Soft. 14, 1 24(2005).
   ## Parameter vec is a vector defining the raw data vector, n.knots defines the
   ## number of knots in the GAM.
   N         <- length(vec)
   x.time    <- c(vec)
   zFE       <- cbind(rep(1,N), x.time)
   x.knots   <- quantile(unique(x.time), seq(0, 1, length = (n.knots+2))[-c(1,
                                                                            (n.knots+2))], na.rm = TRUE)
   z_K       <- (abs(outer(x.time,x.knots,"-")))^3
   OMEGA.all <- (abs(outer(x.knots,x.knots,"-")))^3
   svd.OMEGA.all  <- svd(OMEGA.all)
   sqrt.OMEGA.all <- t(svd.OMEGA.all$v %*% (t(svd.OMEGA.all$u) *
                                               sqrt(svd.OMEGA.all$d)))
   z.out     <- t(solve(sqrt.OMEGA.all, t(z_K)))
   return(z.out)
}


make_ebird_sp <- function(scientific_name, ROI){
   require(spatialEco)
      df <- ebird_filtered %>% 
            group_by(LATITUDE, LONGITUDE, `SAMPLING EVENT IDENTIFIER`, `DURATION MINUTES`, 
                     `EFFORT DISTANCE KM`, `NUMBER OBSERVERS`, `OBSERVATION DATE`, `LOCALITY`) %>% 
            summarise(occurrence = ifelse(scientific_name %in% `SCIENTIFIC NAME`, TRUE, FALSE)) %>% 
            ungroup() %>% 
            group_by(LATITUDE, LONGITUDE, `DURATION MINUTES`, 
                     `EFFORT DISTANCE KM`, `NUMBER OBSERVERS`, `OBSERVATION DATE`, `LOCALITY`) %>% 
            slice_head() %>%    # Where duplicate checklists occurred, keep only the first
            ungroup() %>% 
            rename(duration_minutes = `DURATION MINUTES`,
                   effort_distance_km = `EFFORT DISTANCE KM`,
                   number_observers = `NUMBER OBSERVERS`)
      
      cat('Number of eBird records: ', length(rownames(df[df$occurrence == TRUE, ])))
      
      df_sp <- SpatialPointsDataFrame(
            coords = df[, c("LONGITUDE", "LATITUDE")],
            data = data.frame(presence = df$occurrence, duration_minutes = df$duration_minutes,
                              effort_distance_km = df$effort_distance_km, number_observers = df$number_observers),
            proj4string = crs(proj)
      )
      
      # Only include eBird data points for the region of interest
      # Get intersecting points
      in_sp <- rgeos::gIntersection(df_sp, ROI)
      
      # Only keep intersecting points in original spdf
      df_sp_ROI <- df_sp[in_sp, ]
      cat('Number of eBird records: ', length(rownames(df_sp_ROI@data[df_sp_ROI@data$presence == TRUE, ])))
      
      
      return(df_sp_ROI)
}

get_hyperpars <- function(path){
   out = NULL
   list <- list.files(path = path)
   
   for(i in 1:NROW(list)){
      file = list[i]
      
      load(paste0(path, file))
      
      # Check hyperparameters. sd ideally should be somewhat smaller than mean. Stdev should not be very small
      outdf <- as.data.frame(model[["model"]][["summary.hyperpar"]])
      outdf$par <- c("Range", "Stdev")
      outdf <- cbind(outdf, rbind(data, data))
      out <-  rbind(out, outdf)
   }
   return(out)
}

plot_INLA_predictions <- function(model_out){
   load(model_out)
   cloglog_inv <- function(x){
      1-exp(-exp(x))
   }
   
   proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
   
   Pred <- SpatialPixelsDataFrame( 
      points = stk.pred$predcoords,
      data = model$predictions,
      proj4string = crs(proj)
   )
   
   Pred@data$inv_mean <- cloglog_inv(Pred@data$mean)
   # test <- raster(Pred)
   # writeRaster(test, "TZ_fsl_predictions.tif")
   # Plot of predictions
   plot_title <- str_split(model_out, pattern = "/")[[1]][length(str_split(model_out, pattern = "/")[[1]])]

   p1 <- ggplot() +
      gg(Pred, aes(fill = inv_mean)) +
      coord_equal() +
      theme_void() +
      ggtitle(plot_title) + 
      labs(fill = "Probability") 
   return(p1)
}



plot_INLA_random <- function(model_out, ROI){
   load(model_out)
   # Visualize the random field, the unexplained spatial error
   projector <- inla.mesh.projector(Mesh$mesh)
   projection <- inla.mesh.project(projector,
                                   model[["model"]][["summary.random"]]$i$mean)
   plot_title <- str_split(model_out, pattern = "/")[[1]][length(str_split(model_out, pattern = "/")[[1]])]
   p1 <- INLAutils::ggplot_projection_shapefile(projection, projector) +
      geom_polygon(data = ROI, aes(long, lat, group = group), fill = NA, col = 'black') +
      coord_equal() +
      ggtitle(plot_title)
   
   return(p1)
}

plot_INLA_fit <- function(model_out){
   require(cowplot)
   require(gridGraphics)
   load(model_out)
   
   failure <- sum(model$model$cpo$failure) 
   print(paste0("Number of CPO failures: ", failure))
   data <- data.frame(pit = model$model$cpo$pit)
   p1 <- ggplot(data) +
      geom_histogram(aes(x = pit))
   # hist(model$model$cpo$pit, main="", breaks = 10, xlab = "PIT")
   # p1 <- recordPlot()  
   # 
   # plot.ecdf(model$model$cpo$pit); abline(0, 1)
   # p2 <- recordPlot()  
   # 
   # p3 <- plot_grid(p1, p2, ncol = 2)
   return(p1)
}

prepare_GAM <- function(df, colname){
   ## create sequences for predictions:
   seq <- seq(min(df@data[, colname], na.rm = TRUE), max(df@data[, colname], na.rm = TRUE), length = 100)
   # add prediction vector and scale numeric values to mean 0 and sd 1
   var.s <- c(scale(c(seq, df@data[, colname])))  
   z.var.s <- scale(sortBase(var.s, n.knots = 2))  
   
   df@data[, paste0("z.", colname, "1.s")] <- c(z.var.s[101:NROW(z.var.s), 1])
   df@data[, paste0("z.", colname, "2.s")] <- c(z.var.s[101:NROW(z.var.s), 2])
   
   # separate out actual data and add to df
   df@data[, paste0( colname, ".s")]  <- var.s[-c(1:100)]
   
   ## Separate out prediction vector
   assign(paste0("z.", colname, "1.s"), c(z.var.s[1:100, 1]), envir = .GlobalEnv)
   assign(paste0("z.", colname, "2.s"), c(z.var.s[1:100, 2]), envir = .GlobalEnv)
   
   return(df)
}


