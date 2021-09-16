#' Function to get covariate values nearest the data/integration point
#'
#' @param points Points for which we need covariates. Either a SpatialPoints* object or a matrix with coordinates in columns.
#' @param covs Covariates, with coordinates in first and second columns, or a SpatialPointsDataFrame object.
#' @return A SpatialPointsDataFrame with the original points plus the data from the closest point.
#'
#' If there are ties, this will use the first element.
#'
#' @export
#' 


# Create data set that exists at center of triangles



GetNearestCovariate <- function(points, covs) {
  
  #if covariate and point objects are not spatial objects, convert them into spatial obkjects
  if(class(covs)!="SpatialPointsDataFrame") {
    covs <- SpatialPointsDataFrame(coords = covs[,1:2], data = covs[,-(1:2)])
  }
  
  if(class(points)!="SpatialPointsDataFrame") {
    points <- SpatialPoints(points)
  }

  #calculate distance between each mesh centroid and each covariate location, 
  #and return the row index of the closest covariate
  Nearest <- apply(points@coords, 1, function(pt, pts) {   # points are covariate locations
    dists <- spDistsN1(pts, pt, longlat = TRUE)
    which(dists==min(dists))[1] # need [1] in case min() is not unique. But if it's not, doesn't matter
  }, pts=covs) 

  #create spatial point data frame with coordinates of mesh centroids and associated covariate values for 
  #nearest cpovariate location
  points.df <- SpatialPointsDataFrame(coords=points@coords,
                                      data = covs@data[Nearest, , drop=FALSE],
                                      proj4string = CRS(proj4string(covs)))
  return(points.df)
}
