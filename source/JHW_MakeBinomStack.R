#' Function to create stack for presence only points
#'
#' @param data SpatialPointsDataFrame of covariates.
#' @param observs SpatialPoints object of presences 
#' @param tag Name for tag for the stack (defaults to "points").
#' @param intercept Boolean, should an intercept be added? Defaults to TRUE.
#' @param mesh INLA mesh.
#' @param presname Names of presences column in observs. Defaults to "NPres".
#' Note that this column can also be logical.
#' @param coordnames Names of coordinates.
#' @param InclCoords Boolean, shoiuld coordinates be included in data (defaults to FALSE).
#' @param polynoms If not NULL, a SpatialPolygons object, with (for example) range maps.
#' @param scale Should the distance be scaled by dividing by the mean of the non-zero distances?
#' Defaults to FALSE, either logical or numeric. Ignored if polynoms is NULL.
#'
#' @return An INLA stack with binomial data
#'
#' @export
#' @import INLA

MakeBinomStack=function(data, observs, tag="points", intercept=TRUE, mesh, presname="NPres", 
                        coordnames=NULL, InclCoords=FALSE, polynoms = NULL, scale = FALSE, add_effort = FALSE, effort_names=NULL) {

  if(length(presname)>1) stop("more than one name given for presences column")


  if(is.null(coordnames)) coordnames <- colnames(data@coords)
  if(!is.null(polynoms)) {
    if(class(polynoms) != "SpatialPolygonsDataFrame" & class(polynoms) != "SpatialPolygons")
      stop("polynoms should be a spatial polygon")
  }

  NearestCovs <- GetNearestCovariate(points=observs, covs=data)
  
  if(InclCoords) {    NearestCovs@data[,coordnames] <- data@coords  }
  if(intercept){ NearestCovs@data[,paste("int",tag,sep=".")] <- 1 }# add intercept
  if(add_effort){ NearestCovs@data[, effort_names] <- observs@data[, effort_names]}
  if(!is.null(polynoms)) {
    NearestCovs <- AddDistToRangeToSpatialPoints(data = NearestCovs, polynoms = polynoms, scale=scale)
  }

# If presences are Boolean, reformat
  if(is.logical(observs@data[,presname])) {
    observs@data[,presname] <- as.integer(observs@data[,presname])
  }
  # Projector matrix from mesh to data.
  projmat <- inla.spde.make.A(mesh, as.matrix(NearestCovs@coords)) # from mesh to point observations

  stk.binom <- inla.stack(data=list(resp=observs@data[,presname]), A=list(1,projmat), tag=tag,
                          effects=list(NearestCovs@data, list(i=1:mesh$n)))

  stk.binom
}
