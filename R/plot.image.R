#' Plot image with triangulation
#'
#' The plot method for class 'image'. The function plots a series of images of triangulations used for mean estimation, estimated mean functions and simultaneous confidence corridors.
#'
#' @param mfit an "image" object returned from either \code{fit.image} or \code{scc.image}.
#'
#' @importFrom Triangulation TriPlot
#' @importFrom graphics image
#'
#' @details This R package is the implementation program for manuscript entitled "Simultaneous Confidence Corridors for Mean Functions in Functional Data Analysis of Imaging Data" by Yueying Wang, Guannan Wang, Li Wang and R. Todd Ogden.
#'
#' @export
plot.image <- function(mfit){
  # Plot 1. Triangulation
  TriPlot(mfit$V.est.a,mfit$Tr.est.a)
  readline(prompt="Press [enter] to continue")

  if(!is.null(mfit$V.est.b) & !is.null(mfit$Tr.est.b)){
    TriPlot(mfit$V.est.b,mfit$Tr.est.b)
    readline(prompt="Press [enter] to continue")
  }

  # Plot 2. Estimated mean function
  Z <- matrix(mfit$Z,ncol=2)
  z1 <- unique(Z[,1]); z2 <- unique(Z[,2]);
  n1 <- length(z1); n2 <- length(z2);
  ind.inside <- mfit$ind.inside
  for(ii in 1:nrow(mfit$Yhat)){
    mu.hat <- rep(NA,n1*n2)
    mu.hat[ind.inside] <- mfit$Yhat[ii,]
    mu.hat.mtx <- matrix(mu.hat,nrow=n2,ncol=n1)
    image(z2,z1,mu.hat.mtx)
    readline(prompt="Press [enter] to continue")
  }

  # Plot 3. SCC
  if(!is.null(mfit$scc)){
    Z.band <- matrix(mfit$Z.band,ncol=2)
    z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]);
    n1 <- length(z1); n2 <- length(z2);
    ind.inside.band <- mfit$ind.inside.cover
    nalpha <- (dim(mfit$scc))[3]
    for(ii in 1:nalpha){
      scc <- matrix(NA,n1*n2,2)
      scc[ind.inside.band,] <- mfit$scc[,,ii]
      scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1)
      image(z2,z1,scc.l.mtx)
      title(paste("Lower SCC when alpha =",mfit$alpha[ii]))
      readline(prompt="Press [enter] to continue")
      scc.u.mtx <- matrix(scc[,2],nrow=n2,ncol=n1)
      image(z2,z1,scc.u.mtx)
      title(paste("Upper SCC when alpha =",mfit$alpha[ii]))
      readline(prompt="Press [enter] to continue")
    }
  }
}
