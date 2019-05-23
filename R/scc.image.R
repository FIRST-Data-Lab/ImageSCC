#' Simultaneous confidence corridors for mean function of imaging data
#'
#' The function is used to construct SCC for mean function of one group of images or
#' difference between mean functions of two sets of images.
#'
#' @param Ya a matrix of imaging data, each row corresponding to one subject/image.
#' @param Yb an optional matrix containing the second group of imaging data. Default is \code{NULL}. When \code{Yb} is \code{NULL}, a one-group SCC is constructed for the mean function of \code{Ya}, otherwise, a two-group SCC is constructed for the difference between the mean functions of \code{Yb} and \code{Ya}.
#' @param Z a 2-column matrix specifying locations of each pixel/voxel.
#' @param Z.band an optional matrix specifying locations for constructing SCC. Default is \code{NULL}. When \code{Z.band} is \code{NULL}, the SCC is evaluated on sample locations provided by matrix \code{Z}.
#' @param d.est degree of bivariate spline for estimating mean function, default is 5.
#' @param d.band degree of bivariate spline for constructing SCC, default is 2.
#' @param r smoothness parameter, default is 1.
#' @param V.est.a the 2-column matrix of vertices' coordinates in the triangulation for estimating mean function of the first set of imaging data.
#' @param Tr.est.a the 3-column matrix of indices of the vertices of triangles in the triangulation.
#' @param V.est.b,Tr.est.b optional information of triangulation used for estimating mean function of the second sample.
#' @param V.band.a,Tr.band.a information of triangultaion for constructing SCC of first set of imaging data.
#' @param V.band.b,Tr.band.b optional information of triangulation for constructing SCC of second set of imaging data.
#' @param penalty logical value indicating whether bivariate penalize spline should be implemented. Default is TRUE.
#' @param lambda the vector of the candidates of penalty parameter.
#' @param alpha.grid vector of confidence levels. Default is \code{c(0.1,0.05,0.01)}.
#' @param adjust.sigma a logical value indicating whether \eqn{\sigma(z)} is adjusted when constructing SCC. Default is TRUE.
#'
#' @importFrom BPST basis
#' @importFrom RSpectra eigs
#' @importFrom pracma Real
#'
#' @details This R package is the implementation program for manuscript entitled "Simultaneous Confidence Corridors for Mean Functions in Functional Data Analysis of Imaging Data" by Yueying Wang, Guannan Wang, Li Wang and R. Todd Ogden.
#'
#' @examples
#' # Triangulation information;
#' data(Brain.V1); data(Brain.Tr1); # triangulation No. 1;
#' data(Brain.V2); data(Brain.Tr2); # triangulation No. 2;
#' V.est=Brain.V2; Tr.est=Brain.Tr2;
#' V.band=Brain.V1; Tr.band=Brain.Tr1;
#' # Location information;
#' n1=40; n2=40;
#' npix=n1*n2
#' u1=seq(0,1,length.out=n1)
#' v1=seq(0,1,length.out=n2)
#' uu=rep(u1,each=n2)
#' vv=rep(v1,times=n1)
#' Z=as.matrix(cbind(uu,vv))
#' ind.inside=inVT(V.est,Tr.est,Z[,1],Z[,2])$ind.inside
#' # Parameters for bivariate spline over triangulation;
#' d.est=5; d.band=2; r=1;
#'
#' # Example 1. One-group SCC;
#' # simulation parameters
#' n=50; lam1=0.5; lam2=0.2; mu.func=2; noise.type='Func';
#' lambda=10^{seq(-6,3,0.5)}; alpha.grid=c(0.1,0.05,0.01);
#' dat=data1g.image(n,Z,ind.inside,mu.func,noise.type,lam1,lam2)
#' Y=dat$Y; beta.true=dat$beta.true;
#' out1=scc.image(Ya=Y,Z=Z,V.est.a=V.est,Tr.est.a=Tr.est,V.band.a=V.band,Tr.band.a=Tr.band,d.est=d.est,d.band=d.band,r=r,penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)
#' scc=out1$scc
#' sum((scc[,1,2]<beta.true[ind.inside]) & (scc[,2,2]>beta.true[ind.inside]))/length(ind.inside)
#' plot(out1)
#'
#' # Example 2. Two-group SCC;
#' # simulation parameters
#' na=50; nb=60; lam1=0.5; lam2=0.2; mu1.func=1; delta=0.3;
#' noise.type='Func'; lambda=10^{seq(-6,3,0.5)}; alpha.grid=c(0.10,0.05,0.01);
#' dat=data2g.image(na,nb,Z,ind.inside,mu1.func,noise.type,lam1,lam2,delta)
#' Ya=dat$Ya; Yb=dat$Yb; beta.true=dat$beta.true;
#' beta.diff=beta.true[,2]-beta.true[,1]
#' V.est.a=V.est.b=V.est;
#' Tr.est.a=Tr.est.b=Tr.est;
#' V.band.a=V.band.b=V.band;
#' Tr.band.a=Tr.band.b=Tr.band;
#' out2=scc.image(Ya=Ya,Yb=Yb,Z=Z,V.est.a=V.est.a,Tr.est.a=Tr.est.a,V.band.a=V.band.a,Tr.band.a=Tr.band.a,V.est.b=V.est.b,Tr.est.b=Tr.est.b,V.band.b=V.band.b,Tr.band.b=Tr.band.b,d.est=d.est,d.band=d.band,r=r,penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)
#' scc=out2$scc
#' sum((scc[,1,2]<beta.diff[ind.inside]) & (scc[,2,2]>beta.diff[ind.inside]))/length(ind.inside)
#' plot(out2)
#'
#' @export
scc.image <- function(Ya,Yb=NULL,Z,Z.band=NULL,d.est=5,d.band=2,r,V.est.a,Tr.est.a,V.band.a,Tr.band.a,
                      V.est.b=NULL,Tr.est.b=NULL,V.band.b=NULL,Tr.band.b=NULL,
                      penalty=TRUE,lambda,alpha.grid=c(0.1,0.05,0.01),adjust.sigma=TRUE){
  try(if(is.null(Ya)) stop("The response images are missing."))
  try(if(nrow(Z)!=ncol(Ya)) stop("The numbers of pixels/voxels in the response images and the location matrix are different."))

  this.call <- match.call()
  mfit <- list()
  # One Group SCC
  if(is.null(Yb)){
    # Basis for estimating mu: Tr.est,V.est,d.est,r
    Bfull.est <- basis(V.est.a,Tr.est.a,d.est,r,Z)
    B.est <- Bfull.est$B
    ind.inside.est <- Bfull.est$Ind.inside
    Q2.est <- Bfull.est$Q2
    K.est <- Bfull.est$K

    # Basis for estimating eta: Tr.band,V.band,d.band,r
    Bfull.band <- basis(V.band.a,Tr.band.a,d.band,r,Z)
    B.band <- Bfull.band$B
    ind.inside.band <- Bfull.band$Ind.inside
    Q2.band <- Bfull.band$Q2
    K.band <- Bfull.band$K

    ind.inside <- intersect(ind.inside.est,ind.inside.band)
    B.est <- B.est[match(ind.inside,ind.inside.est),]
    B.band <- B.band[match(ind.inside,ind.inside.band),]

    if(is.null(Z.band)){
      Z.band <- Z
      B.cover1 <- B.est; B.cover2 <- B.band;
      ind.inside.cover <- ind.inside
    }else{
      Bfull.cover1 <- basis(V.est.a,Tr.est.a,d.est,r,Z.band)
      B.cover1 <- Bfull.cover1$B
      ind.inside.cover1 <- Bfull.cover1$Ind.inside
      Bfull.cover2 <- basis(V.band.a,Tr.band.a,d.band,r,Z.band)
      B.cover2 <- Bfull.cover2$B
      ind.inside.cover2 <- Bfull.cover2$Ind.inside
      ind.inside.cover <- intersect(ind.inside.cover1,ind.inside.cover2)
      B.cover1 <- B.cover1[match(ind.inside.cover1,ind.inside.cover),]
      B.cover2 <- B.cover2[match(ind.inside.cover2,ind.inside.cover),]
    }

    if(!all(paste(round(Z.band[,1],6),round(Z.band[,2],6)) %in%
           paste(round(Z[,1],6),round(Z[,2],6)))){
      adjust.sigma <- FALSE
      warning("Not all the pixels/voxels of constructing SCC are observed, thus the sigma are not adjusted.")
    }
    out <- scc1g.image(Ya,Z,Z.band,d.est,d.band,r,B.est,Q2.est,K.est,B.band,Q2.band,K.band,ind.inside,
                    B.cover1,B.cover2,ind.inside.cover,penalty,lambda,alpha.grid,adjust.sigma)
    mfit$Yhat <- out$mu.hat
    mfit$scc <- out$scc; mfit$cover.zero <- out$cover.zero; mfit$bw <- out$bw;
  }

  if (!is.null(Yb)){
    if(is.null(V.est.b) | is.null(Tr.est.b)){
      V.est.b <- V.est.a
      Tr.est.b <- Tr.est.a
    }
    if(is.null(V.band.b) | is.null(Tr.band.b)){
      V.band.b <- V.band.a
      Tr.band.b <- Tr.band.a
    }

    # Basis for estimating mu: Tr.est,V.est,d.est,r
    Bfull.est.a <- basis(V.est.a,Tr.est.a,d.est,r,Z)
    B.est.a <- Bfull.est.a$B
    ind.inside.est.a <- Bfull.est.a$Ind.inside
    Q2.est.a <- Bfull.est.a$Q2
    K.est.a <- Bfull.est.a$K

    Bfull.est.b <- basis(V.est.b,Tr.est.b,d.est,r,Z)
    B.est.b <- Bfull.est.b$B
    ind.inside.est.b <- Bfull.est.b$Ind.inside
    Q2.est.b <- Bfull.est.b$Q2
    K.est.b <- Bfull.est.b$K

    # Basis for estimation of eta: Tr.band,V.band,d.band,r
    Bfull.band.a <- basis(V.band.a,Tr.band.a,d.band,r,Z)
    B.band.a <- Bfull.band.a$B
    ind.inside.band.a <- Bfull.band.a$Ind.inside
    Q2.band.a <- Bfull.band.a$Q2
    K.band.a <- Bfull.band.a$K

    Bfull.band.b <- basis(V.band.b,Tr.band.b,d.band,r,Z)
    B.band.b <- Bfull.band.b$B
    ind.inside.band.b <- Bfull.band.b$Ind.inside
    Q2.band.b <- Bfull.band.b$Q2
    K.band.b <- Bfull.band.b$K

    ind.inside.est <- intersect(ind.inside.est.a,ind.inside.est.b)
    ind.inside.band <- intersect(ind.inside.band.a,ind.inside.band.b)
    ind.inside <- intersect(ind.inside.est,ind.inside.band)
    B.est.a <- B.est.a[match(ind.inside,ind.inside.est.a),]
    B.est.b <- B.est.b[match(ind.inside,ind.inside.est.b),]
    B.band.a <- B.band.a[match(ind.inside,ind.inside.band.a),]
    B.band.b <- B.band.b[match(ind.inside,ind.inside.band.b),]

    if(is.null(Z.band)){
      Z.band <- Z
      B.cover1.a <- B.est.a; B.cover1.b <- B.est.b;
      B.cover2.a <- B.band.a; B.cover2.b <- B.band.b;
      ind.inside.cover <- ind.inside
    }else{
      Bfull.cover1.a <- basis(V.est.a,Tr.est.a,d.est,r,Z.band)
      B.cover1.a <- Bfull.cover1.a$B
      ind.inside.cover1.a <- Bfull.cover1.a$Ind.inside
      Bfull.cover1.b <- basis(V.est.b,Tr.est.b,d.est,r,Z.band)
      B.cover1.b <- Bfull.cover1.b$B
      ind.inside.cover1.b <- Bfull.cover1.b$Ind.inside
      ind.inside.cover1 <- intersect(ind.inside.cover1.a,ind.inside.cover1.b)

      Bfull.cover2.a <- basis(V.band.a,Tr.band.a,d.band,r,Z.band)
      B.cover2.a <- Bfull.cover2.a$B
      ind.inside.cover2.a <- Bfull.cover2.a$Ind.inside
      Bfull.cover2.b <- basis(V.band.b,Tr.band.b,d.band,r,Z.band)
      B.cover2.b <- Bfull.cover2.b$B
      ind.inside.cover2.b <- Bfull.cover2.b$Ind.inside
      ind.inside.cover2 <- intersect(ind.inside.cover2.a,ind.inside.cover2.b)
      ind.inside.cover <- intersect(ind.inside.cover1,ind.inside.cover2)

      B.cover1.a <- B.cover1.a[match(ind.inside.cover,ind.inside.cover1.a),]
      B.cover1.b <- B.cover1.b[match(ind.inside.cover,ind.inside.cover1.b),]
      B.cover2.a <- B.cover2.a[match(ind.inside.cover,ind.inside.cover2.a),]
      B.cover2.b <- B.cover2.b[match(ind.inside.cover,ind.inside.cover2.b),]
    }

    if(!all(paste(round(Z.band[,1],6),round(Z.band[,2],6)) %in%
            paste(round(Z[,1],6),round(Z[,2],6)))){
      adjust.sigma <- FALSE
      warning("Not all the pixels/voxels of constructing SCC are observed, thus the sigma are not adjusted.")
    }
    out <- scc2g.image(Ya,Yb,Z,Z.band,d.est,d.band,r,
                    B.est.a,Q2.est.a,K.est.a,B.est.b,Q2.est.b,K.est.b,
                    B.band.a,Q2.band.a,K.band.a,B.band.b,Q2.band.b,K.band.b,ind.inside,
                    B.cover1.a,B.cover2.a,B.cover1.b,B.cover2.b,ind.inside.cover,
                    penalty,lambda,alpha.grid,adjust.sigma)
      mfit$Yhat <- rbind(out$mu.hat.a,out$mu.hat.b)
      mfit$scc <- out$scc; mfit$cover.zero <- out$cover.zero; mfit$bw <- out$bw;
  }
  mfit$Ya <- Ya; mfit$V.est.a <- V.est.a; mfit$Tr.est.a <- Tr.est.a;
  mfit$Yb <- Yb; mfit$V.est.b <- V.est.b; mfit$Tr.est.b <- Tr.est.b;
  mfit$Z <- Z; mfit$Z.band <- Z.band;
  mfit$d.est <- d.est; mfit$r <- r; mfit$ind.inside <- ind.inside;
  mfit$ind.inside.cover <- ind.inside.cover
  mfit$alpha <- alpha.grid
  mfit$call <- this.call;
  class(mfit) <- "image"
  return(mfit)
}



