#' Triangulation selection for contructing the SCC via bootstrap
#'
#' The function selects triangulation used for constructing SCC by bootstrap method.
#'
#' @param Ya a matrix of data with each row corresponding to one subject/image.
#' @param Yb an optional matrix containing the second group of imaging data. When \code{Yb} is NULL, triangulation selection focuses on one sample SCC, otherwise it focuses on the SCC for mean difference (\eqn{\mu_2-\mu_1}) between two sets of images.
#' @param Z a 2-column matrix specifying locations of information.
#' @param d.est degree of bivariate spline for mean estimation.
#' @param d.band degree of bivariate spline for SCC.
#' @param r smoothness parameter.
#' @param V.est.a the 2-column matrix of vertices' coordinates in the triangulation for estimating mean in the first sample.
#' @param Tr.est.a the 3-column matrix specifying triangles in the triangulation.
#' Each row contains 3 indices of vertices corresponding to one triangle in the triangulation.
#' @param V.est.b,Tr.est.b optional information of triangulation used for estimating mean in the second sample.
#' @param V.bands,Tr.bands lists of candidates for triangulations used to construct SCC.
#' @param lambda the vector of the candidates for penalty parameter when estimating mean function.
#' @param nboot number of bootstrap iterations. Default is 50.
#' @param alpha0 a value specifying confidence level of SCC.
#' @param adjust.sigma a logical value indicating whether \eqn{\sigma(z)} is adjusted when constructing SCC.
#' Default is TRUE.
#'
#' @details This R package is the implementation program for manuscript entitled ``Simultaneous Confidence Corridors for Mean Functions in Functional Data Analysis of Imaging Data" by Yueying Wang, Guannan Wang, Li Wang and R. Todd Ogden.
#'
#' @examples
#' # Triangulation information;
#' data(Brain.V1); data(Brain.Tr1); # triangulation No. 1;
#' data(Brain.V2); data(Brain.Tr2); # triangulation No. 2;
#' data(Brain.V3); data(Brain.Tr3); # triangulation No. 3;
#' #' V.est=Brain.V2; Tr.est=Brain.Tr2;
#' V.bands=list(V1=Brain.V1,V2=Brain.V2,V3=Brain.V3);
#' Tr.bands=list(Tr1=Brain.Tr1,Tr2=Brain.Tr2,Tr3=Brain.Tr3);
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
#' tri.band=boot.image(Ya=Y,Z=Z,d.est=d.est,d.band=d.band,r=r,V.est.a=V.est,Tr.est.a=Tr.est,V.bands=V.bands,Tr.bands=Tr.bands,lambda=lambda)
#' tri.band$tri.band
#' V.band=tri.band$V.band.a; Tr.band=tri.band$Tr.band.a;
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
#' tri.band=boot.image(Ya=Ya,Yb=Yb,Z=Z,d.est=d.est,d.band=d.band,r=r,V.est.a=V.est.a,Tr.est.a=Tr.est.a,V.est.b=V.est.b,Tr.est.b=Tr.est.b,V.bands=V.bands,Tr.bands=Tr.bands,lambda=lambda)
#' tri.band$tri.band
#' V.band.a=tri.band$V.band.a; Tr.band.a=tri.band$Tr.band.a;
#' V.band.b=tri.band$V.band.b; Tr.band.b=tri.band$Tr.band.b;
#'
#' @export
boot.image <- function(Ya,Yb=NULL,Z,d.est,d.band,r,V.est.a,Tr.est.a,
                       V.est.b=NULL,Tr.est.b=NULL,V.bands,Tr.bands,
                       lambda,nboot=50,alpha0=0.05,adjust.sigma=TRUE){
  try(if(missing("Ya")) stop("The response images are missing."))
  try(if(nrow(Z)!=ncol(Ya)) stop("The numbers of pixels/voxels in the response images and the location matrix are different."))
  Bfull <- basis.image(Z,V.bands,Tr.bands,d.band,r)
  B.bands <- Bfull$B.ests
  Q2.bands <- Bfull$Q2.ests
  K.bands <- Bfull$K.ests
  ind.inside.band <- Bfull$ind.inside

  Bfull.est.a <- basis(V.est.a,Tr.est.a,d.est,r,Z)
  B.est.a <- Bfull.est.a$B
  ind.inside.est.a <- Bfull.est.a$Ind.inside
  Q2.est.a <- Bfull.est.a$Q2
  K.est.a <- Bfull.est.a$K
  ind.inside <- intersect(ind.inside.est.a,ind.inside.band)

  # One group
  if(is.null(Yb)){
    Ya <- Ya[,ind.inside]
    tri.band <- boot1g.image(Y=Ya,Z=Z,d.band=d.band,r=r,
                          B.est=B.est.a,Q2.est=Q2.est.a,K.est=K.est.a,
                          B.bands=B.bands,Q2.bands=Q2.bands,K.bands=K.bands,ind.inside=ind.inside,
                          lambda=lambda,alpha0=alpha0,adjust.sigma=adjust.sigma)
    B.band.a <- B.bands[[tri.band]]
    Q2.band.a <- Q2.bands[[tri.band]]
    K.band.a <- K.bands[[tri.band]]
    V.band.a <- V.bands[[tri.band]]
    Tr.band.a <- Tr.bands[[tri.band]]
    B.band.b <- NA
    Q2.band.b <- NA
    K.band.b <- NA
    V.band.b <- NA
    Tr.band.b <- NA
  }
  if (!is.null(Yb)){
    if(is.null(V.est.b) | is.null(Tr.est.b)){
      V.est.b <- V.est.a
      Tr.est.b <- Tr.est.a
      B.est.b <- B.est.a
      Q2.est.b <- Q2.est.a
      K.est.b <- K.est.a
    } else {
      Bfull.est.b <- basis(V.est.b,Tr.est.b,d.est,r,Z)
      B.est.b <- Bfull.est.b$B
      ind.inside.est.b <- Bfull.est.b$Ind.inside
      Q2.est.b <- Bfull.est.b$Q2
      K.est.b <- Bfull.est.b$K
      ind.inside <- intersect(ind.inside,ind.inside.est.b)
    }

    Ya <- Ya[,ind.inside]
    Yb <- Yb[,ind.inside]
    tri.band <- boot2g.image(Ya=Ya,Yb=Yb,Z=Z,d.band=d.band,r=r,
                     B.est.a=B.est.a,Q2.est.a=Q2.est.a,K.est.a=K.est.a,
                     B.est.b=B.est.b,Q2.est.b=Q2.est.b,K.est.b=K.est.b,
                     B.bands=B.bands,Q2.bands=Q2.bands,K.bands=K.bands,
                     ind.inside=ind.inside,lambda=lambda,nboot=nboot,
                     alpha0=alpha0,adjust.sigma=adjust.sigma)
    tri.band.a <- tri.band[1]; tri.band.b <- tri.band[2];
    B.band.a <- B.bands[[tri.band.a]]
    Q2.band.a <- Q2.bands[[tri.band.a]]
    K.band.a <- K.bands[[tri.band.a]]
    V.band.a <- V.bands[[tri.band.a]]
    Tr.band.a <- Tr.bands[[tri.band.a]]
    B.band.b <- B.bands[[tri.band.b]]
    Q2.band.b <- Q2.bands[[tri.band.b]]
    K.band.b <- K.bands[[tri.band.b]]
    V.band.b <- V.bands[[tri.band.b]]
    Tr.band.b <- Tr.bands[[tri.band.b]]
  }
  list(tri.band=tri.band,V.band.a=V.band.a,V.band.b=V.band.b,
       Tr.band.a=Tr.band.a,Tr.band.b=Tr.band.b,
       B.band.a=B.band.a,B.band.b=B.band.b,
       Q2.band.a=Q2.band.a,Q2.band.b=Q2.band.b,
       K.band.a=K.band.a,K.band.b=K.band.b)
}
