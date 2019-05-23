#' Estimate the mean function via bivariate penalized spline over triangulation
#'
#' This function is used to fit the mean function of the imaging data via  bivariate penalized spline over triangulation. The tuning parameter is selected by generalized cross validation.
#'
#' @param Y a matrix of imaging data, each row corresponding to one subject/image.
#' @param Z a 2-column matrix specifying locations of each pixel/voxel.
#' @param V.est the 2-column matrix of vertices' coordinates in the triangulation for estimating mean function of the first set of imaging data.
#' @param Tr.est the 3-column matrix of indices of the vertices of triangles in the triangulation.
#' @param d.est degree of bivariate spline, default is 5.
#' @param r smoothness parameter. Default is 1.
#' @param lambda candidate of the penalty parameter. Default is \code{10^(-6:3)}.
#' @param proj.matrix a logical value indicating whether the projection matrix
#' will be returned for adjusting \eqn{\sigma(z)} in the construction of SCC, default is FALSE.
#'
#' @importFrom BPST basis
#'
#' @details This R package is the implementation program for manuscript entitled "Simultaneous Confidence Corridors for Mean Functions in Functional Data Analysis of Imaging Data" by Yueying Wang, Guannan Wang, Li Wang and R. Todd Ogden.
#'
#' @examples
#' # Triangulation information;
#' data(Brain.V2); data(Brain.Tr2); # triangulation No. 2;
#' V.est=Brain.V2; Tr.est=Brain.Tr2;
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
#' d.est=5; r=1;
#' # simulation parameters
#' n=50; lam1=0.5; lam2=0.2; mu.func=2; noise.type='Func';
#' lambda=10^{seq(-6,3,0.5)}
#' dat=data1g.image(n,Z,ind.inside,mu.func,noise.type,lam1,lam2)
#' Y=dat$Y
#' Ym=matrix(apply(Y,2,mean),nrow=1)
#' mfit=fit.image(Ym,Z,V.est,Tr.est,d.est,r,lambda)
#' plot(mfit)
#'
#' @export
fit.image <- function(Y,Z,V.est,Tr.est,d.est=5,r=1,lambda=10^(-6:3),proj.matrix=FALSE){

  n <- nrow(Y)
  # Basis for estimating mu: Tr.est,V.est,d.est,r
  Bfull.est <- basis(V.est,Tr.est,d.est,r,Z)
  B <- Bfull.est$B
  ind.inside <- Bfull.est$Ind.inside
  Q2 <- Bfull.est$Q2
  K <- Bfull.est$K
  Y <- matrix(Y[,ind.inside],nrow=n)
  lambda <- as.matrix(lambda)

  this.call <- match.call()
  n <- nrow(Y)
  npix <- ncol(Y)
  J <- ncol(Q2)

  W <- as.matrix(B%*%Q2)
  WW <- crossprod(W,W)
  rhs <- crossprod(W,t(Y))
  D <- crossprod(t(crossprod(Q2,as.matrix(K))),Q2)
  D <- as.matrix(D)

  flag <- (rankMatrix(WW)<J)
  if(!flag){
    Ainv <- chol(WW,pivot=TRUE)
    A <- solve(t(Ainv))
    ADA <- A%*%D%*%t(A)
    eigs <- eigen(ADA)
    Cval <- eigs$values
  }

  nl <- length(lambda)
  gcv_all <- sapply(lambda,FUN=function(Lam){  # eventually a n*nl matrix
      Dlam <- Lam*D
      lhs <- WW+Dlam
      lhs.inv <- chol2inv(chol(lhs));
      theta <- crossprod(t(lhs.inv),rhs)
      gamma <- crossprod(t(Q2),theta)
      Yhat <- crossprod(t(W),theta)
      res <- t(Y)-Yhat
      sse <- apply(res^2,2,sum)
      if(!flag){
        df <- sum(1/(1+Cval*Lam))
      }
      if(flag){
        Hmtx <- crossprod(t(crossprod(t(W),lhs.inv)),t(W))
        df <- sum(diag(Hmtx))
      }
      gcv <- npix*sse/(npix-df)^2
  })
  gcv_all <- matrix(gcv_all,nrow=n)
  lam.ind <- apply(gcv_all,1,which.min)
  lambdac <- lambda[lam.ind]
  theta <- c()
  gamma <- c()
  Yhat <- c()
  df <- c()
  for (i in 1:n){
    lamc.tmp <- lambdac[i]
    Dlam <- lamc.tmp*D
    lhs <- WW+Dlam
    lhs.inv <- chol2inv(chol(lhs));
    rhs.tmp <- as.matrix(rhs[,i],ncol=1)
    theta.tmp <- crossprod(t(lhs.inv),rhs.tmp)
    theta <- cbind(theta,theta.tmp)
    gamma.tmp <- crossprod(t(Q2),theta.tmp)
    gamma <- cbind(gamma,gamma.tmp)
    Yhat.tmp <- crossprod(t(W),theta.tmp)
    Yhat <- cbind(Yhat,Yhat.tmp)
    if(!flag){
      df.tmp <- sum(1/(1+Cval*lamc.tmp))
    }
    if(flag){
      Hmtx <- crossprod(t(crossprod(t(W),lhs.inv)),t(W))
      df.tmp <- sum(diag(Hmtx))
    }
    df <- c(df,df.tmp)
  }
  if(proj.matrix==TRUE){   # projection matrix
    if(n==1){
      Dlam <- lambdac*D
      lhs <- WW+Dlam
      lhs.inv <- chol2inv(chol(lhs))
      Slam <- crossprod(t(W),lhs.inv%*%t(W))
    }
    if(n>1){
      Slam <- lapply(lambdac,FUN=function(lamc){
        Dlam <- lamc*D
        lhs <- WW+Dlam
        lhs.inv <- chol2inv(chol(lhs))
        Slam <- crossprod(t(W),lhs.inv%*%t(W))
        return(Slam)
      })
    }
  }else{Slam <- NA}
  Yhat <- t(Yhat)
  res <- Y-Yhat
  sse <- apply(res^2,1,sum)
  gcv <- npix*sse/(npix-df)^2
  bic <- log(sse/npix)+df*log(npix)/npix
  mfit <- list()
  mfit$theta <- theta; mfit$gamma <- gamma;
  mfit$Yhat <- Yhat; mfit$res <- res; mfit$sse <- sse;
  mfit$gcv <- gcv; mfit$bic <- bic;
  mfit$lambdac <- lambdac; mfit$Slam <- Slam;
  mfit$scc <- NULL;
  mfit$Y <- Y; mfit$Z <- Z; mfit$V.est.a <- V.est; mfit$Tr.est.a <- Tr.est;
  mfit$V.est.b <- NULL; mfit$Tr.est.b <- NULL;
  mfit$d.est <- d.est; mfit$r <- r; mfit$ind.inside <- ind.inside;
  mfit$call <- this.call;
  class(mfit) <- "image"
  return(mfit)
}
