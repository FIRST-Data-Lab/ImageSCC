#' Penalized least square estimation with GCV
#'
#' The function selects penalty parameter via generalized
#' cross validation and solves penalized least square problem
#' in estimating mean functions with bivariate penalized splines.
#'
#' @param B bivariate spline basis matrix.
#' @param Q2 qr decomposition of the smoothness matrix.
#' @param K energy matrix.
#' @param lambda candidate of the penalty parameter.
#' @param Y a matrix of data with each row corresponding to one subject/image.
#' @param proj.matrix a logical value indicating whether the projection matrix
#' will be returned for adjusting \eqn{\sigma(z)} in the construction of SCC.
#'
#' @details
#' @export
fit.mean <- function(B,Q2,K,lambda=10^(-6:3),Y,proj.matrix=FALSE){
  if(!is.matrix(Y)){
    warning("The response variable, Y, should be a matrix with each row represents an image.")
    Y <- as.matrix(Y)
  }
  lambda <- as.matrix(lambda)

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

  list(theta=theta,gamma=gamma,Yhat=Yhat,res=res,sse=sse,
       gcv=gcv,bic=bic,lambdac=lambdac,Slam=Slam)
}
