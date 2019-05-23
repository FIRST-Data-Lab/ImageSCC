#' Simultaneous confidence corridors for mean function of imaging data: one sample case
#' @importFrom RSpectra eigs
#' @importFrom pracma Real
#'
scc1g.image <- function(Y,Z,Z.band,d.est,d.band,r,B.est,Q2.est,K.est,B.band,Q2.band,K.band,ind.inside,
                        B.cover1,B.cover2,ind.inside.cover,
                        penalty=FALSE,lambda,alpha.grid=c(0.1,0.05,0.01),adjust.sigma=TRUE){
  n <- nrow(Y)
  npix <- length(ind.inside)
  npix.cover <- length(ind.inside.cover)
  nalpha <- length(alpha.grid)

  # Step 1. Estimate mu.hat
  Ym <- matrix(apply(Y,2,mean),nrow=1)
  Y <- Y[,ind.inside]
  Ym <- matrix(Ym[,ind.inside],nrow=1)
  mfit0 <- fit.mean(B.est,Q2.est,K.est,lambda,Ym,proj.matrix=TRUE)

  gamma.hat0 <- mfit0$gamma
  mu.hat <- mfit0$Yhat
  mu.hat.mtx <- matrix(rep(drop(mu.hat),times=n),nrow=n,byrow=TRUE)
  mu.hat.cover <- B.cover1%*%gamma.hat0

  R0 <- Y-mu.hat.mtx

  # Step 2.1. Estimate eta.hat & Geta.hat;
  if(d.band==-1){# using piecewise constant
    BB.band <- crossprod(B.band)
    BB.band.inv <- chol2inv(chol(BB.band))
    BR.band <- crossprod(B.band,t(R0))
    mfit1 <- B.band%*%BB.band.inv%*%BR.band
    eta.hat <- t(mfit1)

    eta.cover <- B.cover2%*%BB.band.inv%*%BR.band
    eta.cover <- t(eta.cover)
  }
  if(d.band>1){# using higher order spline
    if(penalty==FALSE){
      BQ2.band <- as.matrix(B.band%*%Q2.band)
      BB.band <- crossprod(BQ2.band)
      BB.band.inv <- chol2inv(chol(BB.band))
      BR.band <- crossprod(BQ2.band,t(R0))
      mfit1 <- BQ2.band%*%BB.band.inv%*%BR.band
      eta.hat <- t(mfit1)

      BQ2.cover2 <- as.matrix(B.cover2%*%Q2.band)
      eta.cover <- BQ2.cover2%*%BB.band.inv%*%BR.band
      eta.cover <- t(eta.cover)
    }
    if(penalty==TRUE){
      mfit1 <- fit.mean(B.band,Q2.band,K.band,lambda,R0)
      eta.hat <- mfit1$Yhat
      gamma.hat <- mfit1$gamma
      Beta.hat <- tcrossprod(gamma.hat)/n

      eta.cover <- as.matrix(B.cover2%*%gamma.hat)
      eta.cover <- t(eta.cover)
    }
  }
  Geta.hat <- crossprod(eta.hat)/n
  diag.Geta.hat <- diag(Geta.hat)
  diag.Geta.cover <- apply(eta.cover^2,2,mean)

  # Step 2.2. Estimate sigma.hat
  eps.hat <- R0-eta.hat
  sigma.hat <- apply(eps.hat^2,2,mean)
  sigma.hat <- matrix(sqrt(sigma.hat),ncol=1)

  # Step 3. Get lambda.hat & phi.hat
  Geig <- eigs(Geta.hat,k=10)
  evalues <- Real(Geig$values)
  Kappa <- sum(cumsum(evalues)/sum(evalues)<0.99)+1
  lambdaK.hat <- evalues[1:Kappa]
  psi.hat <- matrix((Geig$vectors)[,1:Kappa],ncol=Kappa) #eigenfunction
  phi.hat <- psi.hat%*%diag(sqrt(lambdaK.hat),nrow=Kappa) #sqrt{lambdaK}*eigenfunction

  # Step 4. Generate qunatile
  nb <- 1000
  if (adjust.sigma==FALSE){
    Zb <- matrix(rnorm(Kappa*nb,0,1),nb,Kappa)
    tmp.hat <- diag(1/sqrt(diag(Geta.hat)))
    zeta.hat <- phi.hat%*%t(Zb)
    zeta.hat <- tmp.hat%*%zeta.hat
    Qalpha.hat <- quantile(apply(abs(zeta.hat),2,max),c(1-alpha.grid))
    Sigma.hat <- NA
  }else if (adjust.sigma==TRUE){
    Proj.est <- mfit0$Slam
    Zk <- matrix(rnorm(Kappa*nb,0,1),nb,Kappa)
    Zeps <- matrix(rnorm(npix*nb,0,1),nb,npix)
    Sigma.hat <- Geta.hat+
      Proj.est%*%diag(as.numeric(sigma.hat)^2)%*%t(Proj.est)
    tmp.hat <- diag(1/sqrt(diag(Sigma.hat)))
    zeta.hat <- phi.hat%*%t(Zk)+ Proj.est%*%t(Zeps%*%diag(as.numeric(sigma.hat)))
    zeta.hat <- tmp.hat%*%zeta.hat
    Qalpha.hat <- quantile(apply(abs(zeta.hat),2,max),probs=c(1-alpha.grid))
  }

  # Step 5. Find the SCC
  if (adjust.sigma==FALSE){
    ME <- matrix(rep(Qalpha.hat,each=npix.cover),nrow=npix.cover,ncol=nalpha)*
      matrix(rep(sqrt(diag.Geta.cover)/sqrt(n),times=nalpha),nrow=npix.cover,ncol=nalpha)
    mu.hat.mtx <- matrix(rep(mu.hat.cover,times=nalpha),nrow=npix.cover,ncol=nalpha)
    scc.l <- mu.hat.mtx-ME
    scc.u <- mu.hat.mtx+ME
    bw <- 2*apply(ME,2,mean)
    eps.cover <- eps.hat
  }else if (adjust.sigma==TRUE){
    match.band <- match(paste(round(Z.band[ind.inside.cover,1],6),round(Z.band[ind.inside.cover,2],6)),
                     paste(round(Z[ind.inside,1],6),round(Z[ind.inside,2],6)))
    eps.cover <- R0[,match.band]-eta.cover
    diag.Sigma.cover <- diag(Sigma.hat)[match.band]
    ME <- matrix(rep(Qalpha.hat,each=npix.cover),nrow=npix.cover,ncol=nalpha)*
      matrix(rep(sqrt(diag.Sigma.cover)/sqrt(n),times=nalpha),nrow=npix.cover,ncol=nalpha)
    mu.hat.mtx <- matrix(rep(mu.hat.cover,times=nalpha),nrow=npix.cover,ncol=nalpha)
    scc.l <- mu.hat.mtx-ME
    scc.u <- mu.hat.mtx+ME
    bw <- 2*apply(ME,2,mean)
  }
  scc <- array(NA,dim=c(npix.cover,2,nalpha))
  for(j in 1:nalpha){
    scc[,,j] <- cbind(scc.l[,j],scc.u[,j])
  }
  # Results
  list(V.est=V.est,Tr.est=Tr.est,V.band=V.band,Tr.band=Tr.band,
       d.est=d.est,d.band=d.band,r=r,mfit=mfit0,
       mu.hat=mu.hat,eta.hat=eta.hat,eps.hat=eps.hat,scc=scc,bw=bw,
       mu.hat.cover=mu.hat.cover,eta.cover=eta.cover,eps.cover=eps.cover)
}


