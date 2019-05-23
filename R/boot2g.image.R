#' Select triangulation for two sample SCC construction
#'
boot2g.image <- function(Ya,Yb,Z,d.band,r,B.est.a,Q2.est.a,K.est.a,
                      B.est.b,Q2.est.b,K.est.b,tri.boot=NULL,
                      B.bands,Q2.bands,K.bands,ind.inside,lambda,
                      nboot=50,alpha0=0.05,adjust.sigma=TRUE){
  na <- nrow(Ya)
  nb <- nrow(Yb)
  rho <- na/nb
  npix <- length(ind.inside)
  nalpha <- 10
  alpha.grid <- seq((alpha0-0.005),(alpha0+0.005),length.out=nalpha)
  ntri <- length(B.bands)

  # Choose triangultion for group 1 and 2 respectively
  Yam <- matrix(apply(Ya,2,mean),nrow=1)
  Ybm <- matrix(apply(Yb,2,mean),nrow=1)

  # Step 1. Estimate mu.hat
  mfit0.a <- fit.mean(B.est.a,Q2.est.a,K.est.a,lambda,Yam,proj.matrix=TRUE)
  lamc.a <- mfit0.a$lambdac
  theta.a <- mfit0.a$theta
  gamma.a <- mfit0.a$gamma
  mu.hat.a <- mfit0.a$Yhat
  mu.hat.mtx.a <- matrix(rep(drop(mu.hat.a),times=na),nrow=na,byrow=TRUE)
  R0.a <- Ya-mu.hat.mtx.a
  Slam.a <- mfit0.a$Slam

  mfit0.b <- fit.mean(B.est.b,Q2.est.b,K.est.b,lambda,Ybm,proj.matrix=TRUE)
  lamc.b <- mfit0.b$lambdac
  theta.b <- mfit0.b$theta
  gamma.b <- mfit0.b$gamma
  mu.hat.b <- mfit0.b$Yhat
  mu.hat.mtx.b <- matrix(rep(drop(mu.hat.b),times=nb),nrow=nb,byrow=TRUE)
  R0.b <- Yb-mu.hat.mtx.b
  Slam.b <- mfit0.b$Slam

  if(is.null(tri.boot)){
    tri.boot <- cbind(a=rep(1:ntri,each=ntri),
                   b=rep(1:ntri,times=ntri))
  }
  ntri.boot <- nrow(tri.boot)
  cover.all <- vector('list',length=ntri.boot)
  if(nboot>=50){
    nblock.boot <- floor(nboot/10)
    nboot0 <- nboot
  }
  for(itri in 1:ntri.boot){
    cover.tri <- c()
    for(i.boot in 1:nblock.boot){
      if(i.boot<nblock.boot){
        nboot <- 10
      }else{
        nboot <- nboot0-(nblock.boot-1)*10
      }
      itri.a <- tri.boot[itri,'a']
      itri.b <- tri.boot[itri,'b']

      B.band.a <- B.bands[[itri.a]]
      Q2.band.a <- Q2.bands[[itri.a]]
      K.band.a <- K.bands[[itri.a]]

      B.band.b <- B.bands[[itri.b]]
      Q2.band.b <- Q2.bands[[itri.b]]
      K.band.b <- K.bands[[itri.b]]

      if(d.band==-1){
        BB.band.a <- crossprod(as.matrix(B.band.a))
        BB.band.inv.a <- chol2inv(chol(BB.band.a))
        HB.band.a <- B.band.a%*%BB.band.inv.a%*%t(B.band.a)

        BB.band.b <- crossprod(as.matrix(B.band.b))
        BB.band.inv.b <- chol2inv(chol(BB.band.b))
        HB.band.b <- B.band.b%*%BB.band.inv.b%*%t(B.band.b)
      }
      if(d.band>1){
        BQ2.band.a <- as.matrix(B.band.a%*%Q2.band.a)
        BB.band.a <- crossprod(BQ2.band.a)
        BB.band.inv.a <- chol2inv(chol(BB.band.a))
        HB.band.a <- BQ2.band.a%*%BB.band.inv.a%*%t(BQ2.band.a)

        BQ2.band.b <- as.matrix(B.band.b%*%Q2.band.b)
        BB.band.b <- crossprod(BQ2.band.b)
        BB.band.inv.b <- chol2inv(chol(BB.band.b))
        HB.band.b <- BQ2.band.b%*%BB.band.inv.b%*%t(BQ2.band.b)
      }
      # Construct SCC
      # Step 2.1. Estimate eta.hat
      eta.hat.a <- R0.a%*%HB.band.a
      eta.hat.b <- R0.b%*%HB.band.b

      # Step 2.2. Estimate sigma.hat
      eps.hat.a <- R0.a-eta.hat.a
      sigma.hat.a <- apply(eps.hat.a^2,2,mean)
      sigma.hat.a <- matrix(sqrt(sigma.hat.a),ncol=1)

      eps.hat.b <- R0.b-eta.hat.b
      sigma.hat.b <- apply(eps.hat.b^2,2,mean)
      sigma.hat.b <- matrix(sqrt(sigma.hat.b),ncol=1)

      # Bootstrap
      # Generate Bootstrap samples
      Ya.boot <- lapply(1:nboot, function(iter){
        tmp1 <- sample(c(-1,1),na,replace=TRUE)
        tmp2 <- sample(c(-1,1),na*npix,replace=TRUE)
        delta1 <- matrix(rep(tmp1,times=npix),nrow=na,ncol=npix)
        delta2 <- matrix(tmp2,nrow=na,ncol=npix)
        return(mu.hat.mtx.a+delta1*eta.hat.a+delta2*eps.hat.a)})

      Yb.boot <- lapply(1:nboot, function(iter){
        tmp1 <- sample(c(-1,1),nb,replace=TRUE)
        tmp2 <- sample(c(-1,1),nb*npix,replace=TRUE)
        delta1 <- matrix(rep(tmp1,times=npix),nrow=nb,ncol=npix)
        delta2 <- matrix(tmp2,nrow=nb,ncol=npix)
        return(mu.hat.mtx.b+delta1*eta.hat.b+delta2*eps.hat.b)})

      Yam.boot <- lapply(1:nboot, function(iter){matrix(apply(Ya.boot[[iter]],2,mean),ncol=1)})
      Ybm.boot <- lapply(1:nboot, function(iter){matrix(apply(Yb.boot[[iter]],2,mean),ncol=1)})

      # Estimate mean function
      mu.hat.boot.a <- lapply(1:nboot, function(iter){
        Slam.a%*%Yam.boot[[iter]]})
      R0.boot.a <- lapply(1:nboot, function(iter){
        mu.boot.mtx.a <- matrix(rep(t(mu.hat.boot.a[[iter]]),times=na),nrow=na,byrow=TRUE)
        R.boot.a <- Ya.boot[[iter]]-mu.boot.mtx.a})

      mu.hat.boot.b <- lapply(1:nboot, function(iter){
        Slam.b%*%Ybm.boot[[iter]]})
      R0.boot.b <- lapply(1:nboot, function(iter){
        mu.boot.mtx.b <- matrix(rep(t(mu.hat.boot.b[[iter]]),times=nb),nrow=nb,byrow=TRUE)
        R.boot.b <- Yb.boot[[iter]]-mu.boot.mtx.b})

      # Estimate Geta and fpca
      eta.boot.a <- lapply(1:nboot, function(iter){
        R0.boot.a[[iter]]%*%HB.band.a})
      Geta.boot.a <- lapply(1:nboot, function(iter){
        crossprod(eta.boot.a[[iter]])/na})
      Geig.boot.a <- lapply(1:nboot, function(iter){
        eigs(Geta.boot.a[[iter]],k=10)})
      evalues.boot.a <- lapply(1:nboot, function(iter){
        Real(Geig.boot.a[[iter]]$values)})
      Kappa.boot.a <- lapply(1:nboot, function(iter){
        sum(cumsum(evalues.boot.a[[iter]])/sum(evalues.boot.a[[iter]])<0.99)+1})
      lambdaK.boot.a <- lapply(1:nboot, function(iter){
        (evalues.boot.a[[iter]])[1:Kappa.boot.a[[iter]]]})
      psi.boot.a <- lapply(1:nboot, function(iter){
        matrix((Geig.boot.a[[iter]]$vectors)[,1:Kappa.boot.a[[iter]]],
               ncol=Kappa.boot.a[[iter]])})
      phi.boot.a <- lapply(1:nboot, function(iter){
        psi.boot.a[[iter]]%*%diag(sqrt(lambdaK.boot.a[[iter]]),
                                  nrow=Kappa.boot.a[[iter]])})
      eps.boot.a <- lapply(1:nboot, function(iter){
        R0.boot.a[[iter]]-eta.boot.a[[iter]]})
      sigma.boot.a <- lapply(1:nboot, function(iter){
        matrix(sqrt(apply(eps.boot.a[[iter]]^2,2,mean)),ncol=1)})

      eta.boot.b <- lapply(1:nboot, function(iter){
        R0.boot.b[[iter]]%*%HB.band.b})
      Geta.boot.b <- lapply(1:nboot, function(iter){
        crossprod(eta.boot.b[[iter]])/nb})
      Geig.boot.b <- lapply(1:nboot, function(iter){
        eigs(Geta.boot.b[[iter]],k=10)})
      evalues.boot.b <- lapply(1:nboot, function(iter){
        Real(Geig.boot.b[[iter]]$values)})
      Kappa.boot.b <- lapply(1:nboot, function(iter){
        sum(cumsum(evalues.boot.b[[iter]])/sum(evalues.boot.b[[iter]])<0.99)+1})
      lambdaK.boot.b <- lapply(1:nboot, function(iter){
        (evalues.boot.b[[iter]])[1:Kappa.boot.b[[iter]]]})
      psi.boot.b <- lapply(1:nboot, function(iter){
        matrix((Geig.boot.b[[iter]]$vectors)[,1:Kappa.boot.b[[iter]]],
               ncol=Kappa.boot.b[[iter]])})
      phi.boot.b <- lapply(1:nboot, function(iter){
        psi.boot.b[[iter]]%*%diag(sqrt(lambdaK.boot.b[[iter]]),
                                  nrow=Kappa.boot.b[[iter]])})
      eps.boot.b <- lapply(1:nboot, function(iter){
        R0.boot.b[[iter]]-eta.boot.b[[iter]]})
      sigma.boot.b <- lapply(1:nboot, function(iter){
        matrix(sqrt(apply(eps.boot.b[[iter]]^2,2,mean)),ncol=1)})

      Veta.boot <- lapply(1:nboot,function(iter){
        Geta.boot.a[[iter]]+rho*Geta.boot.b[[iter]]})
      diag.Veta.boot <- lapply(1:nboot,function(iter){diag(Veta.boot[[iter]])})

      # Generate quantile
      nboot.q <- 1000
      if (adjust.sigma==FALSE){
        Zk.boot.a <- lapply(1:nboot,function(iter){
          matrix(rnorm(Kappa.boot.a[[iter]]*nboot.q,0,1),nboot.q,Kappa.boot.a[[iter]])}) # dim=(1000,Ka)
        Zk.boot.b <- lapply(1:nboot,function(iter){
          matrix(rnorm(Kappa.boot.b[[iter]]*nboot.q,0,1),nboot.q,Kappa.boot.b[[iter]])}) # dim=(1000,Ka)

        tmp1.boot <- lapply(1:nboot,function(iter){diag(1/sqrt(diag.Veta.boot[[iter]]))})
        tmp2.boot <- lapply(1:nboot,function(iter){
          phi.boot.a[[iter]]%*%t(Zk.boot.a[[iter]])-sqrt(rho)*phi.boot.b[[iter]]%*%t(Zk.boot.b[[iter]])})
        zeta.boot <- lapply(1:nboot,function(iter){tmp1.boot[[iter]]%*%tmp2.boot[[iter]]})
        Qalpha.boot <- lapply(1:nboot,function(iter){quantile(apply(abs(zeta.boot[[iter]]),2,max),c(1-alpha0))})
        Sigma.boot.a <- NA
        Sigma.boot.b <- NA
        Veta.boot <- NA
      }else if (adjust.sigma==TRUE){
        Sig.boot.a <- lapply(1:nboot,function(iter){
          Geta.boot.a[[iter]]+Slam.a%*%diag(as.numeric(sigma.boot.a[[iter]])^2)%*%t(Slam.a)})

        Sig.boot.b <- lapply(1:nboot,function(iter){
          Geta.boot.b[[iter]]+Slam.b%*%diag(as.numeric(sigma.boot.b[[iter]])^2)%*%t(Slam.b)})

        Veta.sigma.boot <- lapply(1:nboot,function(iter){Sig.boot.a[[iter]]+rho*Sig.boot.b[[iter]]})
        diag.Veta.sigma.boot <- lapply(1:nboot,function(iter){diag(Veta.sigma.boot[[iter]])})
        tmp1.boot <- lapply(1:nboot,function(iter){diag(1/sqrt(diag.Veta.sigma.boot[[iter]]))})

        Zk.boot.a <- lapply(1:nboot,function(iter){
          matrix(rnorm(Kappa.boot.a[[iter]]*nboot.q,0,1),nboot.q,Kappa.boot.a[[iter]])}) # dim=(1000,Kappa.a)
        Zk.boot.b <- lapply(1:nboot,function(iter){
          matrix(rnorm(Kappa.boot.b[[iter]]*nboot.q,0,1),nboot.q,Kappa.boot.b[[iter]])}) # dim=(1000,Kappa.b)

        Zeps.boot.a <- lapply(1:nboot,function(iter){matrix(rnorm(npix*nboot.q,0,1),nboot.q,npix)})
        Zeps.boot.b <- lapply(1:nboot,function(iter){matrix(rnorm(npix*nboot.q,0,1),nboot.q,npix)})

        tmp2.boot.a <- lapply(1:nboot,function(iter){
          phi.boot.a[[iter]]%*%t(Zk.boot.a[[iter]])+Slam.a%*%t(Zeps.boot.a[[iter]]%*%diag(as.numeric(sigma.boot.a[[iter]])))})
        tmp2.boot.b <- lapply(1:nboot,function(iter){
          phi.boot.b[[iter]]%*%t(Zk.boot.b[[iter]])+Slam.b%*%t(Zeps.boot.b[[iter]]%*%diag(as.numeric(sigma.boot.b[[iter]])))})

        zeta.boot <- lapply(1:nboot,function(iter){
          tmp1.boot[[iter]]%*%(tmp2.boot.a[[iter]]-sqrt(rho)*tmp2.boot.b[[iter]])})
        Qalpha.boot <- lapply(1:nboot,function(iter){
          quantile(apply(abs(zeta.boot[[iter]]),2,max),probs=c(1-alpha0))})
      }
      if(adjust.sigma==FALSE){
        ME.boot <- lapply(1:nboot, function(iter){
          matrix(rep(Qalpha.boot[[iter]],each=npix),nrow=npix,ncol=nalpha)*
            matrix(rep(sqrt(diag.Veta.boot[[iter]])/sqrt(na),times=nalpha),nrow=npix,ncol=nalpha)})

      }else if(adjust.sigma==TRUE){
        ME.boot <- lapply(1:nboot, function(iter){
          matrix(rep(Qalpha.boot[[iter]],each=npix),nrow=npix,ncol=nalpha)*
            matrix(rep(sqrt(diag.Veta.sigma.boot[[iter]])/sqrt(na),times=nalpha),nrow=npix,ncol=nalpha)})
      }
      mudiff.boot.mtx <- lapply(1:nboot, function(iter){
        matrix(rep(mu.hat.boot.b[[iter]]-mu.hat.boot.a[[iter]],times=nalpha),nrow=npix,ncol=nalpha)})
      scc.l.boot <- lapply(1:nboot, function(iter){
        mudiff.boot.mtx[[iter]]-ME.boot[[iter]]})
      scc.u.boot <- lapply(1:nboot, function(iter){
        mudiff.boot.mtx[[iter]]+ME.boot[[iter]]})

      mudiff.nalpha.mtx <- matrix(rep(mu.hat.b-mu.hat.a,times=nalpha),nrow=npix,ncol=nalpha)
      cover.ib <- lapply(1:nboot, function(iter){
        apply(((mudiff.nalpha.mtx<=scc.u.boot[[iter]]) &
                 (mudiff.nalpha.mtx>=scc.l.boot[[iter]])),2,mean)})
      cover.tmp <- do.call(rbind,cover.ib)
      cover.tri <- rbind(cover.tri,cover.tmp)
    }
    cover.all[[itri]] <- cover.tri
  }
  # Results
  cover.avg <- do.call(cbind,(lapply(cover.all,function(x) apply(x==1,2,mean))))
  cover.diff <- apply(cover.avg,2,'-',(1-alpha.grid))
  tri.band <- tri.boot[which.min(apply(cover.diff^2,2,sum)),]
  return(tri.band)
}
