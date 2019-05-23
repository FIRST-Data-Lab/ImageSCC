#' Select triangulation for one sample SCC construction using bootstrap
#'
boot1g.image<-function(Y,Z,d.band,r,B.est,Q2.est,K.est,
                     B.bands,Q2.bands,K.bands,ind.inside,lambda,nboot=50,
                     alpha0,adjust.sigma=TRUE){
  Ym <- matrix(apply(Y,2,mean),nrow=1)
  n <- nrow(Y)
  npix <- length(ind.inside)
  nalpha <- 10
  alpha.grid <- seq((alpha0-0.005),(alpha0+0.005),length.out=nalpha)
  ntri <- length(B.bands)

  # Step 1. Estimate mu.hat
  mfit0 <- fit.mean(B.est,Q2.est,K.est,lambda,Ym,proj.matrix=TRUE)
  lamc <- mfit0$lambdac
  theta <- mfit0$theta
  gamma <- mfit0$gamma
  mu.hat <- mfit0$Yhat
  Slam <- mfit0$Slam
  mu.hat.mtx <- matrix(rep(drop(mu.hat),times=n),nrow=n,byrow=TRUE)
  R0 <- Y-mu.hat.mtx

  cover.all <- vector('list',length=ntri)
  for(itri in 1:ntri){
    B.band <- B.bands[[itri]]
    Q2.band <- Q2.bands[[itri]]
    K.band <- K.bands[[itri]]

    if(d.band==-1){
      BB.band <- crossprod(as.matrix(B.band))
      BB.band.inv <- chol2inv(chol(BB.band))
      HB.band <- B.band%*%BB.band.inv%*%t(B.band)
    }
    if(d.band>1){
      BQ2.band <- as.matrix(B.band%*%Q2.band)
      BB.band <- crossprod(BQ2.band)
      BB.band.inv <- chol2inv(chol(BB.band))
      HB.band <- BQ2.band%*%BB.band.inv%*%t(BQ2.band)
    }

    # Construct SCC
    # Step 2.1. Estimate eta.hat & Geta.hat;
    eta.hat <- R0%*%HB.band
    Geta.hat <- crossprod(eta.hat)/n

    # Step 2.2. Estimate sigma.hat
    eps.hat <- R0-eta.hat
    sigma.hat <- apply(eps.hat^2,2,mean)
    sigma.hat <- matrix(sqrt(sigma.hat),ncol=1)

    # Step 3. Bootstrap
    # Step 3.1. Generate bootstrap sample;
    Y.boot <- lapply(1:nboot, function(iter){
      tmp1 <- sample(c(-1,1),n,replace=TRUE)
      tmp2 <- sample(c(-1,1),n*npix,replace=TRUE)
      delta1 <- matrix(rep(tmp1,times=npix),nrow=n,ncol=npix)
      delta2 <- matrix(tmp2,nrow=n,ncol=npix)
      return(mu.hat.mtx+delta1*eta.hat+delta2*eps.hat)})
    Ym.boot <- lapply(1:nboot, function(iter){matrix(apply(Y.boot[[iter]],2,mean),ncol=1)})
    # Step 3.2. Estimate bootstrap sample;
    mu.hat.boot <- lapply(1:nboot, function(iter){
      Slam%*%Ym.boot[[iter]]})
    R0.boot <- lapply(1:nboot, function(iter){
      mu.boot.mtx <- matrix(rep(t(mu.hat.boot[[iter]]),times=n),nrow=n,byrow=TRUE)
      R.boot <- Y.boot[[iter]]-mu.boot.mtx})
    eta.boot <- lapply(1:nboot, function(iter){
      R0.boot[[iter]]%*%HB.band})
    Geta.boot <- lapply(1:nboot, function(iter){
      crossprod(eta.boot[[iter]])/n})
    Geig.boot <- lapply(1:nboot, function(iter){
      eigs(Geta.boot[[iter]],k=10)})
    evalues.boot <- lapply(1:nboot, function(iter){
      Real(Geig.boot[[iter]]$values)})
    Kappa.boot <- lapply(1:nboot, function(iter){
      sum(cumsum(evalues.boot[[iter]])/sum(evalues.boot[[iter]])<0.99)+1})
    lambdaK.boot <- lapply(1:nboot, function(iter){
      (evalues.boot[[iter]])[1:Kappa.boot[[iter]]]})
    psi.boot <- lapply(1:nboot, function(iter){
      matrix((Geig.boot[[iter]]$vectors)[,1:Kappa.boot[[iter]]],
             ncol=Kappa.boot[[iter]])})
    phi.boot <- lapply(1:nboot, function(iter){
      psi.boot[[iter]]%*%diag(sqrt(lambdaK.boot[[iter]]),
                              nrow=Kappa.boot[[iter]])})
    eps.boot <- lapply(1:nboot, function(iter){
      R0.boot[[iter]]-eta.boot[[iter]]})
    sigma.boot <- lapply(1:nboot, function(iter){
      matrix(sqrt(apply(eps.boot[[iter]]^2,2,mean)),ncol=1)})
    # Step 3.3. Construct bootstrap SCC
    nb <- 1000
    if(adjust.sigma==FALSE){
      Zb <- lapply(1:nboot, function(iter){
        matrix(rnorm(Kappa.boot[[iter]]*nb,0,1),nb,Kappa.boot[[iter]])})
      tmp1.boot <- lapply(1:nboot, function(iter){diag(1/sqrt(diag(Geta.boot[[iter]])))})
      tmp2.boot <- lapply(1:nboot, function(iter){phi.boot[[iter]]%*%t(Zb[[iter]])})
      zeta.boot <- lapply(1:nboot, function(iter){tmp1.boot[[iter]]%*%tmp2.boot[[iter]]})
      Qalpha.boot <- lapply(1:nboot, function(iter){
        quantile(apply(abs(zeta.boot[[iter]]),2,max),c(1-alpha.grid))})
      Sig.boot <- NA
    }else if(adjust.sigma==TRUE){
      Zb <- lapply(1:nboot, function(iter){
        matrix(rnorm(Kappa.boot[[iter]]*nb,0,1),nb,Kappa.boot[[iter]])})
      Zeps <- lapply(1:nboot, function(iter){
        matrix(rnorm(npix*nb,0,1),nb,npix)})
      Sig.boot <- lapply(1:nboot, function(iter){
        Geta.boot[[iter]]+Slam%*%
          diag(as.numeric(sigma.boot[[iter]])^2)%*%t(Slam)})
      tmp1.boot <- lapply(1:nboot, function(iter){
        diag(1/sqrt(diag(Sig.boot[[iter]])))})
      tmp2.boot <- lapply(1:nboot, function(iter){
        phi.boot[[iter]]%*%t(Zb[[iter]])+
          Slam%*%t(Zeps[[iter]]%*%diag(as.numeric(sigma.boot[[iter]])))})
      zeta.boot <- lapply(1:nboot, function(iter){tmp1.boot[[iter]]%*%tmp2.boot[[iter]]})
      Qalpha.boot <- lapply(1:nboot, function(iter){
        quantile(apply(abs(zeta.boot[[iter]]),2,max),c(1-alpha.grid))})
    }

    if(adjust.sigma==FALSE){
      ME.boot <- lapply(1:nboot, function(iter){
        matrix(rep(Qalpha.boot[[iter]],each=npix),nrow=npix,ncol=nalpha)*
          matrix(rep(sqrt(diag(Geta.boot[[iter]]))/sqrt(n),times=nalpha),nrow=npix,ncol=nalpha)})
    }else if(adjust.sigma==TRUE){
      ME.boot <- lapply(1:nboot, function(iter){
        matrix(rep(Qalpha.boot[[iter]],each=npix),nrow=npix,ncol=nalpha)*
          matrix(rep(sqrt(diag(Sig.boot[[iter]]))/sqrt(n),times=nalpha),nrow=npix,ncol=nalpha)})
    }

    mu.boot.mtx <- lapply(1:nboot, function(iter){
      matrix(rep(mu.hat.boot[[iter]],times=nalpha),nrow=npix,ncol=nalpha)})
    scc.l.boot <- lapply(1:nboot, function(iter){
      mu.boot.mtx[[iter]]-ME.boot[[iter]]})
    scc.u.boot <- lapply(1:nboot, function(iter){
      mu.boot.mtx[[iter]]+ME.boot[[iter]]})

    mu.nalpha.mtx <- matrix(rep(mu.hat,times=nalpha),nrow=npix,ncol=nalpha)
    cover.ib <- lapply(1:nboot, function(iter){
      apply(((mu.nalpha.mtx<=scc.u.boot[[iter]]) &
               (mu.nalpha.mtx>=scc.l.boot[[iter]])),2,mean)})
    cover.tmp <- do.call(rbind,cover.ib)
    cover.all[[itri]] <- cover.tmp
  }
  # Results
  cover.avg <- do.call(cbind,(lapply(cover.all,function(x) apply(x==1,2,mean))))
  cover.diff <- apply(cover.avg,2,'-',(1-alpha.grid))
  tri.band <- (1:ntri)[which.min(apply(cover.diff^2,2,sum))]
  return(tri.band)
}
