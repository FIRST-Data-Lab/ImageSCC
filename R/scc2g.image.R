#' Simultaneous confidence corridors for mean function of imaging data: two sample case
#' @importFrom BPST basis
#' @importFrom RSpectra eigs
#' @importFrom pracma Real
#'
scc2g.image <- function(Ya,Yb,Z,Z.band,d.est,d.band,r,
                        B.est.a,Q2.est.a,K.est.a,B.est.b,Q2.est.b,K.est.b,
                        B.band.a,Q2.band.a,K.band.a,B.band.b,Q2.band.b,K.band.b,ind.inside,
                        B.cover1.a,B.cover2.a,B.cover1.b,B.cover2.b,ind.inside.cover,
                        penalty=TRUE,lambda,alpha.grid=c(0.1,0.05,0.01),adjust.sigma=TRUE){
  npix <- length(ind.inside)
  npix.cover <- length(ind.inside.cover)
  nalpha <- length(alpha.grid)
  na <- nrow(Ya); nb <- nrow(Yb); rho <- na/nb;

  Ya <- Ya[,ind.inside]; Yb <- Yb[,ind.inside];
  Yam <- matrix(apply(Ya,2,mean),nrow=1)
  Ybm <- matrix(apply(Yb,2,mean),nrow=1)

  # Step 1. Estimate Mean Function for Each Group
  mfit0.a <- fit.mean(B.est.a,Q2.est.a,K.est.a,lambda,Yam,proj.matrix=TRUE)
  gamma.hat0.a <- mfit0.a$gamma
  mu.hat.a <- mfit0.a$Yhat
  mu.hat.mtx.a <- matrix(rep(drop(mu.hat.a),times=na),nrow=na,byrow=TRUE)

  mfit0.b <- fit.mean(B.est.b,Q2.est.b,K.est.b,lambda,Ybm,proj.matrix=TRUE)
  gamma.hat0.b <- mfit0.b$gamma
  mu.hat.b <- mfit0.b$Yhat
  mu.hat.mtx.b <- matrix(rep(drop(mu.hat.b),times=nb),nrow=nb,byrow=TRUE)

  mu.hat.cover.a <- B.cover1.a%*%gamma.hat0.a
  mu.hat.cover.b <- B.cover1.b%*%gamma.hat0.b

  R0.a <- Ya-mu.hat.mtx.a
  R0.b <- Yb-mu.hat.mtx.b

  # Step 2.1. Estimate eta.hat & Geta.hat;
  if(d.band==-1){# using piecewise constant
    BB.band.a <- crossprod(B.band.a)
    BB.band.a.inv <- chol2inv(chol(BB.band.a))
    BR.band.a <- crossprod(B.band.a,t(R0.a))
    mfit1.a <- B.band.a%*%BB.band.a.inv%*%BR.band.a
    eta.hat.a <- t(mfit1.a)

    eta.cover.a <- B.cover2.a%*%BB.band.a.inv%*%BR.band.a
    eta.cover.a <- t(eta.cover.a)

    BB.band.b <- crossprod(B.band.b)
    BB.band.b.inv <- chol2inv(chol(BB.band.b))
    BR.band.b <- crossprod(B.band.b,t(R0.b))
    mfit1.b <- B.band.b%*%BB.band.b.inv%*%BR.band.b
    eta.hat.b <- t(mfit1.b)

    eta.cover.b <- B.cover2.b%*%BB.band.b.inv%*%BR.band.b
    eta.cover.b <- t(eta.cover.b)
  }
  if(d.band>1){# using higher order spline
    if(penalty==FALSE){
      BQ2.band.a <- as.matrix(B.band.a%*%Q2.band.a)
      BB.band.a <- crossprod(BQ2.band.a)
      BB.band.a.inv <- chol2inv(chol(BB.band.a))
      BR.band.a <- crossprod(BQ2.band.a,t(R0.a))

      BQ2.band.b <- as.matrix(B.band.b%*%Q2.band.b)
      BB.band.b <- crossprod(BQ2.band.b)
      BB.band.b.inv <- chol2inv(chol(BB.band.b))
      BR.band.b <- crossprod(BQ2.band.b,t(R0.b))

      mfit1.a <- BQ2.band.a%*%BB.band.a%*%BR.band.a
      eta.hat.a <- t(mfit1.a)
      mfit1.b <- BQ2.band.b%*%BB.band.b%*%BR.band.b
      eta.hat.b <- t(mfit1.b)

      BQ2.cover2.a <- as.matrix(B.cover2.a%*%Q2.band.a)
      eta.cover.a <- BQ2.cover2.a%*%BB.band.a%*%BR.band.a
      eta.cover.a <- t(eta.cover.a)

      BQ2.cover2.b <- as.matrix(B.cover2.b%*%Q2.band.b)
      eta.cover.b <- BQ2.cover2.b%*%BB.band.b%*%BR.band.b
      eta.cover.b <- t(eta.cover.b)
    }
    if(penalty==TRUE){
      mfit1.a <- fit.mean(B.band.a,Q2.band.a,K.band.a,lambda,R0.a)
      eta.hat.a <- mfit1.a$Yhat
      gamma.hat.a <- mfit1.a$gamma  # dim #basis * #images   gamma=Q2%*%theta
      Beta.hat.a <- crossprod(gamma.hat.a)/na

      eta.cover.a <- as.matrix(B.cover2.a%*%gamma.hat.a)
      eta.cover.a <- t(eta.cover.a)

      mfit1.b <- fit.mean(B.band.b,Q2.band.b,K.band.b,lambda,R0.b)
      eta.hat.b <- mfit1.b$Yhat
      gamma.hat.b <- mfit1.b$gamma  # dim #basis * #images   gamma=Q2%*%theta
      Beta.hat.b <- crossprod(gamma.hat.b)/nb

      eta.cover.b <- as.matrix(B.cover2.b%*%gamma.hat.b)
      eta.cover.b <- t(eta.cover.b)
    }
  }

  Geta.hat.a <- crossprod(eta.hat.a)/na
  diag.Geta.hat.a <- diag(Geta.hat.a)
  diag.Geta.cover.a <- apply(eta.cover.a^2,2,mean)

  Geta.hat.b <- crossprod(eta.hat.b)/nb
  diag.Geta.hat.b <- diag(Geta.hat.b)
  diag.Geta.cover.b <- apply(eta.cover.b^2,2,mean)

  Veta.hat <- Geta.hat.a+rho*Geta.hat.b
  diag.Veta.hat <- diag.Geta.hat.a+rho*diag.Geta.hat.b
  diag.Veta.cover <- diag.Geta.cover.a+rho*diag.Geta.cover.b

  # Step 2.2. Estimate sigma.hat
  eps.hat.a <- R0.a-eta.hat.a
  sigma.hat.a <- apply(eps.hat.a^2,2,mean)
  sigma.hat.a <- matrix(sqrt(sigma.hat.a),ncol=1)

  eps.hat.b <- R0.b-eta.hat.b
  sigma.hat.b <- apply(eps.hat.b^2,2,mean)
  sigma.hat.b <- matrix(sqrt(sigma.hat.b),ncol=1)

  # Step 3. Get lambda.hat & phi.hat
  Geig.a <- eigs(Geta.hat.a,k=10)
  evalues.a <- Real(Geig.a$values)
  Kappa.a <- sum(cumsum(evalues.a)/sum(evalues.a)<0.99)+1
  lambdaK.hat.a <- evalues.a[1:Kappa.a]
  psi.hat.a <- matrix((Geig.a$vectors)[,1:Kappa.a],ncol=Kappa.a) #eigenfunction
  phi.hat.a <- psi.hat.a%*%diag(sqrt(lambdaK.hat.a),nrow=Kappa.a) #sqrt{lambdaK}*eigenfunction

  Geig.b <- eigs(Geta.hat.b,k=10)
  evalues.b <- Real(Geig.b$values)
  Kappa.b <- sum(cumsum(evalues.b)/sum(evalues.b)<0.99)+1
  lambdaK.hat.b <- evalues.b[1:Kappa.b]
  psi.hat.b <- matrix((Geig.b$vectors)[,1:Kappa.b],ncol=Kappa.b) #eigenfunction
  phi.hat.b <- psi.hat.b%*%diag(sqrt(lambdaK.hat.b),nrow=Kappa.b) #sqrt{lambdaK}*eigenfunction

  # Step 4. Generate qunatile
  nboot <- 1000
  if (adjust.sigma==FALSE){
    Zk.a <- matrix(rnorm(Ka*nboot,0,1),nboot,Kappa.a) # dim=(1000,Ka)
    Zk.b <- matrix(rnorm(Kb*nboot,0,1),nboot,Kappa.b) # dim=(1000,Kb)
    tmp.hat <- diag(1/sqrt(diag.Veta.hat))
    zeta.hat <- tmp.hat%*%(phi.hat.a%*%t(Zk.a)-sqrt(rho)*phi.hat.b%*%t(Zk.b))
    Qalpha.hat <- quantile(apply(abs(zeta.hat),2,max),c(1-alpha.grid))
    Sigma.hat.a <- NA
    Sigma.hat.b <- NA
    Veta.sigma <- NA
  }else if (adjust.sigma==TRUE){
    Proj.est.a <- mfit0.a$Slam
    Sigma.hat.a <- Geta.hat.a+
      Proj.est.a%*%diag(as.numeric(sigma.hat.a)^2)%*%t(Proj.est.a)

    Proj.est.b <- mfit0.b$Slam
    Sigma.hat.b <- Geta.hat.b+
      Proj.est.b%*%diag(as.numeric(sigma.hat.b)^2)%*%t(Proj.est.b)

    Veta.sigma <- Sigma.hat.a+rho*Sigma.hat.b
    diag.Veta.sigma <- diag(Veta.sigma)
    tmp.hat <- diag(1/sqrt(diag.Veta.sigma))

    Zk.a <- matrix(rnorm(Kappa.a*nboot,0,1),nboot,Kappa.a)
    Zk.b <- matrix(rnorm(Kappa.b*nboot,0,1),nboot,Kappa.b)

    Zeps.a <- matrix(rnorm(npix*nboot,0,1),nboot,npix)
    Zeps.b <- matrix(rnorm(npix*nboot,0,1),nboot,npix)

    zeta.hat.a <- phi.hat.a%*%t(Zk.a)+Proj.est.a%*%t(Zeps.a%*%diag(as.numeric(sigma.hat.a)))
    zeta.hat.b <- phi.hat.b%*%t(Zk.b)+Proj.est.b%*%t(Zeps.b%*%diag(as.numeric(sigma.hat.b)))
    zeta.hat <- tmp.hat%*%(zeta.hat.a-sqrt(rho)*zeta.hat.b)
    Qalpha.hat <- quantile(apply(abs(zeta.hat),2,max),probs=c(1-alpha.grid))
  }

  # Step 5. Generate SCC
  if(adjust.sigma==FALSE){
    ME <- matrix(rep(Qalpha.hat,each=npix.cover),nrow=npix.cover,ncol=nalpha)*
      matrix(rep(sqrt(diag.Veta.cover)/sqrt(na),times=nalpha),nrow=npix.cover,ncol=nalpha)
    mu.hat.diff <- matrix(rep(mu.hat.cover.b-mu.hat.cover.a,times=nalpha),nrow=npix.cover,ncol=nalpha)
    scc.l <- mu.hat.diff-ME
    scc.u <- mu.hat.diff+ME
    bw <- 2*apply(ME,2,mean)
  }else if(adjust.sigma==TRUE){
    match.band <- match(paste(round(Z.band[ind.inside.cover,1],6),round(Z.band[ind.inside.cover,2],6)),
                     paste(round(Z[ind.inside,1],6),round(Z[ind.inside,2],6)))
    diag.Veta_sigma.cover <- diag.Veta.sigma[match.band]
    ME <- matrix(rep(Qalpha.hat,each=npix.cover),nrow=npix.cover,ncol=nalpha)*
      matrix(rep(sqrt(diag.Veta_sigma.cover)/sqrt(na),times=nalpha),nrow=npix.cover,ncol=nalpha)
    mu.hat.diff <- matrix(rep(mu.hat.cover.b-mu.hat.cover.a,times=nalpha),nrow=npix.cover,ncol=nalpha)
    scc.l <- mu.hat.diff-ME
    scc.u <- mu.hat.diff+ME
    bw <- 2*apply(ME,2,mean)
  }

  scc <- array(NA,dim=c(npix.cover,2,nalpha))
  cover.zero <- matrix(NA,nrow=npix.cover,ncol=nalpha)
  for(j in 1:nalpha){
    scc[,,j] <- cbind(scc.l[,j],scc.u[,j])
    cover.zero[,j] <- matrix(as.numeric((0>scc.u[,j]))-as.numeric((0<scc.l[,j])))
  }
  list(na,nb,mfit.a=mfit0.a,mfit.b=mfit0.b,
       mu.hat.a=mu.hat.a,mu.hat.b=mu.hat.b,
       mu.hat.cover.a=mu.hat.cover.a,mu.hat.cover.b=mu.hat.cover.b,
       eta.hat.a=eta.hat.a,eta.hat.b=eta.hat.b,
       eta.cover.a=eta.cover.a,eta.cover.b=eta.cover.b,
       Geta.hat.a=Geta.hat.a,Geta.hat.b=Geta.hat.b,
       Sigma.hat.a=Sigma.hat.a,Sigma.hat.b=Sigma.hat.b,
       scc=scc,cover.zero=cover.zero,bw=bw,
       ind.inside=ind.inside,ind.inside.cover=ind.inside.cover)
}

