#' Generate simulation examples: one sample case
#'
#' The function generates one-sample simulation data using four different mean functions
#' and scheme described in the paper.
#'
#' @param n number of images in the simulated data.
#' @param Z a 2-column matrix specifying locations of information.
#' @param ind.inside a vector of indices specifying which locations in matrix Z fall inside
#' the irregular domain.
#' @param mu.func a integer in \code{c(1,2,3,4)} that specifies what mean function is used for
#' data generation. Value \code{1,2,3,4} correspond to quadratic, exponential, cubic and sine
#' functions described in the paper respectively.
#' @param noise.type 'Func'/'Const' indicating heterogeneous or homogeneous variance for measurement error.
#' @param lam1,lam2 the eigenvalues used to adjust subject-level variation in simulation data.
#' @param iter random seed.
#'
#' @details  This R package is the implementation program for manuscript entitled "Simultaneous Confidence Corridors for Mean Functions in Functional Data Analysis of Imaging Data" by Yueying Wang, Guannan Wang, Li Wang and R. Todd Ogden.
#'
#' @export
data1g.image <- function(n,Z,ind.inside,mu.func=c(1,2,3,4),noise.type=c('Func','Const'),
                         lam1,lam2,iter=2019){
  xx <- Z[,1]
  uu <- unique(xx)
  n1 <- length(uu)
  yy <- Z[,2]
  vv <- unique(yy)
  n2 <- length(vv)
  npix <- n1*n2;

  # Generate true mean function
  if(mu.func==1){
    beta.true <- 20*((xx-0.5)^2+(yy-0.5)^2)
  }else if(mu.func==2){
    beta.true <- 5*exp(-15*((xx-0.5)^2+(yy-0.5)^2))+0.5
  }else if(mu.func==3){
    beta.true <- 3.2*(-xx^3+yy^3)+2.4
  }else if(mu.func==4){
    # beta.true <- -10*sin(7*pi*(xx+0.09))+10*sin(7*pi*(yy-0.14))+2.5
    beta.true <- -10*sin(5*pi*(xx+0.22))+10*sin(5*pi*(yy-0.18))+2.8
  }
  ind.outside <- setdiff(1:npix,ind.inside)
  beta.true[ind.outside] <- NA
  beta.true <- as.matrix(beta.true)

  set.seed(iter)
  psi1 <- 0.9878812*sin(pi*xx)+0.5
  psi1.mtx <- matrix(psi1,ncol=n1,nrow=n2)
  psi2 <- 2.157205*(cos(pi*yy)-0.03903665)
  psi2.mtx <- matrix(psi2,ncol=n1,nrow=n2)
  psi.true <- cbind(psi1,psi2)
  phi.true <- psi.true%*%diag(sqrt(c(lam1,lam2)),nrow=2)
  phi1.mtx <- matrix(phi.true[,1],ncol=n1,nrow=n2)
  phi2.mtx <- matrix(phi.true[,2],ncol=n1,nrow=n2)
  phi.mtx <- phi1.mtx+phi2.mtx

  Y <- c()
  signal <- c()
  noise <- c()
  eta <- c()
  eps1 <- rnorm(n,mean=0,sd=sqrt(lam1))
  eps2 <- rnorm(n,mean=0,sd=sqrt(lam2))
  eps <- c();
  if(noise.type=='Func'){
    sigma.func <- 0.25*(1-(xx-0.5)^2-(yy-0.5)^2)
  }else{
    sigma.func <- 0.2
  }
  for(i in 1:n){
    epsi <- rnorm(n1*n2,mean=0,sd=sigma.func)
    eps <- rbind(eps,epsi)
    signali <- t(beta.true)
    etai <- eps1[i]*psi1+eps2[i]*psi2
    noisei <- etai+epsi
    noisei <- matrix(noisei,nrow=1)
    Yi <- signali+noisei
    signal <- rbind(signal,signali)
    eta <- rbind(eta,etai)
    noise <- rbind(noise,noisei)
    Y <- rbind(Y,Yi)
  }
  Y[,ind.outside] <- NA
  phi.true[ind.outside,] <- NA
  signal[,ind.outside] <- NA
  noise[ind.outside] <- NA
  list(Y=Y,Z=Z,beta.true=beta.true,psi.true=psi.true,phi.true=phi.true,eta.true=eta,
       signal=signal,noise=noise,eps=eps)
}
