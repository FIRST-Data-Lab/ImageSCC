#' Generate simulation examples: two sample case
#'
#' The function generates two-sample simulation data using four different mean functions
#' and eigenfunctions described in the paper.
#'
#' @param na number of images/subjects in the first sample.
#' @param nb number of images/subjects in the second sample.
#' @param Z a 2-column matrix specifying locations of information.
#' @param ind.inside a vector of indices specifying which locations in matrix Z fall inside
#' the irregular domain.
#' @param mu1.func a integer in \code{c(1,2,3,4)} that specifies what mean function is used for
#' data generation of the first sample. Value \code{1,2,3,4} correspond to quadratic, exponential, cubic and sine
#' functions described in the paper respectively.
#' @param noise.type 'Func'/'Const' indicating heterogeneous or homogeneous variance for measurement error.
#' @param lam1,lam2 the eigenvalues used to adjust subject-level variation in simulation data.
#' @param delta a parameter indicating the scale of difference between two samples' mean functions.
#' @param iter random seed.
#'
#' @details  This R package is the implementation program for manuscript entitled "Simultaneous Confidence Corridors for Mean Functions in Functional Data Analysis of Imaging Data" by Yueying Wang, Guannan Wang, Li Wang and R. Todd Ogden.
#'
#' @export
data2g.image <- function(na,nb,Z,ind.inside,mu1.func,noise.type=c('Func','Const'),
                         lam1,lam2,delta,iter=2019){
  xx <- Z[,1]
  uu <- unique(xx)
  n1 <- length(uu)
  yy <- Z[,2]
  vv <- unique(yy)
  n2 <- length(vv)
  npix <- n1*n2;

  # Generate true mean function
  if(mu1.func==1){
    beta1.true <- 20*((xx-0.5)^2+(yy-0.5)^2)
  }else if(mu1.func==2){
    beta1.true <- 5*exp(-15*((xx-0.5)^2+(yy-0.5)^2))+0.5
  }else if(mu1.func==3){
    beta1.true <- 3.2*(-xx^3+yy^3)+2.4
  }else if(mu1.func==4){
    # beta.true <- -10*sin(7*pi*(xx+0.09))+10*sin(7*pi*(yy-0.14))+2.5
    beta1.true <- -10*sin(5*pi*(xx+0.22))+10*sin(5*pi*(yy-0.18))+2.8
  }

  ind.outside <- setdiff(1:npix,ind.inside)
  beta1.true[ind.outside] <- NA
  beta1.mtx <- matrix(beta1.true,ncol=n1,nrow=n2)

  beta2.true <- beta1.true+delta*(-xx^3+yy^3)
  beta2.true[ind.outside] <- NA
  beta2.mtx <- matrix(beta2.true,ncol=n1,nrow=n2)
  beta.true <- cbind(beta1.true,beta2.true)

  set.seed(iter)
  psi1 <- 0.9878812*sin(pi*xx)+0.5
  psi1.mtx <- matrix(psi1,ncol=n1,nrow=n2)
  psi2 <- 2.157205*(cos(pi*yy)-0.03903665)
  psi2.mtx <- matrix(psi2,ncol=n1,nrow=n2)
  psi.true <- cbind(psi1,psi2)
  phi.true <- psi.true%*%diag(sqrt(c(lam1,lam2)),nrow=2)
  phi.true[ind.inside,] <- NA

  Y1 <- Y2 <- c()
  signal1 <- signal2 <- c()
  noise1 <- noise2 <- c()
  eta1 <- eta2 <- c()
  eps1 <- list(rnorm(na,mean=0,sd=sqrt(lam1)),
            rnorm(nb,mean=0,sd=sqrt(lam1)))
  eps2 <- list(rnorm(na,mean=0,sd=sqrt(lam2)),
            rnorm(nb,mean=0,sd=sqrt(lam2)))
  if(noise.type=='Func'){
    sigma.func <- 0.25*(1-(xx-0.5)^2-(yy-0.5)^2)
  }else{
    sigma.func <- 0.1
  }
  for (i in 1:na){
    eps1i <- rnorm(n1*n2,mean=0,sd=sigma.func)
    signal1i <- t(beta.true[,1])
    eta1i <- eps1[[1]][i]*psi1+eps2[[1]][i]*psi2
    noise1i <- eta1i+eps1i
    Y1i <- signal1i+noise1i
    signal1 <- rbind(signal1,signal1i)
    eta1 <- rbind(eta1,eta1i)
    noise1 <- rbind(noise1,noise1i)
    Y1 <- rbind(Y1,Y1i)
  }
  for (i in 1:nb){
    eps2i <- rnorm(n1*n2,mean=0,sd=sigma.func)
    signal2i <- t(beta.true[,2])
    eta2i <- eps1[[2]][i]*psi1+eps2[[2]][i]*psi2
    noise2i <- eta2i+eps2i
    Y2i <- signal2i+noise2i
    signal2 <- rbind(signal2,signal2i)
    eta2 <- rbind(eta2,eta2i)
    noise2 <- rbind(noise2,noise2i)
    Y2 <- rbind(Y2,Y2i)
  }
  Y1[,ind.outside] <- NA
  Y2[,ind.outside] <- NA
  signal1[,ind.outside] <- NA
  signal2[,ind.outside] <- NA
  phi.true[ind.outside,] <- NA
  noise1[,ind.outside] <- NA
  noise2[,ind.outside] <- NA
  list(Ya=Y1,Yb=Y2,Z=Z,beta.true=beta.true,psi.true=psi.true,phi.true=phi.true,
       signal1=signal1,signal2=signal2,eta1=eta1,eta2=eta2,delta=delta,
       noise1=noise1,noise2=noise2)
}
