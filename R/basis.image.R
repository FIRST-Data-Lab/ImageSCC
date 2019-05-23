#' Basis genration of all the triangulation
#'
#' @export
basis.image <- function(Z,V.ests,Tr.ests,d.est,r){
  ntr <- length(Tr.ests)
  B.ests <- vector('list',length=ntr)
  Q2.ests <- vector('list',length=ntr)
  K.ests <- vector('list',length=ntr)
  ind.inside <- 1:nrow(Z)

  for(itr in 1:ntr){
    V.est <- V.ests[[itr]]
    Tr.est <- Tr.ests[[itr]]
    Bfull.est <- basis(V.est,Tr.est,d.est,r,Z)
    B.est <- Bfull.est$B
    ind.inside.est <- Bfull.est$Ind.inside
    Q2.est <- Bfull.est$Q2
    K.est <- Bfull.est$K

    B.ests[[itr]] <- B.est
    Q2.ests[[itr]] <- Q2.est
    K.ests[[itr]] <- K.est
    ind.inside <- intersect(ind.inside,ind.inside.est)
  }
  list(B.ests=B.ests,Q2.ests=Q2.ests,K.ests=K.ests,
       ind.inside=ind.inside)
}
