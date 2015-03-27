# Fisher's method to combine p-values:

fishers <- function(p1,p2){
    chi.sq <- -2*(log(p1)+log(p2))
    p.merge <- sapply(chi.sq,pchisq,df=4,lower=FALSE)
    return(p.merge)
  }
