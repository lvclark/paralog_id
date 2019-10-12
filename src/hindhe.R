# wrapper function to estimate Hind/He for each locus in a RADdata object.
Rcpp::sourceCpp("simpson.cpp")

HindHe <- function(object, ...){
  UseMethod("HindHe", object)
}
HindHe.RADdata <- function(object, omitTaxa = GetBlankTaxa(object), ...){
  taxa <- GetTaxa(object)[!GetTaxa(object) %in% omitTaxa]
  
  hindhe <- HindHeMat(object$alleleDepth[taxa,, drop = FALSE],
                      object$depthRatio[taxa,, drop = FALSE],
                      object$alleles2loc, nLoci(object), numeric(0))
  colnames(hindhe) <- GetLoci(object)
  rownames(hindhe) <- taxa
  return(hindhe)
}
