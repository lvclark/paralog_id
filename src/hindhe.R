# wrapper function to estimate Hind/He for each locus in a RADdata object.
Rcpp::sourceCpp("simpson.cpp")

HindHe <- function(object, ...){
  UseMethod("HindHe", object)
}
HindHe.RADdata <- function(object, omitTaxa = GetBlankTaxa(object), ...){
  taxa <- GetTaxa(object)[!GetTaxa(object) %in% omitTaxa]
  
  hindhe <- HindHeByLoc(object$alleleDepth[taxa,, drop = FALSE],
                        object$depthRatio[taxa,, drop = FALSE],
                        object$alleles2loc, nLoci(object))
  names(hindhe) <- GetLoci(object)
  return(hindhe)
}
