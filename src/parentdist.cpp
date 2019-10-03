#include <Rcpp.h>
using namespace Rcpp;

// Function to get a per-locus estimate of probability of sampling two different
// alleles (without replacement) from a parental genotype.
// [[Rcpp::export]]
NumericVector HoOneParent(IntegerVector genotypes, IntegerVector alleles2loc,
                          IntegerVector keeploc, double ploidy){
  int nloc = keeploc.size();
  int L;
  IntegerVector thisgen;
  NumericVector out(nloc, 1.0);
  
  for(int i = 0; i < nloc; i++){
    L = keeploc[i];
    thisgen = genotypes[alleles2loc == L];
    for(int a = 0; a < thisgen.size(); a++){
      out[i] -= (thisgen[a] / ploidy) * ((thisgen[a] - 1)/(ploidy - 1));
    }
  }
  
  return out;
}

// Function to get a per-locus estimate of the probability of sampling two
// different alleles if one is selected from each of two parental genotypes.
// [[Rcpp::export]]
NumericVector HoTwoParents(IntegerVector genotypes1, IntegerVector genotypes2,
                           IntegerVector alleles2loc, IntegerVector keeploc,
                           double ploidy){
  int nloc = keeploc.size();
  int L;
  IntegerVector thisgen1;
  IntegerVector thisgen2;
  NumericVector out(nloc, 1.0);
  
  for(int i = 0; i < nloc; i++){
    L = keeploc[i];
    thisgen1 = genotypes1[alleles2loc == L];
    thisgen2 = genotypes2[alleles2loc == L];
    for(int a = 0; a < thisgen1.size(); a++){
      out[i] -= (thisgen1[a] / ploidy) * (thisgen2[a]/ploidy);
    }
  }
  
  return out;
}
