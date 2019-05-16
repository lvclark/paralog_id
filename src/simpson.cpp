#include <Rcpp.h>
using namespace Rcpp;

// Function to get the Simpson index, for a group of sequencing reads or
// anything else.
// [[Rcpp::export]]
double Simpson(IntegerVector counts) {
  double N = sum(counts);
  double out = 1.0;
  
  for(int i = 0; i < counts.size(); i++){
    out -= pow(counts[i] / N, 2);
  }
  
  return out;
}

// Function to get the Simpson index on allelic read depth as encoded in polyRAD.
// Returns a matrix of taxa x loci, as well as a vector for the total across taxa.
// [[Rcpp::export]]
List HReadDepth(IntegerMatrix alleleDepth, IntegerVector alleles2loc, int nLoci){
  int nTaxa = alleleDepth.nrow();
  IntegerVector alleles = seq(0, alleles2loc.size() - 1);
  IntegerVector thesecol;
  IntegerVector thesecounts;
  IntegerVector thesetotcounts;
  int thesenal;
  NumericMatrix Hind(nTaxa, nLoci);
  NumericVector He(nLoci);
  
  for(int L = 0; L < nLoci; L++){
    thesecol = alleles[alleles2loc == L + 1];
    thesenal = thesecol.size();
    thesecounts = IntegerVector(thesenal);
    thesetotcounts = IntegerVector(thesenal);
    
    for(int t = 0; t < nTaxa; t++){
      // For this taxon and locus, copy to temporary vector
      // and add to vector for total across taxa.
      for(int a = 0; a < thesenal; a++){
        thesecounts[a] = alleleDepth(t, thesecol[a]);
        thesetotcounts[a] += thesecounts[a];
      }
      Hind(t, L) = Simpson(thesecounts); // H for t x L
    }
    He[L] = Simpson(thesetotcounts); // H for locus
  }
  
  List outlist = List::create(Hind, He);
  return outlist;
}
