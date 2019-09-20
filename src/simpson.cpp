#include <Rcpp.h>
using namespace Rcpp;

// Function to get the Simpson index, for a group of sequencing reads or
// anything else.
// [[Rcpp::export]]
double GiniSimpson(NumericVector counts) {
  double N = sum(counts);
  double out = 1.0;
  
  for(int i = 0; i < counts.size(); i++){
    out -= pow(counts[i] / N, 2);
  }
  
  return out;
}

// Function to get Hind/He on allelic read depth as encoded in polyRAD.
// Returns a vector of values, individual x locus
// [[Rcpp::export]]
NumericMatrix HindHeMat(IntegerMatrix alleleDepth, NumericMatrix depthRatio,
                        IntegerVector alleles2loc, int nLoci){
  int nTaxa = alleleDepth.nrow();
  IntegerVector alleles = seq(0, alleles2loc.size() - 1);
  IntegerVector thesecol;
  NumericVector thesecounts;
  int thesenal;
  int thisdepth;
  NumericVector thesefreq;
  double thisHind;
  double thisHe;
  NumericMatrix HindHe(nTaxa, nLoci);
  
  for(int L = 0; L < nLoci; L++){
    thesecol = alleles[alleles2loc == L + 1];
    thesenal = thesecol.size();
    thesecounts = NumericVector(thesenal);
    thesefreq = NumericVector(thesenal);
    // fill in allele frequencies from depth ratios
    for(int a = 0; a < thesenal; a++){
      thesefreq[a] = mean(na_omit(depthRatio( _ , thesecol[a])));
    }
    thisHe = GiniSimpson(thesefreq);
    
    for(int t = 0; t < nTaxa; t++){
      // For this taxon and locus, copy to temporary vector
      // and add to vector for total across taxa.
      for(int a = 0; a < thesenal; a++){
        thesecounts[a] = alleleDepth(t, thesecol[a]);
      }
      thisdepth = sum(thesecounts);
      if(thisdepth < 2){ // can't calculate for depth below 2
        HindHe(t, L) = R_NaN;
      } else {
        thisHind = GiniSimpson(thesecounts); // H for t x L
        HindHe(t, L) = thisHind / thisHe * thisdepth / (thisdepth - 1);
      }
    }
  }
  
  return HindHe;
}
