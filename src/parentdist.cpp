#include <Rcpp.h>
using namespace Rcpp;

// Function to take vectors of read counts from two different individuals and
// return the probability that a read sampled from one would be different from
// a read sampled from the other.  Used for examining the parents of a mapping
// population.
// If there is backcrossing, counts1 is assumed to come from the recurrent
// parent.
// [[Rcpp::export]]
double ParentDist(IntegerVector counts1, IntegerVector counts2,
                  int gen_backcrossing = 0, int gen_selfing = 0){
  double tot1 = sum(counts1);
  double tot2 = sum(counts2);
  double out = 1;
  double freq1;
  double freq2;
  
  if(tot1 == 0 || tot2 == 0){
    return R_NaN;
  }
  
  for(int i = 0; i < counts1.size(); i++){
    freq1 = counts1[i]/tot1;
    freq2 = counts2[i]/tot2;
    for(int g = 0; g < gen_backcrossing; g++){
      freq2 = (freq1 + freq2)/2;
    }
    out -= freq1 * freq2;
  }
  for(int g = 0; g < gen_selfing; g++){
    out /= 2;
  }
  
  return out;
}

// Function to determine prob of sampling two alleles from recurrent parent,
// two alleles from donor parent, or two alleles from different parents, with
// replacment, in a progeny.
// [[Rcpp::export]]
NumericVector ProgAlProbs(int ploidy = 2, int gen_backcrossing = 0, int gen_selfing = 0){
  NumericVector probs(3);
  int hap = ploidy / 2;
  double shiftval;
  int i;
  
  // Probabilities for F1
  probs[0] = 0.5 * (hap - 1) / (ploidy - 1); // both drawn from recurrent parent
  probs[1] = probs[0];                       // both drawn from donor parent
  probs[2] = 1 - probs[0] - probs[1];        // drawn from different parents
  
  // Backcross
  for(i = 0; i < gen_backcrossing; i++){
    shiftval = probs[2] / 2;
    probs[2] = shiftval;
    probs[0] += shiftval;
  }
  
  return probs;
}