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
