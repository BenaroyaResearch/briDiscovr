#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

//' Compute jaccard coefficient between nearest-neighbor sets
//'
//' Weights of both i->j and j->i are recorded if they have intersection. In this case
//' w(i->j) should be equal to w(j->i). In some case i->j has weights while j<-i has no
//' intersections, only w(i->j) is recorded. This is determinded in code `if(u>0)`.
//' In this way, the undirected graph is symmetrized by halfing the weight
//' in code `weights(r, 2) = u/(2.0*ncol - u)/2`.
//'
//' Rows of the returned matrix beyond the last recorded link are left as zeros.
//'
//' Author: Chen Hao, Date: 25/09/2015
//' Copied from https://github.com/JinmiaoChenLab/Rphenograph/blob/master/src/jaccard_coeff.cpp
//' for use in briDiscovr on 200611 under the Artistic-2.0 license.
//' Rewritten in briDiscovr to count intersections with a marker array instead of
//' Rcpp::intersect(); output (values and row order) is unchanged.
//' @param idx NumericMatrix of 1-based neighbor indices (one row per event, self excluded)
//' @return NumericMatrix
//' @export
// [[Rcpp::export]]
NumericMatrix jaccard_coeff(NumericMatrix idx) {
  const int nrow = idx.nrow(), ncol = idx.ncol();
  NumericMatrix weights(nrow * ncol, 3);

  // copy neighbor indices to a row-major, 0-based integer array for cache-friendly access
  std::vector<int> nn((size_t) nrow * ncol);
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      const int v = (int) idx(i, j) - 1;
      if (v < 0 || v >= nrow) stop("neighbor indices must be between 1 and nrow(idx)");
      nn[(size_t) i * ncol + j] = v;
    }
  }

  // mark[x] == i  <=>  x is a neighbor of event i
  std::vector<int> mark(nrow, -1);
  size_t r = 0;
  for (int i = 0; i < nrow; i++) {
    const int *nodei = &nn[(size_t) i * ncol];
    for (int j = 0; j < ncol; j++) mark[nodei[j]] = i;
    for (int j = 0; j < ncol; j++) {
      const int k = nodei[j];
      const int *nodek = &nn[(size_t) k * ncol];
      int u = 0;  // count intersection number
      for (int l = 0; l < ncol; l++) if (mark[nodek[l]] == i) u++;
      if (u > 0) {
        weights(r, 0) = i + 1;
        weights(r, 1) = k + 1;
        weights(r, 2) = u / (2.0 * ncol - u) / 2;  // symmetrize the graph
        r++;
      }
    }
  }

  return weights;
}
