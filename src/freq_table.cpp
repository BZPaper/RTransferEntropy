#include <Rcpp.h>

// Calculates the frequency for a given input of strings
// [[Rcpp::export]]
Rcpp::NumericVector freq_table(std::vector<int> x) {
  std::map<int, int> counts;

  // load the values into the count-map
  for (int i = 0; i < x.size(); ++i) {
    ++counts[x[i]];
  }

  // flatten the map to a vector
  // contains the cluster names, i.e., "1 2 2", "2 2 2" etc
  std::vector<int> cluster_names;
  cluster_names.reserve(counts.size());

  // contains the cluster counts
  std::vector<int> cluster_counts;
  cluster_counts.reserve(counts.size());

  double total = 0;
  for (auto it = counts.begin(); it != counts.end(); ++it) {
    total += it->second;
    cluster_names.push_back(it->first);
    cluster_counts.push_back(it->second);
  }
  Rcpp::NumericVector vec(cluster_names.size());
  for (int i = 0; i < cluster_names.size(); ++i) {
    vec[i] = cluster_counts[i] / total;
  }

  vec.attr("names") = cluster_names;
  return vec;
}
