#include <Rcpp.h>
using namespace Rcpp;

// Wrapper function that calculates the transition probabilities for a given
// vector
// [[Rcpp::export]]
List calculate_transition_probabilities(CharacterVector x, int lx = 1) {

  std::string cluster_id; // i.e., 1, 2, 3 (the first letter of the cluster_val)
  std::string cluster_val; // i.e., 21, 12, 11, 22, 231 (lx = 2) etc.
  std::map<std::string, std::map<std::string, int>> counts;

  // Fill the count map map
  for (int i = 0; i < x.size() - lx; ++i) {
    cluster_id = x[i];
    // construct the clusterValue
    cluster_val = "";
    for (int l = 0; l < lx + 1; ++l) {
      cluster_val += x[i + l];
    }
    counts[cluster_id][cluster_val] += 1;
  }

  // flatten map map to a vector of vectors
  std::vector<std::string> names; // contains the names of the clusters
  names.reserve(counts.size());
  std::vector< NumericVector > res; // contains the count of clusters
  res.reserve(counts.size()) ;

  // loop through the outer map (containing the names and the inner maps)
  for (auto id_ptr = counts.begin(); id_ptr != counts.end(); ++id_ptr) {
    // fill the transition probabilities from the respective ID
    NumericVector tmp_vec;
    CharacterVector name_vec;

    // calculate the total number of counts for this ID (current x-value)
    int n = 0;
    for (auto val_ptr = id_ptr->second.begin();
         val_ptr != id_ptr->second.end();
         ++ val_ptr) {
      n += val_ptr->second;
    }

    // calculate the transition probabilities
    for (auto val_ptr = id_ptr->second.begin();
         val_ptr != id_ptr->second.end();
         ++val_ptr) {
      tmp_vec.push_back((double) val_ptr->second / n);
      name_vec.push_back(val_ptr->first);
    }

    names.push_back(id_ptr->first);
    tmp_vec.attr("names") = name_vec;
    res.push_back(tmp_vec);
  }
  List ret_list = wrap(res);
  ret_list.names() = names;

  return ret_list;
}
