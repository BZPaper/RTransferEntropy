#include <Rcpp.h>
#include <sstream>
using namespace Rcpp;

// Wrapper function that calculates the transition probabilities for a given
// vector
// [[Rcpp::export]]
List calculate_transition_probabilities(IntegerVector x, int lx = 1) {
  if (lx > 20) {
    Rprintf("in calculate_transition_probabilities lx > 20 - Aborting calculation\n");
    List rr;
    return rr;
  }

  int cluster_id; // i.e., 1, 2, 3 (the first letter of the cluster_val)

  std::map<int, std::map<std::vector<int>, int>> counts;

  // Fill the count map map
  for (int i = 0; i < x.size() - lx; ++i) {
    cluster_id = x[i];
    // construct the clusterValue
    std::vector<int> cluster_val; // i.e., 21, 12, 11, 22, 231 (lx = 2) etc.
    cluster_val.reserve(lx + 1);
    for (int l = 0; l < lx + 1; ++l) {
      cluster_val.push_back(x[i + l]);
    }
    counts[cluster_id][cluster_val] += 1;
  }

  // check that all value in x are covered; no values between N - lx and N are missed
  std::vector<int> unique_vals;
  for (auto it = counts.begin(); it != counts.end(); ++it)
    unique_vals.push_back(it->first);
  std::vector<int> miss;

  // check the remaining values so that we dont get an out of bounds error later
  for (int i = x.size() - lx; i < x.size(); ++i) {
    cluster_id = x[i];
    // if the cluster_id is not found in the transition probs; create an empty one
    if (counts.find(cluster_id) == counts.end()) {
      unique_vals.push_back(cluster_id);
      miss.push_back(cluster_id);
    }
  }

  // add values that are between N - lx and N to the counts
  for (int v : miss) {
    counts[v][unique_vals] = 0;
  }

  // add missing elements, i.e., if there is no 2 in the series,
  // replace the element with an empty map
  for (int i = counts.begin()->first; i <= counts.rbegin()->first; ++i) {
    auto it = counts.find(i);
    if (it == counts.end()) {
      counts[i] = std::map<std::vector<int>, int>();
    }
  }



  // flatten map map to a vector of vectors
  std::vector<int> names; // contains the names of the clusters
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

      if (n == 0) {
        tmp_vec.push_back(0.0);
      } else {
        tmp_vec.push_back((double) val_ptr->second / n);
      }

      std::stringstream name_val;
      for (auto it = val_ptr->first.begin();
           it != val_ptr->first.end();
           ++it) {
        name_val << *it << ' ';
      }
      name_vec.push_back(name_val.str().substr(0, name_val.str().size() - 1));
    }

    names.push_back(id_ptr->first);
    tmp_vec.attr("names") = name_vec;
    res.push_back(tmp_vec);
  }
  List ret_list = wrap(res);
  ret_list.names() = names;

  return ret_list;
}
