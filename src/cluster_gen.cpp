#include <Rcpp.h>
/// Function that generates clusters of states and calculates frequencies for
/// each cluster.
///
/// @param x a vector of coded values
/// @param y a vector of coded values
/// @param prog if TRUE, lag of x (Markov order) is increased by one
/// @param lx Markov order of x
/// @param ly Markov order of y
///
/// @return returns a list with clusters and associated frequencies
/// @keywords internal
/// @export
///
/// @examples
///
// [[Rcpp::export]]
Rcpp::NumericVector cluster_gen(Rcpp::IntegerVector x, int lx = 1,
                       Rcpp::Nullable<Rcpp::IntegerVector> y = R_NilValue,
                       Rcpp::Nullable<int> ly = R_NilValue,
                       bool prog = true) {

  // treat a possible missing y-series
  Rcpp::IntegerVector yy;
  int ly_ = 0;
  bool no_y = y.isNull();

  if (!no_y) {
    yy = Rcpp::as<Rcpp::IntegerVector>(y);
    ly_ = Rcpp::as<int>(ly);
  }

  std::map<std::vector<int>, int> counts;
  int nclust = x.size() - std::max(lx, ly_);

  if (no_y) ly = 0;

  int x_pre;
  if (no_y) {
    x_pre = 0;
  } else {
    x_pre = std::max(ly_ - lx, 0);
  }
  int x_post = lx;
  if (!no_y) x_post += std::max(ly_ - lx, 0);
  if (!prog) x_post -= 1;

  int y_pre = std::max(lx - ly_, 0);
  int y_post = std::max(lx, ly_) - 1;

  // iterate over the clusters and generate the cluster-sequences
  if (no_y) {
    for (int i = 0; i < nclust; ++i) {
      std::vector<int> v;
      v.reserve(x_post - x_pre + 1);

      for (int x_ = x_pre; x_ <= x_post; ++x_) {
        v.push_back(x[i + x_]);
      }
      counts[v]++;
    }
  } else {
    for (int i = 0; i < nclust; ++i) {
      std::vector<int> v;
      v.reserve(x_post - x_pre + y_post - y_pre + 2);

      for (int x_ = x_pre; x_ <= x_post; ++x_) {
        v.push_back(x[i + x_]);
      }
      for (int y_ = y_pre; y_ <= y_post; ++y_) {
        v.push_back(yy[i + y_]);
      }
      counts[v]++;
    }
  }


  // flatten the map to a vector
  // contains the cluster names, i.e., "1 2 2", "2 2 2" etc
  std::vector<std::string> cluster_names;
  cluster_names.reserve(counts.size());

  // contains the cluster counts
  std::vector<int> cluster_counts;
  cluster_counts.reserve(counts.size());

  double total = 0;
  for (auto it = counts.begin(); it != counts.end(); ++it) {
    total += it->second;
    std::string cluster_name;
    for (auto ii = it->first.begin(); ii != it->first.end(); ++ii) {
      cluster_name += std::to_string(*ii) + ' ';
    }
    cluster_names.push_back(cluster_name.substr(0, cluster_name.size() - 1));
    cluster_counts.push_back(it->second);
  }
  Rcpp::NumericVector vec(cluster_names.size());
  for (int i = 0; i < cluster_names.size(); ++i) {
    vec[i] = cluster_counts[i] / total;
  }

  vec.attr("names") = cluster_names;

  return vec;
}
