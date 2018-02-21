#include <Rcpp.h>

/// Creates the text-clusters using x, lx, y, ly, and prog,
/// i.e., "1 2 2", "2 2 2", etc
std::vector<std::string> generate_clusters(Rcpp::IntegerVector x,
                                           int lx,
                                           Rcpp::IntegerVector y,
                                           int ly,
                                           bool prog = true,
                                           bool no_y = false) {

  int nclust = x.size() - std::max(lx, ly);
  std::vector<std::string> vec(nclust);
  int lxi, lyi;

  if (no_y) ly = 0;

  int x_pre;
  if (no_y) {
    x_pre = 0;
  } else {
    x_pre = std::max(ly - lx, 0);
  }
  int x_post = lx;
  if (!no_y) x_post += std::max(ly - lx, 0);
  if (!prog) x_post -= 1;

  int y_pre = std::max(lx - ly, 0);
  int y_post = std::max(lx, ly) - 1;

  // iterate over the clusters and generate the cluster-sequences
  if (no_y) {
    for (int i = 0; i < nclust; ++i) {
      std::string str = "";

      for (int x_ = x_pre; x_ <= x_post; ++x_) {
        str += std::to_string(x[i + x_]);
        str += " ";
      }

      vec[i] = str.substr(0, str.size() - 1);
    }
  } else {
    for (int i = 0; i < nclust; ++i) {
      std::string str = "";

      for (int x_ = x_pre; x_ <= x_post; ++x_) {
        str += std::to_string(x[i + x_]);
        str += " ";
      }
      for (int y_ = y_pre; y_ <= y_post; ++y_) {
        str += std::to_string(y[i + y_]);
        str += " ";
      }

      vec[i] = str.substr(0, str.size() - 1);
    }
  }

  return vec;
}

// Calculates the frequency for a given input of strings
Rcpp::NumericVector freq_table(std::vector<std::string> x) {
  std::map<std::string, int> counts;

  // load the values into the count-map
  for (int i = 0; i < x.size(); ++i) {
    ++counts[x[i]];
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
Rcpp::List cluster_gen(Rcpp::IntegerVector x, int lx = 1,
                       Rcpp::Nullable<Rcpp::IntegerVector> y = R_NilValue,
                       Rcpp::Nullable<int> ly = R_NilValue,
                       bool prog = true) {

  // treat a possible missing y-series
  Rcpp::IntegerVector y_;
  int ly_ = 0;
  bool no_y = y.isNull();

  if (!no_y) {
    y_ = Rcpp::as<Rcpp::IntegerVector>(y);
    ly_ = Rcpp::as<int>(ly);
  }

  std::vector<std::string> clusters;
  clusters = generate_clusters(x, lx, y_, ly_, prog, no_y);

  Rcpp::NumericVector freq;
  freq = freq_table(clusters);

  return Rcpp::List::create(Rcpp::Named("cluster") = clusters,
                            Rcpp::Named("frequency") = freq);
}
