// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat c_ecdf(arma::vec dat,
                 arma::vec pts) {
  int n_pts = pts.size();
  arma::mat results(n_pts, 2);
  arma::uvec t_vec(dat.size());

  
  for(int d = 0; d < n_pts; d++){
    // for(int t = 0; t < t_vec.size(); t++){
    //   t_vec[t] = dat[t] < pts[d];
    // }

    t_vec = (dat <= pts[d]);
    results(d,0) = pts[d];
    results(d,1) = std::accumulate(t_vec.begin(), t_vec.end(), 0.0) / dat.size();
  }
  return results;
}

// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat c_ecdf2(arma::vec dat,
                  arma::vec pts) {
  int n_pts = pts.size();
  arma::mat results(n_pts, 2);
  arma::uvec t_vec(dat.size());
  int d = 0;
  for(arma::vec::iterator it = pts.begin(); it != pts.end(); ++it){
    t_vec = (dat <= *it);
    results(d,0) = *it;
    results(d,1) = std::accumulate(t_vec.begin(), t_vec.end(), 0.0) / dat.size();
    d++;
  }
  return results;
}



// // [[Rcpp::depends(RcppArmadillo)]]
// # include <RcppArmadillo.h>
// 
// // [[Rcpp::export]]
// arma::mat c_wcdf(arma::vec dat,
//                  arma::vec ps,
//                  arma::vec ta,
//                  arma::vec pts) {
//   int n_pts = pts.size();
//   arma::mat results(n_pts, 2);
//   arma::uvec t_vec(dat.size());
//   
//   for(int d = 0; d < n_pts; d++){
//     t_vec = (dat <= pts[d]);
//     results(d,0) = pts[d];
//     results(d,1) = std::accumulate(t_vec.begin(), t_vec.end(), 0.0) / dat.size();
//   }
//   return results;
// }
