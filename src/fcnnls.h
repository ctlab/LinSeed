#include <RcppArmadillo.h>
#include <random>
#include <math.h>
#include <algorithm>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/fisher_f.hpp>
using namespace boost::math;
using namespace Rcpp;
using namespace arma;

arma::mat cssls(const arma::mat& CtC, const arma::mat& CtA, bool pseudo);
arma::mat cssls(const arma::mat& CtC, const arma::mat& CtA, const arma::umat& Pset, bool pseudo);
arma::mat fcnnls_c(const arma::mat& C, const arma::mat& A);
arma::mat fcnnls_sum_to_one(const arma::mat& C, const arma::mat& A, double delta);
