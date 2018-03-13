// You should have received a copy of the GNU General Public License
// along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.
// [[Rcpp::depends(BH)]]

#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat pairwiseR2(const arma::mat& X) {
  int genes_count = X.n_cols;
  int sample_count = X.n_rows;
  
  mat r2(genes_count, genes_count, fill::zeros);
  vec tsss(genes_count);
  
  vec sxx(genes_count);
  
  for (int i = 0; i < genes_count; i++) {
    sxx[i] = dot(X.col(i), X.col(i));
    
    vec x_int(sample_count);
    double mean_x_intercept = mean(X.col(i));
    for (int j = 0; j < sample_count; j++) {
      x_int[j] = X.at(j, i) - mean_x_intercept;
    }
    double tss = dot(x_int, x_int);
    tsss[i] = tss;
  }
  
  for (int i = 0; i < genes_count; i++) {
    rowvec sxy = X.col(i).t() * X;
    
    vec rss(genes_count);
    for (int j = 0; j < genes_count; j++) {
      vec z = X.col(j) - X.col(i);
      rss[j] = dot(z, z);
    }
    r2.col(i) += (ones(genes_count) - rss / tsss) / 2.0;
    r2.row(i) += ((ones(genes_count) - rss / tsss) / 2.0).t();
  }
  
  return r2;
}

