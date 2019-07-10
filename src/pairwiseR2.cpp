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
  
  mat Y = square(X);
  rowvec Ysum = sum(Y);
  mat Z = 2 * X.t() * X;
  
  mat r2(genes_count, genes_count, fill::zeros);
  vec tsss(genes_count);
  
  for (int i = 0; i < genes_count; i++) {
    double mean_x_intercept = mean(X.col(i));
    vec x_int = X.col(i) - mean_x_intercept;
    double tss = dot(x_int, x_int);
    tsss[i] = tss;
  }
  
  for (int i = 0; i < genes_count; i++) {
    vec rss = Ysum[i] + (Ysum - Z.row(i)).t();
    r2.col(i) += (ones(genes_count) - rss / tsss) / 2.0;
    r2.row(i) += ((ones(genes_count) - rss / tsss) / 2.0).t();
  }
  
  return r2;
}

