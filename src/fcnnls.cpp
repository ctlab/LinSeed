#include "fcnnls.h"

arma::mat cssls(const arma::mat& CtC,
                const arma::mat& CtA,
                bool pseudo) {
    mat K(CtC.n_rows, CtC.n_cols);
    K.fill(0);
    if (pseudo) {
        K = pinv(CtC) * CtA;
    } else {
        K = solve(CtC, diagmat(ones(CtC.n_rows))) * CtA;
    }
    return K;
}

arma::mat cssls(const arma::mat& CtC,
                const arma::mat& CtA,
                const arma::umat& Pset,
                bool pseudo) {
    mat K(CtA.n_rows, CtA.n_cols);
    K.fill(0);

    int l = Pset.n_rows;
    int RHSc = Pset.n_cols;
    rowvec lvec(l);
    for (int i = 1; i <= l; i++) {
        lvec[i - 1] = 1 << (l - i);
    }

    mat codedPset = lvec * Pset;
    vec sortedPset = arma::sort(codedPset.row(0).t());
    uvec sortedEset = sort_index(codedPset.row(0).t());
    vec breaks = diff(sortedPset);
    uvec bbreaks = find(breaks) + 1;
    uvec breakIdx(bbreaks.n_elem + 2);
    breakIdx[0] = 0;
    if (bbreaks.n_elem > 0) {
        breakIdx.subvec(1, bbreaks.n_elem) = bbreaks;
    }

    breakIdx[bbreaks.n_elem + 1] = RHSc;

    for (int k = 0; k < breakIdx.n_elem - 1; k++) {
        uvec cols2solve = sortedEset.subvec(breakIdx[k], breakIdx[k + 1] - 1);
        uvec vars = find(Pset.col(sortedEset[breakIdx[k]]));
        if (vars.n_elem == 0) break;
        //K.submat(vars, cols2solve) = solve(CtC.submat(vars, vars), diagmat(ones(vars.n_elem))) * CtA.submat(vars, cols2solve);
        K.submat(vars, cols2solve) = pinv(CtC.submat(vars, vars)) * CtA.submat(vars, cols2solve);
    }
    return K;
}
// [[Rcpp::export]]
arma::mat fcnnls_c(const arma::mat& C,
                   const arma::mat& A) {
    int n_obs = C.n_rows;
    int l_var = C.n_cols;
    if (C.n_rows != A.n_rows) {

        exit(1);
    }
    int pRHS = A.n_cols;
    mat W = zeros(l_var, pRHS);

    int iter = 0;
    int max_iter = 3 * l_var;

    mat CtC = C.t() * C;
    mat CtA = C.t() * A;
    mat K = cssls(CtC, CtA, false);
    umat Pset = K > zeros(K.n_rows, K.n_cols);
    K.elem(find(K <= zeros(K.n_rows, K.n_cols))).zeros();
    mat D = K;
    uvec Fset = find(any(Pset < ones(Pset.n_rows, Pset.n_cols)));
    while (Fset.n_elem > 0) {

        K.cols(Fset) = cssls(CtC, CtA.cols(Fset), Pset.cols(Fset), false);
        uvec Hset = Fset.elem( find(any(K.cols(Fset) < zeros(K.n_rows, Fset.n_elem))) );

        if (Hset.n_elem > 0) {
            int nHset = Hset.n_elem;
            mat alpha(l_var, nHset, fill::zeros);
            while (Hset.n_elem > 0 && iter < max_iter) {
                iter++;
                alpha.cols(0, nHset - 1).fill(datum::inf);
                mat twos(Pset.n_rows, Hset.n_elem);
                twos.fill(2);
                uvec ij = find(Pset.cols(Hset) + (K.cols(Hset) < 0) == twos);
                uvec is(ij.n_elem), js(ij.n_elem);
                for (int i = 0; i < ij.n_elem; i++) {
                    is[i] = ij[i] % Pset.n_rows;
                    js[i] = ij[i] / Pset.n_rows;
                }

                uvec hIdx = ij;
                uvec negIdx = Hset.elem(js) * l_var + is;

                alpha.elem(hIdx) = D.elem(negIdx) / (D.elem(negIdx) - K.elem(negIdx));

                urowvec minIdx = index_min(alpha.cols(0, nHset - 1));
                rowvec alphaMin = arma::min(alpha.cols(0, nHset - 1));
                for (int i = 0; i < l_var; i ++) {
                    alpha.cols(0, nHset-1).row(i) = alphaMin;
                }


                D.cols(Hset) = D.cols(Hset) - alpha.cols(0, nHset - 1) % (D.cols(Hset) - K.cols(Hset));
                uvec idx2zero = Hset * l_var + minIdx.t();
                D.elem(idx2zero).zeros();
                Pset.elem(idx2zero).zeros();
                K.cols(Hset) = cssls(CtC, CtA.cols(Hset), Pset.cols(Hset), false);
                Hset = find(any(K < zeros(K.n_rows, K.n_cols)));
                nHset = Hset.n_elem;
            }
        }

        // if (iter == max_iter) {

        //     return K;
        // }
        W.cols(Fset) = CtA.cols(Fset) - CtC * K.cols(Fset);
        mat tmp = (ones(l_var, Fset.n_elem) - Pset.cols(Fset)) % W.cols(Fset);
        uvec Jset = find(all(tmp <= zeros(l_var, Fset.n_elem)));
        uvec FsetDiff(Fset.n_elem - Jset.n_elem);
        uvec toDiff = Fset.elem(Jset);
        std::set_difference(Fset.begin(), Fset.end(),
                            toDiff.begin(), toDiff.end(),
                            FsetDiff.begin());

        Fset = FsetDiff;
        if (Fset.n_elem > 0) {
            tmp = (ones(l_var, Fset.n_elem) - Pset.cols(Fset)) % W.cols(Fset);
            urowvec mxidx = index_max(tmp);
            Pset.elem(Fset * l_var + mxidx.t()).ones();
            // Pset = K >= zeros(K.n_rows, K.n_cols);
            D.cols(Fset) = K.cols(Fset);
        }

    }
    return K;

}

// [[Rcpp::export]]
arma::mat fcnnls_sum_to_one(const arma::mat& C,
                            const arma::mat& A,
                            double delta = 1) {
    mat C_copy(C.n_rows + 1, C.n_cols);
    mat A_copy(A.n_rows + 1, A.n_cols);
    C_copy.rows(0, C.n_rows - 1) = C * delta;
    A_copy.rows(0, A.n_rows - 1) = A * delta;
    C_copy.row(C.n_rows).fill(1);
    A_copy.row(A.n_rows).fill(1);
    return(fcnnls_c(C_copy, A_copy));
}
