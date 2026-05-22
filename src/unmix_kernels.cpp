#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

arma::mat safe_inverse_cpp(const arma::mat& mat, const double tol = 1e-10) {
    arma::mat inv_mat;
    bool ok = arma::inv(inv_mat, mat);
    if (ok) {
        return inv_mat;
    }
    return arma::pinv(mat, tol);
}

arma::vec nnls_lawson_hanson_cpp(const arma::mat& A,
                                 const arma::vec& b,
                                 const double tol = 1e-10,
                                 const int max_outer = 500,
                                 const int max_inner = 500) {
    const arma::uword p = A.n_cols;
    arma::vec x(p, arma::fill::zeros);
    arma::uvec passive(p, arma::fill::zeros);

    arma::vec w = A.t() * (b - A * x);
    int outer_iter = 0;

    while (outer_iter < max_outer) {
        arma::uvec active_idx = arma::find((passive == 0) && (w > tol));
        if (active_idx.n_elem == 0) {
            break;
        }

        arma::uword t = active_idx(arma::index_max(w.elem(active_idx)));
        passive(t) = 1;

        int inner_iter = 0;
        while (inner_iter < max_inner) {
            arma::uvec passive_idx = arma::find(passive == 1);
            arma::vec s(p, arma::fill::zeros);

            if (passive_idx.n_elem > 0) {
                arma::mat AP = A.cols(passive_idx);
                arma::mat normal_mat = AP.t() * AP;
                arma::vec rhs = AP.t() * b;
                arma::vec z;
                bool solved = arma::solve(z, normal_mat, rhs);
                if (!solved) {
                    z = arma::pinv(normal_mat, tol) * rhs;
                }
                s.elem(passive_idx) = z;
            }

            arma::uvec nonpos_idx = arma::find((passive == 1) && (s <= tol));
            if (nonpos_idx.n_elem == 0) {
                x = s;
                break;
            }

            double alpha = 1.0;
            for (arma::uword k = 0; k < nonpos_idx.n_elem; ++k) {
                arma::uword idx = nonpos_idx(k);
                double denom = x(idx) - s(idx);
                if (denom > 0) {
                    alpha = std::min(alpha, x(idx) / denom);
                }
            }

            x = x + alpha * (s - x);

            arma::uvec reset_idx = arma::find((passive == 1) && (x <= tol));
            if (reset_idx.n_elem > 0) {
                x.elem(reset_idx).zeros();
                passive.elem(reset_idx).zeros();
            }

            ++inner_iter;
        }

        w = A.t() * (b - A * x);
        ++outer_iter;
    }

    x.elem(arma::find(x < 0)).zeros();
    return x;
}

} // namespace

// [[Rcpp::export]]
arma::mat spectreasy_wls_unmix_cpp(const arma::mat& Y,
                                   const arma::mat& M,
                                   const double background_noise = 25,
                                   const double tol = 1e-10) {
    const arma::uword n_cells = Y.n_rows;
    const arma::uword n_markers = M.n_rows;
    arma::mat A(n_cells, n_markers, arma::fill::zeros);

    arma::mat Mt = M.t();
    arma::mat MMt = M * Mt;
    arma::mat MMt_inv = safe_inverse_cpp(MMt, tol);
    arma::mat fallback_linear = MMt_inv * M;

    for (arma::uword i = 0; i < n_cells; ++i) {
        arma::rowvec y = Y.row(i);
        arma::rowvec weights = 1.0 / (arma::clamp(y, 0.0, arma::datum::inf) + background_noise);
        arma::mat M_weighted = M.each_row() % weights;
        arma::mat MWMt = M_weighted * Mt;
        arma::vec rhs = M * (weights.t() % y.t());

        arma::vec coeffs;
        bool solved = arma::solve(coeffs, MWMt, rhs);
        if (!solved || !coeffs.is_finite()) {
            coeffs = fallback_linear * y.t();
        }
        A.row(i) = coeffs.t();
    }

    return A;
}

// [[Rcpp::export]]
arma::mat spectreasy_nnls_unmix_cpp(const arma::mat& Y,
                                    const arma::mat& M,
                                    const double tol = 1e-10,
                                    const int max_outer = 500,
                                    const int max_inner = 500) {
    const arma::uword n_cells = Y.n_rows;
    const arma::uword n_markers = M.n_rows;
    arma::mat A_out(n_cells, n_markers, arma::fill::zeros);
    arma::mat A = M.t();

    for (arma::uword i = 0; i < n_cells; ++i) {
        arma::vec b = Y.row(i).t();
        arma::vec x = nnls_lawson_hanson_cpp(A, b, tol, max_outer, max_inner);
        A_out.row(i) = x.t();
    }

    return A_out;
}

// [[Rcpp::export]]
arma::mat spectreasy_unmix_best_af_cpp(const arma::mat& Y,
                                       const arma::mat& M,
                                       const arma::uvec& fluor_idx,
                                       const arma::uvec& af_idx,
                                       const std::string& method,
                                       const double background_noise = 25,
                                       const double tol = 1e-10,
                                       const int max_outer = 500,
                                       const int max_inner = 500) {
    const arma::uword n_cells = Y.n_rows;
    const arma::uword n_detectors = Y.n_cols;
    const arma::uword n_markers = M.n_rows;
    const arma::uword n_fluors = fluor_idx.n_elem;
    const arma::uword n_af = af_idx.n_elem;

    arma::mat A_out(n_cells, n_markers, arma::fill::zeros);

    // Precompute projection matrices P_k for each candidate AF band k.
    std::vector<arma::mat> P_list(n_af);
    for (arma::uword k = 0; k < n_af; ++k) {
        arma::mat X_k(n_fluors + 1, n_detectors);
        for (arma::uword f = 0; f < n_fluors; ++f) {
            X_k.row(f) = M.row(fluor_idx(f));
        }
        X_k.row(n_fluors) = M.row(af_idx(k));

        arma::mat A_k = X_k.t();
        arma::mat AtA = A_k.t() * A_k;
        arma::mat AtA_inv = safe_inverse_cpp(AtA, tol);
        P_list[k] = A_k * AtA_inv * A_k.t();
    }

    for (arma::uword i = 0; i < n_cells; ++i) {
        arma::vec y = Y.row(i).t();

        // 1. Find the best AF band index k* that minimizes OLS RSS.
        double min_rss = std::numeric_limits<double>::max();
        arma::uword best_k = 0;
        for (arma::uword k = 0; k < n_af; ++k) {
            arma::vec resid = y - P_list[k] * y;
            double rss = arma::dot(resid, resid);
            if (rss < min_rss) {
                min_rss = rss;
                best_k = k;
            }
        }

        // 2. Build the model matrix X for the selected best_k.
        arma::mat X_best(n_fluors + 1, n_detectors);
        for (arma::uword f = 0; f < n_fluors; ++f) {
            X_best.row(f) = M.row(fluor_idx(f));
        }
        X_best.row(n_fluors) = M.row(af_idx(best_k));
        arma::mat A_best = X_best.t();

        // 3. Unmix the cell using the selected model.
        arma::vec coeffs(n_fluors + 1, arma::fill::zeros);
        if (method == "OLS") {
            arma::mat AtA = A_best.t() * A_best;
            coeffs = safe_inverse_cpp(AtA, tol) * A_best.t() * y;
        } else if (method == "WLS") {
            arma::vec weights = 1.0 / (arma::clamp(y, 0.0, arma::datum::inf) + background_noise);
            arma::mat M_weighted = X_best.each_row() % weights.t();
            arma::mat MWMt = M_weighted * A_best;
            arma::vec rhs = X_best * (weights % y);
            bool solved = arma::solve(coeffs, MWMt, rhs);
            if (!solved || !coeffs.is_finite()) {
                arma::mat AtA = A_best.t() * A_best;
                coeffs = safe_inverse_cpp(AtA, tol) * A_best.t() * y;
            }
        } else if (method == "NNLS") {
            coeffs = nnls_lawson_hanson_cpp(A_best, y, tol, max_outer, max_inner);
        }

        // 4. Map the resulting coefficients back to the full reference matrix layout.
        for (arma::uword f = 0; f < n_fluors; ++f) {
            A_out(i, fluor_idx(f)) = coeffs(f);
        }
        A_out(i, af_idx(best_k)) = coeffs(n_fluors);
    }

    return A_out;
}
