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

bool normal_system_is_usable_cpp(const arma::mat& mat, const double tol = 1e-10) {
    if (mat.n_rows == 0 || mat.n_cols == 0 || mat.n_rows != mat.n_cols || !mat.is_finite()) {
        return false;
    }

    double rc = arma::rcond(mat);
    return std::isfinite(rc) && rc >= tol;
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

arma::vec weighted_lsq_coeffs_cpp(const arma::mat& X,
                                  const arma::vec& y,
                                  const arma::vec& weights,
                                  const double tol = 1e-10) {
    arma::mat A = X.t();
    arma::mat X_weighted = X.each_row() % weights.t();
    arma::mat XWXt = X_weighted * A;
    arma::vec rhs = X * (weights % y);

    arma::vec coeffs;
    bool solved = arma::solve(coeffs, XWXt, rhs);
    if (!solved || !coeffs.is_finite()) {
        coeffs = safe_inverse_cpp(XWXt, tol) * rhs;
    }
    return coeffs;
}

double positive_median_cpp(arma::vec x, const double fallback = 1.0) {
    std::vector<double> positive_vals;
    positive_vals.reserve(x.n_elem);
    for (arma::uword i = 0; i < x.n_elem; ++i) {
        if (std::isfinite(x(i)) && x(i) > 0) {
            positive_vals.push_back(x(i));
        }
    }
    if (positive_vals.empty()) {
        return fallback;
    }
    arma::vec vals(positive_vals);
    double med = arma::median(vals);
    if (!std::isfinite(med) || med <= 0) {
        return fallback;
    }
    return med;
}

arma::vec wls_event_weights_cpp(const arma::vec& y,
                                const arma::vec& noise_floor,
                                const arma::vec& signal_scale,
                                const double max_weight_ratio = 1600.0) {
    const arma::uword n_detectors = y.n_elem;
    if (noise_floor.n_elem != n_detectors || signal_scale.n_elem != n_detectors) {
        Rcpp::stop("WLS noise_floor and signal_scale must match the number of detectors.");
    }

    arma::vec denom(n_detectors, arma::fill::zeros);
    for (arma::uword j = 0; j < n_detectors; ++j) {
        double floor_j = noise_floor(j);
        if (!std::isfinite(floor_j) || floor_j <= 0) {
            floor_j = 125.0;
        }
        double scale_j = signal_scale(j);
        if (!std::isfinite(scale_j) || scale_j < 0) {
            scale_j = 1.0;
        }
        double signal_j = y(j);
        if (!std::isfinite(signal_j) || signal_j < 0) {
            signal_j = 0.0;
        }
        denom(j) = floor_j + scale_j * signal_j;
    }

    double fallback = positive_median_cpp(denom, 125.0);
    arma::vec weights(n_detectors, arma::fill::zeros);
    for (arma::uword j = 0; j < n_detectors; ++j) {
        double d = denom(j);
        if (!std::isfinite(d) || d <= 0) {
            d = fallback;
        }
        weights(j) = 1.0 / d;
    }

    double med = positive_median_cpp(weights, 1.0);
    weights /= med;

    if (std::isfinite(max_weight_ratio) && max_weight_ratio > 1.0) {
        const double half_cap = std::sqrt(max_weight_ratio);
        const double lo = 1.0 / half_cap;
        const double hi = half_cap;
        for (arma::uword j = 0; j < n_detectors; ++j) {
            if (!std::isfinite(weights(j)) || weights(j) <= 0) {
                weights(j) = 1.0;
            }
            weights(j) = std::min(std::max(weights(j), lo), hi);
        }
    }

    return weights;
}

} // namespace

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
arma::mat spectreasy_wls_unmix_cpp(const arma::mat& Y,
                                   const arma::mat& M,
                                   const arma::vec& noise_floor,
                                   const arma::vec& signal_scale,
                                   const double max_weight_ratio = 1600.0,
                                   const double tol = 1e-10) {
    const arma::uword n_cells = Y.n_rows;
    const arma::uword n_markers = M.n_rows;
    arma::mat A_out(n_cells, n_markers, arma::fill::zeros);

    for (arma::uword i = 0; i < n_cells; ++i) {
        arma::vec y = Y.row(i).t();
        arma::vec weights = wls_event_weights_cpp(y, noise_floor, signal_scale, max_weight_ratio);
        arma::vec coeffs = weighted_lsq_coeffs_cpp(M, y, weights, tol);
        A_out.row(i) = coeffs.t();
    }

    return A_out;
}

// [[Rcpp::export]]
arma::mat spectreasy_unmix_best_af_cpp(const arma::mat& Y,
                                       const arma::mat& M,
                                       const arma::uvec& fluor_idx,
                                       const arma::uvec& af_idx,
                                       const std::string& method,
                                       const arma::vec& noise_floor,
                                       const arma::vec& signal_scale,
                                       const double max_weight_ratio = 1600.0,
                                       const double tol = 1e-10,
                                       const int max_outer = 500,
                                       const int max_inner = 500) {
    const arma::uword n_cells = Y.n_rows;
    const arma::uword n_detectors = Y.n_cols;
    const arma::uword n_markers = M.n_rows;
    const arma::uword n_fluors = fluor_idx.n_elem;
    const arma::uword n_af = af_idx.n_elem;
    const bool has_wls_noise = noise_floor.n_elem == n_detectors && signal_scale.n_elem == n_detectors;
    if (method == "WLS" && !has_wls_noise) {
        Rcpp::stop("WLS noise_floor and signal_scale must match the number of detectors.");
    }

    arma::mat A_out(n_cells, n_markers, arma::fill::zeros);

    // Precompute candidate model matrices. OLS also gets fixed projection
    // matrices after a condition check because candidate conditioning does not
    // vary by event for unweighted fits.
    std::vector<arma::mat> X_list(n_af);
    std::vector<arma::mat> A_list(n_af);
    std::vector<arma::mat> P_list(n_af);
    std::vector<bool> ols_candidate_ok(n_af, true);
    int skipped_ols_candidates = 0;
    for (arma::uword k = 0; k < n_af; ++k) {
        arma::mat X_k(n_fluors + 1, n_detectors);
        for (arma::uword f = 0; f < n_fluors; ++f) {
            X_k.row(f) = M.row(fluor_idx(f));
        }
        X_k.row(n_fluors) = M.row(af_idx(k));

        arma::mat A_k = X_k.t();
        X_list[k] = X_k;
        A_list[k] = A_k;

        if (method == "OLS") {
            arma::mat AtA = A_k.t() * A_k;
            if (!normal_system_is_usable_cpp(AtA, tol)) {
                ols_candidate_ok[k] = false;
                ++skipped_ols_candidates;
                continue;
            }
            arma::mat AtA_inv = safe_inverse_cpp(AtA, tol);
            P_list[k] = A_k * AtA_inv * A_k.t();
        }
    }

    if (method == "OLS" && skipped_ols_candidates == static_cast<int>(n_af)) {
        Rcpp::stop("No usable AF candidate model for OLS unmixing; all candidate matrices are singular or ill-conditioned.");
    }

    for (arma::uword i = 0; i < n_cells; ++i) {
        arma::vec y = Y.row(i).t();
        arma::vec event_weights;
        if (method == "WLS") {
            event_weights = wls_event_weights_cpp(y, noise_floor, signal_scale, max_weight_ratio);
        }

        // 1. Find the best AF band index k*. For WLS, use weighted RSS;
        // otherwise use the existing OLS projection selection.
        double min_rss = std::numeric_limits<double>::max();
        arma::uword best_k = 0;
        arma::vec best_coeffs(n_fluors + 1, arma::fill::zeros);
        bool have_best_coeffs = false;

        for (arma::uword k = 0; k < n_af; ++k) {
            double rss;
            arma::vec coeffs_k(n_fluors + 1, arma::fill::zeros);
            const arma::mat& X_k = X_list[k];
            const arma::mat& A_k = A_list[k];

            if (method == "WLS") {
                coeffs_k = weighted_lsq_coeffs_cpp(X_k, y, event_weights, tol);
                arma::vec resid = y - A_k * coeffs_k;
                rss = arma::dot(event_weights % resid, resid);
            } else if (method == "NNLS") {
                coeffs_k = nnls_lawson_hanson_cpp(A_k, y, tol, max_outer, max_inner);
                arma::vec resid = y - A_k * coeffs_k;
                rss = arma::dot(resid, resid);
            } else {
                if (!ols_candidate_ok[k]) {
                    continue;
                }
                arma::vec resid = y - P_list[k] * y;
                rss = arma::dot(resid, resid);
            }

            if (rss < min_rss) {
                min_rss = rss;
                best_k = k;
                best_coeffs = coeffs_k;
                have_best_coeffs = (method == "WLS" || method == "NNLS");
            }
        }

        // 2. Build the model matrix X for the selected best_k.
        if (min_rss == std::numeric_limits<double>::max()) {
            Rcpp::stop("No usable AF candidate model for event %d; all candidate matrices are singular or ill-conditioned.", static_cast<int>(i + 1));
        }

        const arma::mat& X_best = X_list[best_k];
        const arma::mat& A_best = A_list[best_k];

        // 3. Unmix the cell using the selected model.
        arma::vec coeffs(n_fluors + 1, arma::fill::zeros);
        if (method == "OLS") {
            arma::mat AtA = A_best.t() * A_best;
            coeffs = safe_inverse_cpp(AtA, tol) * A_best.t() * y;
        } else if (method == "WLS") {
            if (have_best_coeffs) {
                coeffs = best_coeffs;
            } else {
                arma::vec weights = event_weights;
                coeffs = weighted_lsq_coeffs_cpp(X_best, y, weights, tol);
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

    if (method == "OLS" && skipped_ols_candidates > 0) {
        Rcpp::warning("Skipped %d ill-conditioned AF candidate model(s) during OLS AF selection.", skipped_ols_candidates);
    }

    return A_out;
}
