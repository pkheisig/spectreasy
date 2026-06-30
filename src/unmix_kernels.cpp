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

arma::vec robust_weighted_lsq_coeffs_cpp(const arma::mat& X,
                                         const arma::vec& y,
                                         const arma::vec& base_weights,
                                         const double huber_k = 1.345,
                                         const int max_iter = 1,
                                         const double robust_tol = 1e-6,
                                         const double tol = 1e-10) {
    arma::vec weights = base_weights;
    arma::vec coeffs = weighted_lsq_coeffs_cpp(X, y, weights, tol);
    const arma::mat A = X.t();
    const double k = (std::isfinite(huber_k) && huber_k > 0) ? huber_k : 1.345;
    const int n_iter = (max_iter > 0) ? max_iter : 1;
    const double conv_tol = (std::isfinite(robust_tol) && robust_tol > 0) ? robust_tol : 1e-6;

    for (int iter = 0; iter < n_iter; ++iter) {
        arma::vec resid = y - A * coeffs;
        arma::vec std_resid = resid % arma::sqrt(arma::clamp(base_weights, 1e-12, arma::datum::inf));
        double center = arma::median(std_resid);
        arma::vec abs_dev = arma::abs(std_resid - center);
        double scale = 1.4826 * arma::median(abs_dev);
        if (!std::isfinite(scale) || scale <= tol) {
            break;
        }

        arma::vec robust_weights(std_resid.n_elem, arma::fill::ones);
        for (arma::uword j = 0; j < std_resid.n_elem; ++j) {
            double z = std::abs((std_resid(j) - center) / scale);
            if (!std::isfinite(z)) {
                robust_weights(j) = 1.0;
            } else if (z > k) {
                robust_weights(j) = k / z;
            }
        }

        weights = base_weights % robust_weights;
        arma::vec next_coeffs = weighted_lsq_coeffs_cpp(X, y, weights, tol);
        double denom = std::max(1.0, arma::norm(coeffs, 2));
        double delta = arma::norm(next_coeffs - coeffs, 2) / denom;
        coeffs = next_coeffs;
        if (!std::isfinite(delta) || delta < conv_tol) {
            break;
        }
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
arma::mat spectreasy_rwls_unmix_cpp(const arma::mat& Y,
                                    const arma::mat& M,
                                    const arma::vec& noise_floor,
                                    const arma::vec& signal_scale,
                                    const double max_weight_ratio = 1600.0,
                                    const double huber_k = 1.345,
                                    const int max_iter = 1,
                                    const double robust_tol = 1e-6,
                                    const double tol = 1e-10) {
    const arma::uword n_cells = Y.n_rows;
    const arma::uword n_markers = M.n_rows;
    arma::mat A_out(n_cells, n_markers, arma::fill::zeros);

    for (arma::uword i = 0; i < n_cells; ++i) {
        arma::vec y = Y.row(i).t();
        arma::vec weights = wls_event_weights_cpp(y, noise_floor, signal_scale, max_weight_ratio);
        arma::vec coeffs = robust_weighted_lsq_coeffs_cpp(
            M, y, weights, huber_k, max_iter, robust_tol, tol
        );
        A_out.row(i) = coeffs.t();
    }

    return A_out;
}

// [[Rcpp::export]]
Rcpp::IntegerVector spectreasy_assign_af_projection_cpp(const arma::mat& Y,
                                                        const arma::mat& F,
                                                        const arma::mat& AF,
                                                        const double tol = 1e-10) {
    const arma::uword n_cells = Y.n_rows;
    const arma::uword n_detectors = Y.n_cols;
    const arma::uword n_fluors = F.n_rows;
    const arma::uword n_af = AF.n_rows;
    if (n_cells == 0 || n_detectors == 0 || n_fluors == 0 || n_af == 0) {
        Rcpp::stop("Projection AF assignment requires non-empty Y, fluorophore matrix, and AF matrix.");
    }
    if (F.n_cols != n_detectors || AF.n_cols != n_detectors) {
        Rcpp::stop("Y, F, and AF must have the same detector columns.");
    }

    const double eps = (std::isfinite(tol) && tol > 0) ? tol : 1e-10;
    const arma::mat U = safe_inverse_cpp(F * F.t(), eps) * F;
    const arma::mat V = U * AF.t();
    const arma::mat Rlib = AF.t() - F.t() * V;
    arma::rowvec denominator = arma::sum(arma::square(Rlib), 0);
    for (arma::uword k = 0; k < n_af; ++k) {
        if (!std::isfinite(denominator(k)) || denominator(k) <= eps) {
            denominator(k) = eps;
        }
    }

    const arma::mat F0 = Y * U.t();
    Rcpp::IntegerVector assignments(n_cells);
    for (arma::uword i = 0; i < n_cells; ++i) {
        double best_score = arma::datum::inf;
        arma::uword best_k = 0;
        for (arma::uword k = 0; k < n_af; ++k) {
            const double af_scale = arma::dot(Y.row(i), Rlib.col(k).t()) / denominator(k);
            double score = 0.0;
            for (arma::uword f = 0; f < n_fluors; ++f) {
                score += std::abs(F0(i, f) - af_scale * V(f, k));
            }
            if (score < best_score) {
                best_score = score;
                best_k = k;
            }
        }
        assignments(i) = static_cast<int>(best_k + 1);
    }

    return assignments;
}
