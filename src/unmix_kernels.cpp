#include <RcppArmadillo.h>
#include <RcppParallel.h>

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

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

struct JointAfAssignmentPrecomp {
    arma::mat P;
    arma::mat v_library;
    arma::mat r_library;
    arma::vec r_dots;
    arma::vec fluor_weights;
};

JointAfAssignmentPrecomp joint_af_assignment_precomp_cpp(const arma::mat& F,
                                                         const arma::mat& AF,
                                                         const double tol = 1e-10) {
    const arma::uword n_fluors = F.n_rows;
    const arma::uword n_detectors = F.n_cols;
    const arma::uword n_af = AF.n_rows;

    JointAfAssignmentPrecomp out;
    out.P = safe_inverse_cpp(F * F.t(), tol) * F;
    out.v_library = out.P * AF.t();
    out.r_library = AF.t() - F.t() * out.v_library;
    out.r_dots = arma::sum(arma::square(out.r_library), 0).t();
    for (arma::uword k = 0; k < out.r_dots.n_elem; ++k) {
        if (!std::isfinite(out.r_dots(k)) || out.r_dots(k) <= 0) {
            out.r_dots(k) = 1e-10;
        }
    }

    arma::mat af_cov(n_detectors, n_detectors, arma::fill::zeros);
    if (n_af > 1) {
        arma::rowvec af_mean = arma::mean(AF, 0);
        arma::mat centered = AF.each_row() - af_mean;
        af_cov = centered.t() * centered / static_cast<double>(n_af - 1);
    }

    arma::mat fluor_cov = out.P * af_cov * out.P.t();
    out.fluor_weights.set_size(n_fluors);
    for (arma::uword f = 0; f < n_fluors; ++f) {
        double value = fluor_cov(f, f);
        out.fluor_weights(f) = std::sqrt(std::abs(value)) + 1e-8;
    }

    return out;
}

arma::uword select_joint_cov_af_cpp(const arma::vec& y,
                                    const arma::mat& F,
                                    const JointAfAssignmentPrecomp& precomp);

struct BestAfWorker : public RcppParallel::Worker {
    const arma::mat& Y;
    const arma::mat& M;
    const arma::uvec& fluor_idx;
    const arma::uvec& af_idx;
    const arma::mat& F;
    const JointAfAssignmentPrecomp& af_precomp;
    const std::vector<arma::mat>& X_list;
    const std::vector<arma::mat>& A_list;
    const arma::vec& noise_floor;
    const arma::vec& signal_scale;
    const bool use_rwls;
    const double max_weight_ratio;
    const double tol;
    const int rwls_max_iter;
    arma::mat& A_out;

    BestAfWorker(const arma::mat& Y,
                 const arma::mat& M,
                 const arma::uvec& fluor_idx,
                 const arma::uvec& af_idx,
                 const arma::mat& F,
                 const JointAfAssignmentPrecomp& af_precomp,
                 const std::vector<arma::mat>& X_list,
                 const std::vector<arma::mat>& A_list,
                 const arma::vec& noise_floor,
                 const arma::vec& signal_scale,
                 const bool use_rwls,
                 const double max_weight_ratio,
                 const double tol,
                 const int rwls_max_iter,
                 arma::mat& A_out)
        : Y(Y),
          M(M),
          fluor_idx(fluor_idx),
          af_idx(af_idx),
          F(F),
          af_precomp(af_precomp),
          X_list(X_list),
          A_list(A_list),
          noise_floor(noise_floor),
          signal_scale(signal_scale),
          use_rwls(use_rwls),
          max_weight_ratio(max_weight_ratio),
          tol(tol),
          rwls_max_iter(rwls_max_iter),
          A_out(A_out) {}

    void operator()(std::size_t begin, std::size_t end) {
        const arma::uword n_fluors = fluor_idx.n_elem;
        for (std::size_t i = begin; i < end; ++i) {
            arma::vec y = Y.row(i).t();
            arma::vec event_weights = wls_event_weights_cpp(y, noise_floor, signal_scale, max_weight_ratio);
            arma::uword best_k = select_joint_cov_af_cpp(y, F, af_precomp);
            const arma::mat& X_best = X_list[best_k];

            arma::vec coeffs = use_rwls
                ? robust_weighted_lsq_coeffs_cpp(X_best, y, event_weights, 1.345, rwls_max_iter, 1e-6, tol)
                : weighted_lsq_coeffs_cpp(X_best, y, event_weights, tol);

            for (arma::uword f = 0; f < n_fluors; ++f) {
                A_out(i, fluor_idx(f)) = coeffs(f);
            }
            A_out(i, af_idx(best_k)) = coeffs(n_fluors);
        }
    }
};

arma::uword select_joint_cov_af_cpp(const arma::vec& y,
                                    const arma::mat& F,
                                    const JointAfAssignmentPrecomp& precomp) {
    const arma::uword n_af = precomp.r_library.n_cols;
    if (n_af <= 1) {
        return 0;
    }

    arma::vec unmixed = precomp.P * y;
    arma::vec unmixed_nonneg = unmixed;
    for (arma::uword f = 0; f < unmixed_nonneg.n_elem; ++f) {
        if (!std::isfinite(unmixed_nonneg(f)) || unmixed_nonneg(f) < 0) {
            unmixed_nonneg(f) = 0.0;
        }
    }

    arma::vec base_resid = y - F.t() * unmixed_nonneg;
    double base_fluor = arma::sum(precomp.fluor_weights % arma::abs(unmixed)) + 1e-6;
    double base_resid_norm = arma::norm(base_resid, 2) + 1e-6;
    if (!std::isfinite(base_fluor) || base_fluor <= 0) {
        base_fluor = 1e-6;
    }
    if (!std::isfinite(base_resid_norm) || base_resid_norm <= 0) {
        base_resid_norm = 1e-6;
    }

    double best_score = std::numeric_limits<double>::infinity();
    arma::uword best_k = 0;
    for (arma::uword k = 0; k < n_af; ++k) {
        double af_scale = arma::dot(y, precomp.r_library.col(k)) / precomp.r_dots(k);
        if (!std::isfinite(af_scale) || af_scale < 0) {
            af_scale = 0.0;
        }

        arma::vec adjusted_fluors = unmixed - af_scale * precomp.v_library.col(k);
        arma::vec adjusted_resid = base_resid - af_scale * precomp.r_library.col(k);

        double p_fluor = arma::sum(precomp.fluor_weights % arma::abs(adjusted_fluors)) / base_fluor;
        double p_resid = arma::norm(adjusted_resid, 2) / base_resid_norm;
        double score = p_fluor * p_resid;
        if (!std::isfinite(score)) {
            score = std::numeric_limits<double>::infinity();
        }

        if (score < best_score) {
            best_score = score;
            best_k = k;
        }
    }

    return best_k;
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
                                       const int max_inner = 500,
                                       const int rwls_max_iter = 1,
                                       const int n_threads = 1) {
    const arma::uword n_cells = Y.n_rows;
    const arma::uword n_detectors = Y.n_cols;
    const arma::uword n_markers = M.n_rows;
    const arma::uword n_fluors = fluor_idx.n_elem;
    const arma::uword n_af = af_idx.n_elem;
    const bool has_wls_noise = noise_floor.n_elem == n_detectors && signal_scale.n_elem == n_detectors;
    if ((method == "WLS" || method == "RWLS") && !has_wls_noise) {
        Rcpp::stop("WLS/RWLS noise_floor and signal_scale must match the number of detectors.");
    }

    arma::mat A_out(n_cells, n_markers, arma::fill::zeros);

    const arma::mat F = M.rows(fluor_idx);
    const arma::mat AF = M.rows(af_idx);
    JointAfAssignmentPrecomp af_precomp = joint_af_assignment_precomp_cpp(F, AF, tol);

    // Precompute candidate model matrices. AF selection is handled separately
    // by the joint covariance + residual score; these matrices are
    // only used for the final coefficient fit after one AF row is selected.
    std::vector<arma::mat> X_list(n_af);
    std::vector<arma::mat> A_list(n_af);
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
        }
    }

    if (method == "OLS" && skipped_ols_candidates == static_cast<int>(n_af)) {
        Rcpp::stop("No usable AF candidate model for OLS unmixing; all candidate matrices are singular or ill-conditioned.");
    }

    if ((method == "WLS" || method == "RWLS") && n_threads > 1) {
        BestAfWorker worker(
            Y,
            M,
            fluor_idx,
            af_idx,
            F,
            af_precomp,
            X_list,
            A_list,
            noise_floor,
            signal_scale,
            method == "RWLS",
            max_weight_ratio,
            tol,
            rwls_max_iter,
            A_out
        );
        RcppParallel::parallelFor(0, static_cast<std::size_t>(n_cells), worker, 100, n_threads);
        return A_out;
    }

    for (arma::uword i = 0; i < n_cells; ++i) {
        arma::vec y = Y.row(i).t();
        arma::vec event_weights;
        if (method == "WLS" || method == "RWLS") {
            event_weights = wls_event_weights_cpp(y, noise_floor, signal_scale, max_weight_ratio);
        }

        // 1. Select one AF band with the joint covariance + residual
        // assignment score, independent of the final coefficient solver.
        arma::uword best_k = select_joint_cov_af_cpp(y, F, af_precomp);
        if (method == "OLS" && !ols_candidate_ok[best_k]) {
            Rcpp::stop("The selected AF candidate model for event %d is singular or ill-conditioned.", static_cast<int>(i + 1));
        }

        // 2. Unmix the cell using the selected model.
        const arma::mat& X_best = X_list[best_k];
        const arma::mat& A_best = A_list[best_k];

        arma::vec coeffs(n_fluors + 1, arma::fill::zeros);
        if (method == "OLS") {
            arma::mat AtA = A_best.t() * A_best;
            coeffs = safe_inverse_cpp(AtA, tol) * A_best.t() * y;
        } else if (method == "WLS") {
            coeffs = weighted_lsq_coeffs_cpp(X_best, y, event_weights, tol);
        } else if (method == "RWLS") {
            coeffs = robust_weighted_lsq_coeffs_cpp(X_best, y, event_weights, 1.345, rwls_max_iter, 1e-6, tol);
        } else if (method == "NNLS") {
            coeffs = nnls_lawson_hanson_cpp(A_best, y, tol, max_outer, max_inner);
        }

        // 3. Map the resulting coefficients back to the full reference matrix layout.
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
