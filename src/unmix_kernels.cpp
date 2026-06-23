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

double positive_median_std(std::vector<double>& values, const double fallback) {
    if (values.empty()) {
        return fallback;
    }

    const std::size_t n = values.size();
    const std::size_t mid = n / 2;
    std::nth_element(values.begin(), values.begin() + mid, values.end());
    double med = values[mid];
    if (n % 2 == 0) {
        std::nth_element(values.begin(), values.begin() + mid - 1, values.end());
        med = 0.5 * (med + values[mid - 1]);
    }
    if (!std::isfinite(med) || med <= 0) {
        return fallback;
    }
    return med;
}

struct WlsAfWorkspace {
    arma::vec y;
    arma::vec weights;
    arma::vec denom;
    arma::mat fluor_w_fluor_t;
    arma::vec fluor_w_y;
    arma::mat normal;
    arma::vec rhs;
    arma::vec coeffs;
    arma::vec best_coeffs;
    arma::mat selected_model;
    std::vector<double> positive_denom;
    std::vector<double> positive_weights;

    WlsAfWorkspace(const arma::uword n_detectors, const arma::uword n_fluors)
        : y(n_detectors, arma::fill::zeros),
          weights(n_detectors, arma::fill::zeros),
          denom(n_detectors, arma::fill::zeros),
          fluor_w_fluor_t(n_fluors, n_fluors, arma::fill::zeros),
          fluor_w_y(n_fluors, arma::fill::zeros),
          normal(n_fluors + 1, n_fluors + 1, arma::fill::zeros),
          rhs(n_fluors + 1, arma::fill::zeros),
          coeffs(n_fluors + 1, arma::fill::zeros),
          best_coeffs(n_fluors + 1, arma::fill::zeros),
          selected_model(n_fluors + 1, n_detectors, arma::fill::zeros) {
        positive_denom.reserve(n_detectors);
        positive_weights.reserve(n_detectors);
    }
};

void fill_wls_event_weights(const arma::vec& y,
                            const arma::vec& noise_floor,
                            const arma::vec& signal_scale,
                            const double max_weight_ratio,
                            WlsAfWorkspace& ws) {
    const arma::uword n_detectors = y.n_elem;
    ws.positive_denom.clear();
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
        ws.denom(j) = floor_j + scale_j * signal_j;
        if (std::isfinite(ws.denom(j)) && ws.denom(j) > 0) {
            ws.positive_denom.push_back(ws.denom(j));
        }
    }

    const double fallback = positive_median_std(ws.positive_denom, 125.0);
    ws.positive_weights.clear();
    for (arma::uword j = 0; j < n_detectors; ++j) {
        double d = ws.denom(j);
        if (!std::isfinite(d) || d <= 0) {
            d = fallback;
        }
        ws.weights(j) = 1.0 / d;
        if (std::isfinite(ws.weights(j)) && ws.weights(j) > 0) {
            ws.positive_weights.push_back(ws.weights(j));
        }
    }

    const double med = positive_median_std(ws.positive_weights, 1.0);
    ws.weights /= med;

    if (std::isfinite(max_weight_ratio) && max_weight_ratio > 1.0) {
        const double half_cap = std::sqrt(max_weight_ratio);
        const double lo = 1.0 / half_cap;
        const double hi = half_cap;
        for (arma::uword j = 0; j < n_detectors; ++j) {
            if (!std::isfinite(ws.weights(j)) || ws.weights(j) <= 0) {
                ws.weights(j) = 1.0;
            }
            ws.weights(j) = std::min(std::max(ws.weights(j), lo), hi);
        }
    }
}

void fill_fluor_wls_blocks(const arma::mat& F,
                           const arma::vec& y,
                           const arma::vec& weights,
                           WlsAfWorkspace& ws) {
    const arma::uword n_fluors = F.n_rows;
    const arma::uword n_detectors = F.n_cols;
    ws.fluor_w_fluor_t.zeros();
    ws.fluor_w_y.zeros();

    for (arma::uword f1 = 0; f1 < n_fluors; ++f1) {
        double rhs_sum = 0.0;
        for (arma::uword j = 0; j < n_detectors; ++j) {
            rhs_sum += F(f1, j) * weights(j) * y(j);
        }
        ws.fluor_w_y(f1) = rhs_sum;

        for (arma::uword f2 = f1; f2 < n_fluors; ++f2) {
            double normal_sum = 0.0;
            for (arma::uword j = 0; j < n_detectors; ++j) {
                normal_sum += F(f1, j) * weights(j) * F(f2, j);
            }
            ws.fluor_w_fluor_t(f1, f2) = normal_sum;
            ws.fluor_w_fluor_t(f2, f1) = normal_sum;
        }
    }
}

double fit_af_candidate_from_blocks(const arma::mat& F,
                                    const arma::rowvec& af,
                                    const arma::vec& y,
                                    const arma::vec& weights,
                                    const double tol,
                                    WlsAfWorkspace& ws) {
    const arma::uword n_fluors = F.n_rows;
    const arma::uword n_detectors = F.n_cols;
    ws.normal.zeros();
    ws.rhs.zeros();

    for (arma::uword f1 = 0; f1 < n_fluors; ++f1) {
        ws.rhs(f1) = ws.fluor_w_y(f1);
        for (arma::uword f2 = 0; f2 < n_fluors; ++f2) {
            ws.normal(f1, f2) = ws.fluor_w_fluor_t(f1, f2);
        }
    }

    double af_w_af = 0.0;
    double af_w_y = 0.0;
    for (arma::uword f = 0; f < n_fluors; ++f) {
        double fluor_w_af = 0.0;
        for (arma::uword j = 0; j < n_detectors; ++j) {
            fluor_w_af += F(f, j) * weights(j) * af(j);
        }
        ws.normal(f, n_fluors) = fluor_w_af;
        ws.normal(n_fluors, f) = fluor_w_af;
    }
    for (arma::uword j = 0; j < n_detectors; ++j) {
        af_w_af += af(j) * weights(j) * af(j);
        af_w_y += af(j) * weights(j) * y(j);
    }
    ws.normal(n_fluors, n_fluors) = af_w_af;
    ws.rhs(n_fluors) = af_w_y;

    bool solved = arma::solve(ws.coeffs, ws.normal, ws.rhs);
    if (!solved || !ws.coeffs.is_finite()) {
        ws.coeffs = safe_inverse_cpp(ws.normal, tol) * ws.rhs;
    }

    double rss = 0.0;
    for (arma::uword j = 0; j < n_detectors; ++j) {
        double fitted = 0.0;
        for (arma::uword f = 0; f < n_fluors; ++f) {
            fitted += F(f, j) * ws.coeffs(f);
        }
        fitted += af(j) * ws.coeffs(n_fluors);
        const double resid = y(j) - fitted;
        rss += weights(j) * resid * resid;
    }
    return rss;
}

void fill_selected_model(const arma::mat& F,
                         const arma::rowvec& af,
                         WlsAfWorkspace& ws) {
    const arma::uword n_fluors = F.n_rows;
    for (arma::uword f = 0; f < n_fluors; ++f) {
        ws.selected_model.row(f) = F.row(f);
    }
    ws.selected_model.row(n_fluors) = af;
}

void unmix_wls_best_af_event(const arma::mat& Y,
                             const arma::mat& F,
                             const arma::mat& AF,
                             const arma::uvec& fluor_idx,
                             const arma::uvec& af_idx,
                             const arma::vec& noise_floor,
                             const arma::vec& signal_scale,
                             const double max_weight_ratio,
                             const double tol,
                             const bool robust,
                             const int rwls_max_iter,
                             arma::mat& A_out,
                             const arma::uword i,
                             WlsAfWorkspace& ws) {
    const arma::uword n_detectors = Y.n_cols;
    const arma::uword n_fluors = F.n_rows;
    const arma::uword n_af = AF.n_rows;
    for (arma::uword j = 0; j < n_detectors; ++j) {
        ws.y(j) = Y(i, j);
    }

    fill_wls_event_weights(ws.y, noise_floor, signal_scale, max_weight_ratio, ws);
    fill_fluor_wls_blocks(F, ws.y, ws.weights, ws);

    double min_rss = std::numeric_limits<double>::max();
    arma::uword best_k = 0;
    for (arma::uword k = 0; k < n_af; ++k) {
        const double rss = fit_af_candidate_from_blocks(
            F, AF.row(k), ws.y, ws.weights, tol, ws
        );
        if (rss < min_rss) {
            min_rss = rss;
            best_k = k;
            ws.best_coeffs = ws.coeffs;
        }
    }

    if (robust) {
        fill_selected_model(F, AF.row(best_k), ws);
        ws.best_coeffs = robust_weighted_lsq_coeffs_cpp(
            ws.selected_model, ws.y, ws.weights, 1.345, rwls_max_iter, 1e-6, tol
        );
    }

    for (arma::uword f = 0; f < n_fluors; ++f) {
        A_out(i, fluor_idx(f)) = ws.best_coeffs(f);
    }
    A_out(i, af_idx(best_k)) = ws.best_coeffs(n_fluors);
}

struct WlsBestAfWorker : public RcppParallel::Worker {
    const arma::mat& Y;
    const arma::mat& F;
    const arma::mat& AF;
    const arma::uvec& fluor_idx;
    const arma::uvec& af_idx;
    const arma::vec& noise_floor;
    const arma::vec& signal_scale;
    const double max_weight_ratio;
    const double tol;
    const bool robust;
    const int rwls_max_iter;
    arma::mat& A_out;

    WlsBestAfWorker(const arma::mat& Y,
                    const arma::mat& F,
                    const arma::mat& AF,
                    const arma::uvec& fluor_idx,
                    const arma::uvec& af_idx,
                    const arma::vec& noise_floor,
                    const arma::vec& signal_scale,
                    const double max_weight_ratio,
                    const double tol,
                    const bool robust,
                    const int rwls_max_iter,
                    arma::mat& A_out)
        : Y(Y),
          F(F),
          AF(AF),
          fluor_idx(fluor_idx),
          af_idx(af_idx),
          noise_floor(noise_floor),
          signal_scale(signal_scale),
          max_weight_ratio(max_weight_ratio),
          tol(tol),
          robust(robust),
          rwls_max_iter(rwls_max_iter),
          A_out(A_out) {}

    void operator()(std::size_t begin, std::size_t end) {
        WlsAfWorkspace ws(Y.n_cols, F.n_rows);
        for (std::size_t i = begin; i < end; ++i) {
            unmix_wls_best_af_event(
                Y, F, AF, fluor_idx, af_idx, noise_floor, signal_scale,
                max_weight_ratio, tol, robust, rwls_max_iter, A_out,
                static_cast<arma::uword>(i), ws
            );
        }
    }
};

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

    if ((method == "WLS" || method == "RWLS") && n_af > 0) {
        const arma::mat F = M.rows(fluor_idx);
        const arma::mat AF = M.rows(af_idx);
        const bool robust = method == "RWLS";
        const int threads = std::max(1, n_threads);
        if (threads > 1 && n_cells > 1) {
            WlsBestAfWorker worker(
                Y, F, AF, fluor_idx, af_idx, noise_floor, signal_scale,
                max_weight_ratio, tol, robust, rwls_max_iter, A_out
            );
            RcppParallel::parallelFor(0, n_cells, worker, 128, threads);
        } else {
            WlsAfWorkspace ws(n_detectors, n_fluors);
            for (arma::uword i = 0; i < n_cells; ++i) {
                unmix_wls_best_af_event(
                    Y, F, AF, fluor_idx, af_idx, noise_floor, signal_scale,
                    max_weight_ratio, tol, robust, rwls_max_iter, A_out, i, ws
                );
            }
        }
        return A_out;
    }

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
        if (method == "WLS" || method == "RWLS") {
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

            if (method == "WLS" || method == "RWLS") {
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
        } else if (method == "RWLS") {
            coeffs = robust_weighted_lsq_coeffs_cpp(X_best, y, event_weights, 1.345, rwls_max_iter, 1e-6, tol);
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
