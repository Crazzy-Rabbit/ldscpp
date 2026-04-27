#include "regression.hpp"

namespace ldsc_cpp {

RegressionCore run_ld_score_regression_core(const std::vector<double>& y, const Matrix& x_raw,
                                            const std::vector<double>& w_ld,
                                            const std::vector<double>& n_core,
                                            const std::vector<double>& m_annot, int n_blocks,
                                            std::optional<double> intercept,
                                            std::optional<double> twostep,
                                            bool old_weights,
                                            const RegressionOps& ops,
                                            std::optional<std::vector<bool>> step1_mask) {
    const size_t n = y.size();
    const size_t p = x_raw.cols;
    if (x_raw.rows != n || w_ld.size() != n || n_core.size() != n || m_annot.size() != p) {
        throw Error("LD Score regression dimension mismatch.");
    }
    const double m_tot = std::accumulate(m_annot.begin(), m_annot.end(), 0.0);
    const std::vector<double> x_tot = row_sums(x_raw);
    const double intercept_for_agg = intercept.value_or(ops.null_intercept);
    double denom = 0.0;
    for (size_t i = 0; i < n; ++i) denom += x_tot[i] * n_core[i];
    denom /= static_cast<double>(n);
    const double aggregate = m_tot * (mean(y) - intercept_for_agg) / denom;
    std::vector<double> initial_w = ops.initial_weights(x_tot, w_ld, n_core, m_tot, aggregate, intercept);
    const double nbar = mean(n_core);

    Matrix x_scaled(n, p);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < p; ++j) x_scaled(i, j) = n_core[i] * x_raw(i, j) / nbar;
    }

    Matrix x_reg = x_scaled;
    Matrix x_tot_matrix(n, 1);
    for (size_t i = 0; i < n; ++i) x_tot_matrix(i, 0) = x_tot[i];
    std::vector<double> yp = y;
    if (!intercept) {
        x_reg = append_intercept(x_scaled);
        x_tot_matrix = append_intercept(x_tot_matrix);
    } else {
        for (double& v : yp) v -= *intercept;
    }

    RegressionCore out;
    out.nbar = nbar;
    out.constrained = intercept.has_value();
    out.constrained_intercept = intercept;

    if (twostep && intercept) throw Error("twostep is not compatible with constrained intercept.");
    if (twostep && p > 1) throw Error("twostep not compatible with partitioned LD Score yet.");

    if (twostep) {
        std::vector<bool> mask = step1_mask.value_or(std::vector<bool>(n, true));
        std::vector<int> rows;
        for (size_t i = 0; i < n; ++i) {
            if (mask[i]) rows.push_back(static_cast<int>(i));
        }
        if (rows.empty()) throw Error("No SNPs retained for two-step estimator step 1.");
        Matrix x1 = x_reg.select_rows(rows);
        std::vector<double> yp1 = select_vec(yp, rows);
        std::vector<double> w1 = select_vec(w_ld, rows);
        std::vector<double> n1 = select_vec(n_core, rows);
        std::vector<double> initial_w1 = select_vec(initial_w, rows);
        auto update1 = [&](const std::vector<double>& coef, const Matrix& ref, const std::vector<double>& wld,
                           const std::vector<double>& ncore, double mtot, double nb, std::optional<double> intc,
                           const std::vector<int>* src_rows) {
            return ops.update_weights(coef, ref, wld, ncore, mtot, nb, intc, src_rows);
        };
        JackknifeResult step1 = irwls(x1, yp1, update1, n_blocks, initial_w1, std::nullopt, &x1, &w1, &n1,
                                      m_tot, nbar, std::nullopt, &rows);
        const double step1_int = step1.est[p];

        std::vector<double> yp2 = y;
        for (double& v : yp2) v -= step1_int;
        Matrix x2 = x_reg.first_cols(p);
        Matrix x_tot2 = x_tot_matrix.first_cols(1);
        std::vector<int> separators = update_separators(step1.separators, mask);
        JackknifeResult step2 = irwls(x2, yp2, ops.update_weights, n_blocks, initial_w, separators, &x_tot2, &w_ld,
                                      &n_core, m_tot, nbar, step1_int, nullptr);
        double c_num = 0.0;
        double c_den = 0.0;
        for (size_t i = 0; i < n; ++i) {
            c_num += initial_w[i] * x2(i, 0);
            c_den += initial_w[i] * x2(i, 0) * x2(i, 0);
        }
        const double c = c_num / c_den;
        JackknifeResult combined;
        combined.est = step2.est;
        combined.est.push_back(step1_int);
        combined.delete_values = Matrix(step1.delete_values.rows, p + 1);
        for (size_t b = 0; b < combined.delete_values.rows; ++b) {
            combined.delete_values(b, p) = step1.delete_values(b, p);
            combined.delete_values(b, 0) = step2.delete_values(b, 0) -
                                           c * (step1.delete_values(b, p) - step1_int);
        }
        Matrix pseudo = delete_values_to_pseudovalues(combined.delete_values, combined.est);
        jackknife_from_pseudovalues(pseudo, combined.jknife_est, combined.jknife_var, combined.jknife_se,
                                    combined.jknife_cov);
        combined.separators = step2.separators;
        out.jknife = std::move(combined);
    } else if (old_weights) {
        std::vector<double> sqrt_initial(initial_w.size());
        for (size_t i = 0; i < initial_w.size(); ++i) sqrt_initial[i] = std::sqrt(initial_w[i]);
        Matrix xw = weighted_matrix(x_reg, sqrt_initial);
        std::vector<double> yw = weighted_vec(yp, sqrt_initial);
        out.jknife = lstsq_jackknife_fast(xw, yw, n_blocks);
    } else {
        out.jknife = irwls(x_reg, yp, ops.update_weights, n_blocks, initial_w, std::nullopt, &x_tot_matrix, &w_ld,
                           &n_core, m_tot, nbar, intercept, nullptr);
    }
    return out;
}

LdscRegressionResult finish_regression(const RegressionCore& core, const std::vector<double>& m_annot) {
    const size_t p = m_annot.size();
    LdscRegressionResult out;
    out.n_annot = p;
    out.constrained = core.constrained;
    out.constrained_intercept = core.constrained_intercept;
    out.jknife = core.jknife;
    out.m = m_annot;
    out.coef.assign(p, 0.0);
    out.coef_cov = Matrix(p, p);
    out.coef_se.assign(p, 0.0);
    for (size_t j = 0; j < p; ++j) {
        out.coef[j] = out.jknife.est[j] / core.nbar;
        for (size_t k = 0; k < p; ++k) out.coef_cov(j, k) = out.jknife.jknife_cov(j, k) / (core.nbar * core.nbar);
        out.coef_se[j] = std::sqrt(out.coef_cov(j, j));
    }
    out.cat.assign(p, 0.0);
    out.cat_cov = Matrix(p, p);
    out.cat_se.assign(p, 0.0);
    for (size_t j = 0; j < p; ++j) {
        out.cat[j] = m_annot[j] * out.coef[j];
        for (size_t k = 0; k < p; ++k) out.cat_cov(j, k) = m_annot[j] * m_annot[k] * out.coef_cov(j, k);
        out.cat_se[j] = std::sqrt(out.cat_cov(j, j));
    }
    out.tot = std::accumulate(out.cat.begin(), out.cat.end(), 0.0);
    out.tot_cov = 0.0;
    for (size_t j = 0; j < p; ++j) {
        for (size_t k = 0; k < p; ++k) out.tot_cov += out.cat_cov(j, k);
    }
    out.tot_se = std::sqrt(out.tot_cov);

    Matrix numer(out.jknife.delete_values.rows, p);
    Matrix denom(out.jknife.delete_values.rows, p);
    for (size_t b = 0; b < numer.rows; ++b) {
        double d = 0.0;
        for (size_t j = 0; j < p; ++j) {
            numer(b, j) = m_annot[j] * out.jknife.delete_values(b, j) / core.nbar;
            d += numer(b, j);
        }
        for (size_t j = 0; j < p; ++j) denom(b, j) = d;
    }
    std::vector<double> prop_est(p);
    for (size_t j = 0; j < p; ++j) prop_est[j] = out.cat[j] / out.tot;
    RatioJackknifeResult prop_jk = ratio_jackknife(prop_est, numer, denom);
    out.prop = prop_jk.est;
    out.prop_cov = prop_jk.jknife_cov;
    out.prop_se = prop_jk.jknife_se;

    const double m_tot = std::accumulate(m_annot.begin(), m_annot.end(), 0.0);
    out.enrichment.assign(p, 0.0);
    out.m_prop.assign(p, 0.0);
    for (size_t j = 0; j < p; ++j) {
        out.m_prop[j] = m_annot[j] / m_tot;
        out.enrichment[j] = (out.cat[j] / m_annot[j]) / (out.tot / m_tot);
    }
    if (!core.constrained) {
        out.intercept = out.jknife.est[p];
        out.intercept_se = out.jknife.jknife_se[p];
    } else {
        out.intercept = *core.constrained_intercept;
    }
    out.tot_delete_values = Matrix(out.jknife.delete_values.rows, 1);
    out.part_delete_values = Matrix(out.jknife.delete_values.rows, p);
    out.coef_delete_values = Matrix(out.jknife.delete_values.rows, p);
    for (size_t b = 0; b < out.jknife.delete_values.rows; ++b) {
        double total = 0.0;
        for (size_t j = 0; j < p; ++j) {
            const double part = out.jknife.delete_values(b, j) / core.nbar;
            out.part_delete_values(b, j) = part;
            out.coef_delete_values(b, j) = part;
            total += m_annot[j] * part;
        }
        out.tot_delete_values(b, 0) = total;
    }
    return out;
}

HsqResult run_hsq(const std::vector<double>& chisq, const Matrix& ld, const std::vector<double>& w_ld,
                  const std::vector<double>& n, const std::vector<double>& m, int n_blocks,
                  std::optional<double> intercept, std::optional<double> twostep, bool old_weights) {
    RegressionOps ops;
    ops.null_intercept = 1.0;
    ops.initial_weights = [](const std::vector<double>& ld_tot, const std::vector<double>& wld,
                             const std::vector<double>& ncore, double mtot, double aggregate,
                             std::optional<double> intc) {
        return hsq_weights(ld_tot, wld, ncore, mtot, aggregate, intc.value_or(1.0));
    };
    ops.update_weights = [](const std::vector<double>& coef, const Matrix& ref, const std::vector<double>& wld,
                            const std::vector<double>& ncore, double mtot, double nbar, std::optional<double> intc,
                            const std::vector<int>*) {
        const double hsq = mtot * coef[0] / nbar;
        const double intercept_value = intc ? *intc : coef[1];
        return hsq_weights(ref.col(0), wld, ncore, mtot, hsq, intercept_value);
    };
    std::optional<std::vector<bool>> mask;
    if (twostep) {
        std::vector<bool> msk(chisq.size());
        for (size_t i = 0; i < chisq.size(); ++i) msk[i] = chisq[i] < *twostep;
        mask = std::move(msk);
    }
    RegressionCore core = run_ld_score_regression_core(chisq, ld, w_ld, n, m, n_blocks, intercept, twostep,
                                                       old_weights, ops, mask);
    HsqResult out;
    static_cast<LdscRegressionResult&>(out) = finish_regression(core, m);
    out.mean_chisq = mean(chisq);
    out.lambda_gc = median(chisq) / 0.4549;
    if (!out.constrained) {
        if (out.mean_chisq > 1.0) {
            out.ratio = (out.intercept - 1.0) / (out.mean_chisq - 1.0);
            out.ratio_se = out.intercept_se / (out.mean_chisq - 1.0);
        } else {
            out.ratio_na = true;
        }
    }
    return out;
}

double normal_pdf(double x) {
    static const double inv_sqrt_2pi = 0.39894228040143267794;
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

double norm_ppf(double p) {
    if (!(p > 0.0 && p < 1.0)) throw Error("Probability must be in (0,1).");
    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                               1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                               6.680131188771972e+01, -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                               -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                               3.754408661907416e+00};
    const double plow = 0.02425;
    const double phigh = 1.0 - plow;
    double q, r;
    if (p < plow) {
        q = std::sqrt(-2.0 * std::log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    if (p > phigh) {
        q = std::sqrt(-2.0 * std::log(1.0 - p));
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    q = p - 0.5;
    r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
}

double h2_obs_to_liab(double h2_obs, double p, double k) {
    if (!(k > 0.0 && k < 1.0)) throw Error("K must be in the range (0,1)");
    if (!(p > 0.0 && p < 1.0)) throw Error("P must be in the range (0,1)");
    const double thresh = norm_ppf(1.0 - k);
    const double cf = k * k * (1.0 - k) * (1.0 - k) / (p * (1.0 - p) * std::pow(normal_pdf(thresh), 2.0));
    return h2_obs * cf;
}

std::string hsq_summary(const HsqResult& r, const std::vector<std::string>& names,
                        std::optional<double> sample_prev,
                        std::optional<double> pop_prev,
                        bool overlap) {
    const bool liability = sample_prev && pop_prev;
    const std::string scale = liability ? "Liability" : "Observed";
    const double c = liability ? h2_obs_to_liab(1.0, *sample_prev, *pop_prev) : 1.0;
    std::ostringstream out;
    out << "Total " << scale << " scale h2: " << format_num(c * r.tot) << " (" << format_num(c * r.tot_se) << ")\n";
    if (r.n_annot > 1) {
        out << "Categories: ";
        for (size_t i = 0; i < names.size(); ++i) {
            if (i) out << ' ';
            out << names[i];
        }
        out << "\n";
        if (!overlap) {
            std::vector<double> cat = r.cat, cat_se = r.cat_se;
            for (double& x : cat) x *= c;
            for (double& x : cat_se) x *= c;
            out << scale << " scale h2: " << join_nums(cat) << "\n";
            out << scale << " scale h2 SE: " << join_nums(cat_se) << "\n";
            out << "Proportion of SNPs: " << join_nums(r.m_prop) << "\n";
            out << "Proportion of h2g: " << join_nums(r.prop) << "\n";
            out << "Enrichment: " << join_nums(r.enrichment) << "\n";
            out << "Coefficients: " << join_nums(r.coef) << "\n";
            out << "Coefficient SE: " << join_nums(r.coef_se) << "\n";
        }
    }
    out << "Lambda GC: " << format_num(r.lambda_gc) << "\n";
    out << "Mean Chi^2: " << format_num(r.mean_chisq) << "\n";
    if (r.constrained) {
        out << "Intercept: constrained to " << format_num(r.intercept);
    } else {
        out << "Intercept: " << format_num(r.intercept) << " (" << format_num(r.intercept_se) << ")\n";
        if (r.mean_chisq > 1.0) {
            if (r.ratio < 0.0) {
                out << "Ratio < 0 (usually indicates GC correction).";
            } else {
                out << "Ratio: " << format_num(r.ratio) << " (" << format_num(r.ratio_se) << ")";
            }
        } else {
            out << "Ratio: NA (mean chi^2 < 1)";
        }
    }
    return out.str();
}

std::pair<double, double> p_z_norm(double est, double se) {
    double z = std::numeric_limits<double>::infinity();
    if (se != 0.0) z = est / se;
    const double p = std::erfc(std::abs(z) / std::sqrt(2.0));
    return {p, z};
}

GencovResult run_gencov(const std::vector<double>& z1, const std::vector<double>& z2, const Matrix& ld,
                        const std::vector<double>& w_ld, const std::vector<double>& n1,
                        const std::vector<double>& n2, const std::vector<double>& m,
                        double hsq1, double hsq2, double intercept_hsq1, double intercept_hsq2,
                        int n_blocks, std::optional<double> intercept_gencov,
                        std::optional<double> twostep) {
    std::vector<double> y(z1.size());
    std::vector<double> sqrt_n(z1.size());
    for (size_t i = 0; i < z1.size(); ++i) {
        y[i] = z1[i] * z2[i];
        sqrt_n[i] = std::sqrt(n1[i] * n2[i]);
    }
    RegressionOps ops;
    ops.null_intercept = 0.0;
    ops.initial_weights = [&](const std::vector<double>& ld_tot, const std::vector<double>& wld,
                              const std::vector<double>&, double mtot, double aggregate,
                              std::optional<double> intc) {
        return gencov_weights(ld_tot, wld, n1, n2, mtot, hsq1, hsq2, aggregate, intc, intercept_hsq1, intercept_hsq2);
    };
    ops.update_weights = [&](const std::vector<double>& coef, const Matrix& ref, const std::vector<double>& wld,
                             const std::vector<double>&, double mtot, double nbar, std::optional<double> intc,
                             const std::vector<int>* source_rows) {
        const double rho_g = mtot * coef[0] / nbar;
        const double intercept_value = intc ? *intc : coef[1];
        std::vector<double> nn1 = n1;
        std::vector<double> nn2 = n2;
        if (source_rows) {
            nn1 = select_vec(n1, *source_rows);
            nn2 = select_vec(n2, *source_rows);
        }
        return gencov_weights(ref.col(0), wld, nn1, nn2, mtot, hsq1, hsq2, rho_g, intercept_value,
                              intercept_hsq1, intercept_hsq2);
    };
    std::optional<std::vector<bool>> mask;
    if (twostep) {
        std::vector<bool> msk(z1.size());
        for (size_t i = 0; i < z1.size(); ++i) msk[i] = z1[i] * z1[i] < *twostep && z2[i] * z2[i] < *twostep;
        mask = std::move(msk);
    }
    RegressionCore core = run_ld_score_regression_core(y, ld, w_ld, sqrt_n, m, n_blocks, intercept_gencov,
                                                       twostep, false, ops, mask);
    GencovResult out;
    static_cast<LdscRegressionResult&>(out) = finish_regression(core, m);
    auto pz = p_z_norm(out.tot, out.tot_se);
    out.p_value = pz.first;
    out.z_score = pz.second;
    out.mean_z1z2 = mean(y);
    return out;
}

double gencov_obs_to_liab(double gencov_obs, std::optional<double> p1, std::optional<double> p2,
                          std::optional<double> k1, std::optional<double> k2) {
    double c1 = 1.0;
    double c2 = 1.0;
    if (p1 && k1) c1 = std::sqrt(h2_obs_to_liab(1.0, *p1, *k1));
    if (p2 && k2) c2 = std::sqrt(h2_obs_to_liab(1.0, *p2, *k2));
    return gencov_obs * c1 * c2;
}

std::string gencov_summary(const GencovResult& r, const std::vector<std::string>& names,
                           const std::vector<std::optional<double>>& sample_prev,
                           const std::vector<std::optional<double>>& pop_prev) {
    bool liability = sample_prev.size() == 2 && pop_prev.size() == 2 && sample_prev[0] && sample_prev[1] &&
                     pop_prev[0] && pop_prev[1];
    const std::string scale = liability ? "Liability" : "Observed";
    const double c = liability ? gencov_obs_to_liab(1.0, sample_prev[0], sample_prev[1], pop_prev[0], pop_prev[1]) : 1.0;
    std::ostringstream out;
    out << "Total " << scale << " scale gencov: " << format_num(c * r.tot) << " (" << format_num(c * r.tot_se) << ")\n";
    if (r.n_annot > 1) {
        out << "Categories: ";
        for (size_t i = 0; i < names.size(); ++i) {
            if (i) out << ' ';
            out << names[i];
        }
        out << "\n";
        std::vector<double> cat = r.cat, cat_se = r.cat_se;
        for (double& x : cat) x *= c;
        for (double& x : cat_se) x *= c;
        out << scale << " scale gencov: " << join_nums(cat) << "\n";
        out << scale << " scale gencov SE: " << join_nums(cat_se) << "\n";
        out << "Proportion of SNPs: " << join_nums(r.m_prop) << "\n";
        out << "Proportion of gencov: " << join_nums(r.prop) << "\n";
        out << "Enrichment: " << join_nums(r.enrichment) << "\n";
    }
    out << "Mean z1*z2: " << format_num(r.mean_z1z2) << "\n";
    if (r.constrained) {
        out << "Intercept: constrained to " << format_num(r.intercept);
    } else {
        out << "Intercept: " << format_num(r.intercept) << " (" << format_num(r.intercept_se) << ")";
    }
    return out.str();
}

RGResult run_rg_model(const std::vector<double>& z1, const std::vector<double>& z2, const Matrix& ld,
                      const std::vector<double>& w_ld, const std::vector<double>& n1, const std::vector<double>& n2,
                      const std::vector<double>& m, std::optional<double> intercept_hsq1,
                      std::optional<double> intercept_hsq2, std::optional<double> intercept_gencov,
                      int n_blocks, std::optional<double> twostep) {
    std::vector<double> z1sq(z1.size()), z2sq(z2.size());
    for (size_t i = 0; i < z1.size(); ++i) {
        z1sq[i] = z1[i] * z1[i];
        z2sq[i] = z2[i] * z2[i];
    }
    RGResult out;
    out.hsq1 = run_hsq(z1sq, ld, w_ld, n1, m, n_blocks, intercept_hsq1, twostep, false);
    out.hsq2 = run_hsq(z2sq, ld, w_ld, n2, m, n_blocks, intercept_hsq2, twostep, false);
    out.gencov = run_gencov(z1, z2, ld, w_ld, n1, n2, m, out.hsq1.tot, out.hsq2.tot,
                            out.hsq1.intercept, out.hsq2.intercept, n_blocks, intercept_gencov, twostep);
    if (out.hsq1.tot <= 0.0 || out.hsq2.tot <= 0.0) {
        out.negative_hsq = true;
        return out;
    }
    const double denom = std::sqrt(out.hsq1.tot * out.hsq2.tot);
    out.rg_ratio = out.gencov.tot / denom;
    Matrix denom_delete(out.gencov.tot_delete_values.rows, 1);
    for (size_t i = 0; i < denom_delete.rows; ++i) {
        denom_delete(i, 0) = std::sqrt(out.hsq1.tot_delete_values(i, 0) * out.hsq2.tot_delete_values(i, 0));
    }
    RatioJackknifeResult rg_jk = ratio_jackknife({out.rg_ratio}, out.gencov.tot_delete_values, denom_delete);
    out.rg_se = rg_jk.jknife_se[0];
    auto pz = p_z_norm(out.rg_ratio, out.rg_se);
    out.p_value = pz.first;
    out.z_score = pz.second;
    return out;
}

std::string rg_summary(const RGResult& r, bool silly) {
    std::ostringstream out;
    if (r.negative_hsq) {
        out << "Genetic Correlation: nan (nan) (h2  out of bounds) \n";
        out << "Z-score: nan (nan) (h2  out of bounds)\n";
        out << "P: nan (nan) (h2  out of bounds)\n";
        out << "WARNING: One of the h2's was out of bounds.\n";
        out << "This usually indicates a data-munging error or that h2 or N is low.";
    } else if ((r.rg_ratio > 1.2 || r.rg_ratio < -1.2) && !silly) {
        out << "Genetic Correlation: nan (nan) (rg out of bounds) \n";
        out << "Z-score: nan (nan) (rg out of bounds)\n";
        out << "P: nan (nan) (rg out of bounds)\n";
        out << "WARNING: rg was out of bounds.";
    } else {
        out << "Genetic Correlation: " << format_num(r.rg_ratio) << " (" << format_num(r.rg_se) << ")\n";
        out << "Z-score: " << format_num(r.z_score) << "\n";
        out << "P: " << (r.p_value == 0.0 ? std::string("0.") : format_num(r.p_value));
    }
    return out.str();
}

}  // namespace ldsc_cpp
