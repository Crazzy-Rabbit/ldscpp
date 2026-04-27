#pragma once

#include "common.hpp"

namespace ldsc_cpp {

double normal_pdf(double x);
double norm_ppf(double p);
double h2_obs_to_liab(double h2_obs, double p, double k);
std::pair<double, double> p_z_norm(double est, double se);
double gencov_obs_to_liab(double gencov_obs, std::optional<double> p1, std::optional<double> p2,
                          std::optional<double> k1, std::optional<double> k2);
RegressionCore run_ld_score_regression_core(const std::vector<double>& y, const Matrix& x_raw,
                                            const std::vector<double>& w_ld,
                                            const std::vector<double>& n_core,
                                            const std::vector<double>& m_annot, int n_blocks,
                                            std::optional<double> intercept,
                                            std::optional<double> twostep,
                                            bool old_weights,
                                            const RegressionOps& ops,
                                            std::optional<std::vector<bool>> step1_mask = std::nullopt);
LdscRegressionResult finish_regression(const RegressionCore& core, const std::vector<double>& m_annot);
HsqResult run_hsq(const std::vector<double>& chisq, const Matrix& ld, const std::vector<double>& w_ld,
                  const std::vector<double>& n, const std::vector<double>& m, int n_blocks,
                  std::optional<double> intercept, std::optional<double> twostep, bool old_weights);
GencovResult run_gencov(const std::vector<double>& z1, const std::vector<double>& z2, const Matrix& ld,
                        const std::vector<double>& w_ld, const std::vector<double>& n1,
                        const std::vector<double>& n2, const std::vector<double>& m,
                        double hsq1, double hsq2, double intercept_hsq1, double intercept_hsq2,
                        int n_blocks, std::optional<double> intercept_gencov,
                        std::optional<double> twostep);
RGResult run_rg_model(const std::vector<double>& z1, const std::vector<double>& z2, const Matrix& ld,
                      const std::vector<double>& w_ld, const std::vector<double>& n1,
                      const std::vector<double>& n2, const std::vector<double>& m,
                      std::optional<double> intercept_hsq1, std::optional<double> intercept_hsq2,
                      std::optional<double> intercept_gencov, int n_blocks,
                      std::optional<double> twostep);
std::string hsq_summary(const HsqResult& r, const std::vector<std::string>& names,
                        std::optional<double> sample_prev = std::nullopt,
                        std::optional<double> pop_prev = std::nullopt,
                        bool overlap = false);
std::string gencov_summary(const GencovResult& r, const std::vector<std::string>& names,
                           const std::vector<std::optional<double>>& sample_prev = {},
                           const std::vector<std::optional<double>>& pop_prev = {});
std::string rg_summary(const RGResult& r, bool silly = false);

}  // namespace ldsc_cpp
