#pragma once

#include "types.hpp"

namespace ldsc_cpp {

std::string trim(const std::string& s);
std::vector<std::string> split_ws(const std::string& s);
std::vector<std::string> split_char(const std::string& s, char delim);
std::vector<std::optional<double>> split_optional_doubles(const std::optional<std::string>& value, size_t n);
std::string sec_to_str(double seconds);
std::string time_now_string();
std::string masthead();
std::string command_header(int argc, char** argv, const std::string& display_name);
std::string command_header_from_args(const std::vector<std::string>& args, const std::string& display_name);
bool ends_with(const std::string& s, const std::string& suffix);
bool is_compressed_name(const std::string& path);
std::vector<std::string> read_text_lines(const std::string& path);
void write_gzip_text(const std::string& path, const std::string& content);
std::vector<std::vector<std::string>> read_rows_no_header(const std::string& path);
double parse_double(const std::string& token);
std::string format_num(double x, int precision = 4);
std::string format_py_float(double x);
std::string join_nums(const std::vector<double>& xs);
double mean(const std::vector<double>& x);
double median(std::vector<double> x);
Matrix append_intercept(const Matrix& x);
std::vector<double> row_sums(const Matrix& x);
std::vector<double> select_vec(const std::vector<double>& x, const std::vector<int>& idx);
Matrix weighted_matrix(const Matrix& x, const std::vector<double>& w);
std::vector<double> weighted_vec(const std::vector<double>& y, const std::vector<double>& w);
std::vector<double> solve_linear(Matrix a, std::vector<double> b);
std::vector<double> least_squares(const Matrix& x, const std::vector<double>& y);
std::vector<int> get_separators(size_t n, int n_blocks);
void jackknife_from_pseudovalues(const Matrix& pseudovalues, std::vector<double>& jk_est,
                                 std::vector<double>& jk_var, std::vector<double>& jk_se, Matrix& jk_cov);
Matrix delete_values_to_pseudovalues(const Matrix& delete_values, const std::vector<double>& est);
JackknifeResult lstsq_jackknife_fast(const Matrix& x, const std::vector<double>& y, int n_blocks,
                                     std::optional<std::vector<int>> separators = std::nullopt);
RatioJackknifeResult ratio_jackknife(const std::vector<double>& est, const Matrix& numer, const Matrix& denom);
JackknifeResult irwls(const Matrix& x, const std::vector<double>& y, const UpdateFunc& update_func, int n_blocks,
                      const std::vector<double>& initial_inverse_cvf,
                      std::optional<std::vector<int>> separators = std::nullopt,
                      const Matrix* ref_ld_for_update = nullptr,
                      const std::vector<double>* w_ld_for_update = nullptr,
                      const std::vector<double>* n_core_for_update = nullptr,
                      double m_tot = 0.0, double nbar = 1.0,
                      std::optional<double> intercept = std::nullopt,
                      const std::vector<int>* source_rows = nullptr);
std::vector<int> update_separators(const std::vector<int>& masked_separators, const std::vector<bool>& mask);
Table read_table(const std::string& path);
std::string resolve_existing(const std::vector<std::string>& candidates, const std::string& label);
std::string resolve_ldscore_file(const std::string& prefix);
std::string resolve_m_file(const std::string& prefix, bool common);
std::string sub_chr(std::string pattern, int chr);
std::vector<std::string> present_chr_ldscore_files(const std::string& prefix);
std::vector<std::string> present_chr_m_files(const std::string& prefix, bool common);
LdFrame read_ldscore_paths(const std::vector<std::string>& paths);
LdFrame read_ldscore_prefix(const std::string& prefix, bool chr_split);
LdFrame read_ldscore_list(const std::string& arg, bool chr_split);
std::vector<double> read_m_line(const std::string& path);
std::vector<double> read_m_prefix(const std::string& prefix, bool chr_split, bool common);
std::vector<double> read_m(const std::optional<std::string>& manual_m, const std::string& ref_arg,
                           bool chr_split, bool common);
void remove_zero_variance(LdFrame& ref_ld, std::vector<double>& m_annot);
std::vector<SumRecord> read_sumstats(const std::string& path, bool alleles, bool dropna);
std::unordered_map<std::string, size_t> index_by_snp(const std::vector<std::string>& snps);
BaseData read_ld_sumstats(const std::string& sumstats_path, const std::string& ref_arg, bool ref_chr,
                          const std::string& w_arg, bool w_chr, const std::optional<std::string>& manual_m,
                          bool common_m, std::vector<double>& m_annot, Logger& log,
                          bool alleles = false, bool dropna = true);
std::vector<double> hsq_weights(const std::vector<double>& ld, const std::vector<double>& w_ld,
                                const std::vector<double>& n, double m, double hsq,
                                std::optional<double> intercept = std::nullopt);
std::vector<double> gencov_weights(const std::vector<double>& ld, const std::vector<double>& w_ld,
                                   const std::vector<double>& n1, const std::vector<double>& n2,
                                   double m, double h1, double h2, double rho_g,
                                   std::optional<double> intercept_gencov,
                                   std::optional<double> intercept_hsq1,
                                   std::optional<double> intercept_hsq2);
void save_matrix(const std::string& path, const Matrix& m);
void ensure_parent_dir(const std::string& out_prefix);
char complement_base(char c);
bool valid_base(char c);
bool alleles_match(const std::string& a);
bool flip_alleles(const std::string& a);

}  // namespace ldsc_cpp
