#pragma once

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace ldsc_cpp {

inline constexpr double kNan = std::numeric_limits<double>::quiet_NaN();

struct Error : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct Matrix {
    size_t rows = 0;
    size_t cols = 0;
    std::vector<double> data;

    Matrix() = default;
    Matrix(size_t r, size_t c, double v = 0.0) : rows(r), cols(c), data(r * c, v) {}

    double& operator()(size_t r, size_t c) { return data[r * cols + c]; }
    double operator()(size_t r, size_t c) const { return data[r * cols + c]; }

    std::vector<double> col(size_t c) const {
        std::vector<double> out(rows);
        for (size_t i = 0; i < rows; ++i) out[i] = (*this)(i, c);
        return out;
    }

    Matrix select_rows(const std::vector<int>& idx) const {
        Matrix out(idx.size(), cols);
        for (size_t r = 0; r < idx.size(); ++r) {
            for (size_t c = 0; c < cols; ++c) out(r, c) = (*this)(static_cast<size_t>(idx[r]), c);
        }
        return out;
    }

    Matrix first_cols(size_t n) const {
        Matrix out(rows, n);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < n; ++j) out(i, j) = (*this)(i, j);
        }
        return out;
    }
};

struct JackknifeResult {
    std::vector<double> est;
    std::vector<double> jknife_est;
    std::vector<double> jknife_var;
    std::vector<double> jknife_se;
    Matrix jknife_cov;
    Matrix delete_values;
    std::vector<int> separators;
};

struct RatioJackknifeResult {
    std::vector<double> est;
    std::vector<double> jknife_est;
    std::vector<double> jknife_var;
    std::vector<double> jknife_se;
    Matrix jknife_cov;
};

struct Table {
    std::vector<std::string> columns;
    std::unordered_map<std::string, size_t> index;
    std::vector<std::vector<std::string>> rows;

    int col(const std::string& name) const {
        auto it = index.find(name);
        return it == index.end() ? -1 : static_cast<int>(it->second);
    }
};

struct LdFrame {
    std::vector<std::string> snp;
    std::vector<std::string> names;
    Matrix values;
};

struct LdTempRow {
    size_t order = 0;
    std::string snp;
    double chr = 0.0;
    double bp = 0.0;
    std::vector<double> values;
};

struct SumRecord {
    std::string snp;
    std::string a1;
    std::string a2;
    double z = kNan;
    double n = kNan;
};

struct Logger {
    std::ofstream file;
    explicit Logger(const std::string& path) : file(path) {
        if (!file) throw Error("Could not open log file: " + path);
    }
    void log(const std::string& msg) {
        std::cout << msg << "\n";
        file << msg << "\n";
    }
};

struct BaseData {
    std::vector<std::string> snp;
    std::vector<std::string> a1;
    std::vector<std::string> a2;
    std::vector<double> z;
    std::vector<double> n;
    std::vector<double> w_ld;
    Matrix ld;
    std::vector<std::string> ld_names;
};

using UpdateFunc = std::function<std::vector<double>(const std::vector<double>& coef, const Matrix& ref_ld_tot,
                                                     const std::vector<double>& w_ld,
                                                     const std::vector<double>& n_core, double m_tot,
                                                     double nbar, std::optional<double> intercept,
                                                     const std::vector<int>* source_rows)>;

struct RegressionOps {
    double null_intercept = 1.0;
    std::function<std::vector<double>(const std::vector<double>& ld, const std::vector<double>& w_ld,
                                      const std::vector<double>& n_core, double m_tot, double aggregate,
                                      std::optional<double> intercept)>
        initial_weights;
    UpdateFunc update_weights;
};

struct RegressionCore {
    JackknifeResult jknife;
    double nbar = 1.0;
    bool constrained = false;
    std::optional<double> constrained_intercept;
};

struct LdscRegressionResult {
    size_t n_annot = 0;
    bool constrained = false;
    std::optional<double> constrained_intercept;
    double intercept = kNan;
    double intercept_se = kNan;
    std::vector<double> coef;
    Matrix coef_cov;
    std::vector<double> coef_se;
    std::vector<double> cat;
    Matrix cat_cov;
    std::vector<double> cat_se;
    double tot = kNan;
    double tot_cov = kNan;
    double tot_se = kNan;
    std::vector<double> prop;
    Matrix prop_cov;
    std::vector<double> prop_se;
    std::vector<double> enrichment;
    std::vector<double> m_prop;
    Matrix tot_delete_values;
    Matrix part_delete_values;
    Matrix coef_delete_values;
    JackknifeResult jknife;
    std::vector<double> m;
};

struct HsqResult : public LdscRegressionResult {
    double mean_chisq = kNan;
    double lambda_gc = kNan;
    double ratio = kNan;
    double ratio_se = kNan;
    bool ratio_na = false;
};

struct GencovResult : public LdscRegressionResult {
    double p_value = kNan;
    double z_score = kNan;
    double mean_z1z2 = kNan;
};

struct RGResult {
    HsqResult hsq1;
    HsqResult hsq2;
    GencovResult gencov;
    bool negative_hsq = false;
    double rg_ratio = kNan;
    double rg_se = kNan;
    double p_value = kNan;
    double z_score = kNan;
};

struct Args {
    std::string program_name;
    std::string command;
    std::vector<std::string> positionals;
    std::vector<std::string> raw_args;
    std::optional<std::string> h2;
    std::optional<std::string> h2_cts;
    std::optional<std::string> rg;
    std::optional<std::string> ref_ld;
    std::optional<std::string> ref_ld_chr;
    std::optional<std::string> w_ld;
    std::optional<std::string> w_ld_chr;
    std::optional<std::string> manual_m;
    std::optional<double> intercept_h2_scalar;
    std::optional<std::string> intercept_h2_list;
    std::optional<std::string> intercept_gencov_list;
    std::optional<double> twostep;
    std::optional<double> chisq_max;
    std::optional<std::string> samp_prev;
    std::optional<std::string> pop_prev;
    std::string out = "ldsc";
    bool out_set = false;
    int n_blocks = 200;
    bool no_intercept = false;
    bool not_m_5_50 = false;
    bool print_cov = false;
    bool print_delete_vals = false;
    bool no_check_alleles = false;
    bool return_silly_things = false;
    std::optional<std::string> ref_ld_chr_cts;
    bool print_all_cts = false;
    bool overlap_annot = false;
    bool print_coefficients = false;
    std::optional<std::string> frqfile;
    std::optional<std::string> frqfile_chr;

    std::optional<std::string> sumstats;
    std::optional<double> n_value;
    std::optional<double> n_cas;
    std::optional<double> n_con;
    double info_min = 0.9;
    double maf_min = 0.01;
    bool daner = false;
    bool daner_n = false;
    bool no_alleles = false;
    std::optional<std::string> merge_alleles;
    std::optional<double> n_min;
    std::optional<std::string> snp_col;
    std::optional<std::string> n_col;
    std::optional<std::string> n_cas_col;
    std::optional<std::string> n_con_col;
    std::optional<std::string> a1_col;
    std::optional<std::string> a2_col;
    std::optional<std::string> p_col;
    std::optional<std::string> frq_col;
    std::optional<std::string> signed_sumstats;
    std::optional<std::string> info_col;
    std::optional<std::string> info_list;
    std::optional<std::string> nstudy_col;
    std::optional<double> nstudy_min;
    std::optional<std::string> ignore_cols;
    bool a1_inc = false;
    bool keep_maf = false;

    std::optional<std::string> gene_set_file;
    std::string gene_coord_file = "ENSG_coord.txt";
    int windowsize = 0;
    std::optional<std::string> bed_file;
    bool nomerge = false;
    std::optional<std::string> bimfile;
    std::optional<std::string> annot_file;

    std::optional<std::string> bfile;
    bool l2 = false;
    std::optional<std::string> extract;
    std::optional<std::string> keep;
    std::optional<int> ld_wind_snps;
    std::optional<double> ld_wind_kb;
    std::optional<double> ld_wind_cm;
    std::optional<std::string> print_snps;
    std::optional<std::string> annot;
    bool thin_annot = false;
    std::optional<std::string> cts_bin;
    std::optional<std::string> cts_breaks;
    std::optional<std::string> cts_names;
    bool per_allele = false;
    std::optional<double> pq_exp;
    bool no_print_annot = false;
    double maf = 0.0;
    bool maf_set = false;
    int chunk_size = 50;
    int chunksize = 5000000;
    bool pickle = false;
    bool invert_anyway = false;
    bool yes_really = false;
    std::optional<std::string> quantile_annotfile;
    std::optional<std::string> quantile_resultfile;
    std::optional<std::string> quantile_outfile;
    std::optional<std::string> ref_annot_chr;
    std::optional<std::string> annot_header;
    int nb_quantile = 5;
    bool exclude0 = false;
};

struct MungeRow {
    std::string snp;
    std::string a1;
    std::string a2;
    double p = kNan;
    double signed_sumstat = kNan;
    double n = kNan;
    double n_cas = kNan;
    double n_con = kNan;
    double nstudy = kNan;
    double frq = kNan;
    double z = kNan;
};

struct Interval { std::string chr; int start = 0; int end = 0; };
struct BimRecord { std::string chr; std::string snp; double cm = 0.0; int bp = 0; std::string a1; std::string a2; };
struct GenotypeData {
    std::vector<BimRecord> records;
    std::vector<size_t> kept_indices;
    Matrix x;
    std::vector<double> freq;
    std::vector<double> maf;
};
struct RGMergedData { std::vector<double> z1; std::vector<double> z2; std::vector<double> n1; std::vector<double> n2; std::vector<double> w_ld; Matrix ld; };

}  // namespace ldsc_cpp
