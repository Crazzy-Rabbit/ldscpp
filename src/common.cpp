#include "common.hpp"

#include <bzlib.h>
#include <zlib.h>

#include <cstdio>

namespace ldsc_cpp {

std::string trim(const std::string& s) {
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

std::vector<std::string> split_ws(const std::string& s) {
    std::istringstream in(s);
    std::vector<std::string> out;
    std::string tok;
    while (in >> tok) out.push_back(tok);
    return out;
}

std::vector<std::string> split_char(const std::string& s, char delim) {
    std::vector<std::string> out;
    std::string tok;
    std::stringstream ss(s);
    while (std::getline(ss, tok, delim)) out.push_back(trim(tok));
    return out;
}

std::vector<std::optional<double>> split_optional_doubles(const std::optional<std::string>& value, size_t n) {
    std::vector<std::optional<double>> out(n, std::nullopt);
    if (!value) return out;
    std::vector<std::string> parts = split_char(*value, ',');
    if (parts.size() != n) throw Error("Argument must have the same number of values as --rg/--h2.");
    for (size_t i = 0; i < parts.size(); ++i) {
        std::replace(parts[i].begin(), parts[i].end(), 'N', '-');
        out[i] = parse_double(parts[i]);
    }
    return out;
}

std::vector<std::string> split_table_line(const std::string& line) {
    if (line.find('\t') == std::string::npos) return split_ws(line);
    std::vector<std::string> out;
    std::string field;
    std::stringstream ss(line);
    while (std::getline(ss, field, '\t')) out.push_back(trim(field));
    if (!line.empty() && line.back() == '\t') out.emplace_back("");
    return out;
}

std::string sec_to_str(double seconds) {
    double whole_part = std::floor(seconds);
    long long whole = static_cast<long long>(whole_part);
    double frac = seconds - whole_part;
    const long long days = whole / 86400;
    whole %= 86400;
    const long long hours = whole / 3600;
    whole %= 3600;
    const long long minutes = whole / 60;
    const long long secs = whole % 60;
    std::ostringstream out;
    if (days > 0) out << days << "d:";
    if (hours > 0) out << hours << "h:";
    if (minutes > 0) out << minutes << "m:";
    out << format_py_float(static_cast<double>(secs) + frac) << "s";
    return out.str();
}

std::string time_now_string() {
    const std::time_t now = std::time(nullptr);
    std::string s = std::ctime(&now);
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
    return s;
}

std::string masthead() {
    std::ostringstream out;
    out << "*********************************************************************\n";
    out << "* LD Score Regression (LDSC)\n";
    out << "* Version 1.0.1\n";
    out << "* (C) 2026 Lulu Shi and GPT\n";
    out << "* West China Hospital of Sichuan University\n";
    out << "* GNU General Public License v3\n";
    out << "*********************************************************************\n";
    return out.str();
}

std::string command_header(int argc, char** argv, const std::string& display_name) {
    std::vector<std::string> args;
    for (int i = 1; i < argc; ++i) args.emplace_back(argv[i]);
    return command_header_from_args(args, display_name);
}

std::string command_header_from_args(const std::vector<std::string>& args, const std::string& display_name) {
    std::ostringstream out;
    out << masthead();
    out << "Call: \n";
    out << display_name << " \\\n";
    for (size_t i = 0; i < args.size();) {
        std::string line = args[i];
        bool boolean_flag = false;
        if (line.rfind("--", 0) == 0 && i + 1 < args.size() && args[i + 1].rfind("--", 0) != 0) {
            line += " " + args[i + 1];
            i += 2;
        } else {
            boolean_flag = line.rfind("--", 0) == 0;
            ++i;
        }
        out << line;
        if (i < args.size()) {
            if (boolean_flag) out << " ";
            out << " \\\n";
        } else {
            if (boolean_flag) out << "  ";
            out << "\n";
        }
    }
    return out.str();
}

bool ends_with(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool is_compressed_name(const std::string& path) {
    return ends_with(path, ".gz") || ends_with(path, ".bz2");
}

std::vector<std::string> read_text_lines(const std::string& path) {
    if (ends_with(path, ".bz2")) {
        std::vector<std::string> lines;
        FILE* fp = std::fopen(path.c_str(), "rb");
        if (!fp) throw Error("Could not open file: " + path);
        int bzerr = BZ_OK;
        BZFILE* bz = BZ2_bzReadOpen(&bzerr, fp, 0, 0, nullptr, 0);
        if (bzerr != BZ_OK) {
            std::fclose(fp);
            throw Error("Could not open bzip2 file: " + path);
        }
        std::string current;
        char buffer[32768];
        while (true) {
            int n = BZ2_bzRead(&bzerr, bz, buffer, sizeof(buffer));
            if (bzerr != BZ_OK && bzerr != BZ_STREAM_END) {
                BZ2_bzReadClose(&bzerr, bz);
                std::fclose(fp);
                throw Error("Error reading bzip2 file: " + path);
            }
            for (int i = 0; i < n; ++i) {
                const char c = buffer[i];
                if (c == '\n') {
                    if (!current.empty() && current.back() == '\r') current.pop_back();
                    lines.push_back(current);
                    current.clear();
                } else {
                    current.push_back(c);
                }
            }
            if (bzerr == BZ_STREAM_END) break;
        }
        if (!current.empty()) {
            if (!current.empty() && current.back() == '\r') current.pop_back();
            lines.push_back(current);
        }
        BZ2_bzReadClose(&bzerr, bz);
        std::fclose(fp);
        return lines;
    }
    std::vector<std::string> lines;
    if (ends_with(path, ".gz")) {
        gzFile gz = gzopen(path.c_str(), "rb");
        if (!gz) throw Error("Could not open file: " + path);
        std::string current;
        char buffer[32768];
        while (gzgets(gz, buffer, sizeof(buffer)) != nullptr) {
            current += buffer;
            if (!current.empty() && current.back() == '\n') {
                while (!current.empty() && (current.back() == '\n' || current.back() == '\r')) current.pop_back();
                lines.push_back(current);
                current.clear();
            }
        }
        if (!current.empty()) {
            while (!current.empty() && (current.back() == '\n' || current.back() == '\r')) current.pop_back();
            lines.push_back(current);
        }
        const int status = gzclose(gz);
        if (status != Z_OK) throw Error("Error closing gzip file: " + path);
    } else {
        std::ifstream in(path);
        if (!in) throw Error("Could not open file: " + path);
        std::string line;
        while (std::getline(in, line)) {
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
            lines.push_back(line);
        }
    }
    return lines;
}

void write_gzip_text(const std::string& path, const std::string& content) {
    gzFile gz = gzopen(path.c_str(), "wb");
    if (!gz) throw Error("Could not open gzip output file: " + path);
    const char* data = content.data();
    size_t remaining = content.size();
    while (remaining > 0) {
        const unsigned chunk = static_cast<unsigned>(std::min<size_t>(remaining, 1u << 30));
        const int written = gzwrite(gz, data, chunk);
        if (written == 0) {
            int err = Z_OK;
            const char* msg = gzerror(gz, &err);
            gzclose(gz);
            throw Error("Error writing gzip output file " + path + ": " + (msg ? std::string(msg) : std::string("unknown error")));
        }
        data += written;
        remaining -= static_cast<size_t>(written);
    }
    const int status = gzclose(gz);
    if (status != Z_OK) throw Error("Error closing gzip output file: " + path);
}

std::vector<std::vector<std::string>> read_rows_no_header(const std::string& path) {
    std::vector<std::vector<std::string>> rows;
    for (std::string line : read_text_lines(path)) {
        line = trim(line);
        if (line.empty()) continue;
        rows.push_back(split_ws(line));
    }
    return rows;
}

double parse_double(const std::string& token) {
    if (token.empty() || token == "." || token == "NA" || token == "nan" || token == "NaN") return kNan;
    char* end = nullptr;
    const double value = std::strtod(token.c_str(), &end);
    if (end == token.c_str() || *end != '\0') {
        throw Error("Could not parse numeric token: " + token);
    }
    return value;
}

std::string format_num(double x, int precision) {
    if (std::isnan(x)) return "NA";
    if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
    if (precision >= 15) {
        std::ostringstream out;
        out << std::setprecision(precision) << x;
        return out.str();
    }
    std::ostringstream out;
    const double ax = std::abs(x);
    if (ax > 0.0 && (ax < 1e-4 || ax >= 1e5)) {
        out << std::scientific << std::setprecision(precision) << x;
        return out.str();
    }
    out << std::fixed << std::setprecision(precision) << x;
    std::string s = out.str();
    while (s.find('.') != std::string::npos && !s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();
    if (s == "-0") s = "0";
    return s;
}

std::string format_py_float(double x) {
    if (std::isnan(x)) return "nan";
    if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
    std::string s = format_num(x, 12);
    if (s.find('.') == std::string::npos && s.find('e') == std::string::npos && s.find('E') == std::string::npos) {
        s += ".0";
    }
    return s;
}

std::string join_nums(const std::vector<double>& xs) {
    std::ostringstream out;
    for (size_t i = 0; i < xs.size(); ++i) {
        if (i) out << ' ';
        out << format_num(xs[i]);
    }
    return out.str();
}

double mean(const std::vector<double>& x) {
    if (x.empty()) throw Error("Cannot take mean of empty vector.");
    return std::accumulate(x.begin(), x.end(), 0.0) / static_cast<double>(x.size());
}

double median(std::vector<double> x) {
    if (x.empty()) throw Error("Cannot take median of empty vector.");
    std::sort(x.begin(), x.end());
    const size_t n = x.size();
    if (n % 2) return x[n / 2];
    return 0.5 * (x[n / 2 - 1] + x[n / 2]);
}

Matrix append_intercept(const Matrix& x) {
    Matrix out(x.rows, x.cols + 1);
    for (size_t i = 0; i < x.rows; ++i) {
        for (size_t j = 0; j < x.cols; ++j) out(i, j) = x(i, j);
        out(i, x.cols) = 1.0;
    }
    return out;
}

std::vector<double> row_sums(const Matrix& x) {
    std::vector<double> out(x.rows, 0.0);
    for (size_t i = 0; i < x.rows; ++i) {
        for (size_t j = 0; j < x.cols; ++j) out[i] += x(i, j);
    }
    return out;
}

std::vector<double> select_vec(const std::vector<double>& x, const std::vector<int>& idx) {
    std::vector<double> out(idx.size());
    for (size_t i = 0; i < idx.size(); ++i) out[i] = x[static_cast<size_t>(idx[i])];
    return out;
}

Matrix weighted_matrix(const Matrix& x, const std::vector<double>& w) {
    if (x.rows != w.size()) throw Error("Weight vector length mismatch.");
    const double sw = std::accumulate(w.begin(), w.end(), 0.0);
    if (!(sw > 0.0)) throw Error("Weights must sum to a positive value.");
    Matrix out(x.rows, x.cols);
    for (size_t i = 0; i < x.rows; ++i) {
        if (!(w[i] > 0.0)) throw Error("Weights must be > 0.");
        const double wi = w[i] / sw;
        for (size_t j = 0; j < x.cols; ++j) out(i, j) = x(i, j) * wi;
    }
    return out;
}

std::vector<double> weighted_vec(const std::vector<double>& y, const std::vector<double>& w) {
    if (y.size() != w.size()) throw Error("Weight vector length mismatch.");
    const double sw = std::accumulate(w.begin(), w.end(), 0.0);
    if (!(sw > 0.0)) throw Error("Weights must sum to a positive value.");
    std::vector<double> out(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        if (!(w[i] > 0.0)) throw Error("Weights must be > 0.");
        out[i] = y[i] * (w[i] / sw);
    }
    return out;
}

std::vector<double> solve_linear(Matrix a, std::vector<double> b) {
    const size_t n = a.rows;
    if (a.rows != a.cols || b.size() != n) throw Error("Linear solve dimension mismatch.");

    for (size_t k = 0; k < n; ++k) {
        size_t pivot = k;
        double best = std::abs(a(k, k));
        for (size_t i = k + 1; i < n; ++i) {
            const double v = std::abs(a(i, k));
            if (v > best) {
                best = v;
                pivot = i;
            }
        }
        if (best < 1e-18) throw Error("Singular matrix in least squares solve.");
        if (pivot != k) {
            for (size_t j = k; j < n; ++j) std::swap(a(k, j), a(pivot, j));
            std::swap(b[k], b[pivot]);
        }
        for (size_t i = k + 1; i < n; ++i) {
            const double factor = a(i, k) / a(k, k);
            a(i, k) = 0.0;
            for (size_t j = k + 1; j < n; ++j) a(i, j) -= factor * a(k, j);
            b[i] -= factor * b[k];
        }
    }

    std::vector<double> x(n, 0.0);
    for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
        double rhs = b[static_cast<size_t>(i)];
        for (size_t j = static_cast<size_t>(i) + 1; j < n; ++j) {
            rhs -= a(static_cast<size_t>(i), j) * x[j];
        }
        x[static_cast<size_t>(i)] = rhs / a(static_cast<size_t>(i), static_cast<size_t>(i));
    }
    return x;
}

std::vector<double> least_squares(const Matrix& x, const std::vector<double>& y) {
    if (x.rows != y.size()) throw Error("Least squares dimension mismatch.");
    Matrix xtx(x.cols, x.cols, 0.0);
    std::vector<double> xty(x.cols, 0.0);
    for (size_t i = 0; i < x.rows; ++i) {
        for (size_t j = 0; j < x.cols; ++j) {
            xty[j] += x(i, j) * y[i];
            for (size_t k = 0; k < x.cols; ++k) xtx(j, k) += x(i, j) * x(i, k);
        }
    }
    return solve_linear(std::move(xtx), std::move(xty));
}

std::vector<int> get_separators(size_t n, int n_blocks) {
    if (n_blocks <= 1) throw Error("--n-blocks must be an integer > 1.");
    if (static_cast<size_t>(n_blocks) > n) throw Error("More jackknife blocks than data points.");
    std::vector<int> out(static_cast<size_t>(n_blocks) + 1);
    for (int i = 0; i <= n_blocks; ++i) {
        out[static_cast<size_t>(i)] = static_cast<int>(std::floor(static_cast<double>(i) * static_cast<double>(n) /
                                                                  static_cast<double>(n_blocks)));
    }
    return out;
}

void jackknife_from_pseudovalues(const Matrix& pseudovalues,
                                 std::vector<double>& jk_est,
                                 std::vector<double>& jk_var,
                                 std::vector<double>& jk_se,
                                 Matrix& jk_cov) {
    const size_t n = pseudovalues.rows;
    const size_t p = pseudovalues.cols;
    jk_est.assign(p, 0.0);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < p; ++j) jk_est[j] += pseudovalues(i, j);
    }
    for (double& v : jk_est) v /= static_cast<double>(n);

    jk_cov = Matrix(p, p, 0.0);
    for (size_t j = 0; j < p; ++j) {
        for (size_t k = 0; k < p; ++k) {
            double acc = 0.0;
            for (size_t i = 0; i < n; ++i) {
                acc += (pseudovalues(i, j) - jk_est[j]) * (pseudovalues(i, k) - jk_est[k]);
            }
            jk_cov(j, k) = acc / static_cast<double>(n - 1) / static_cast<double>(n);
        }
    }
    jk_var.assign(p, 0.0);
    jk_se.assign(p, 0.0);
    for (size_t j = 0; j < p; ++j) {
        jk_var[j] = jk_cov(j, j);
        jk_se[j] = std::sqrt(jk_var[j]);
    }
}

Matrix delete_values_to_pseudovalues(const Matrix& delete_values, const std::vector<double>& est) {
    if (est.size() != delete_values.cols) throw Error("Delete values and estimate dimension mismatch.");
    Matrix out(delete_values.rows, delete_values.cols);
    const double n_blocks = static_cast<double>(delete_values.rows);
    for (size_t i = 0; i < delete_values.rows; ++i) {
        for (size_t j = 0; j < delete_values.cols; ++j) {
            out(i, j) = n_blocks * est[j] - (n_blocks - 1.0) * delete_values(i, j);
        }
    }
    return out;
}

JackknifeResult lstsq_jackknife_fast(const Matrix& x, const std::vector<double>& y, int n_blocks,
                                     std::optional<std::vector<int>> separators) {
    if (x.rows != y.size()) throw Error("Jackknife dimension mismatch.");
    if (x.cols > x.rows) throw Error("More dimensions than datapoints.");
    JackknifeResult out;
    out.separators = separators ? *separators : get_separators(x.rows, n_blocks);
    if (out.separators.front() != 0 || out.separators.back() != static_cast<int>(x.rows)) {
        throw Error("Invalid jackknife separators.");
    }
    const size_t blocks = out.separators.size() - 1;
    const size_t p = x.cols;
    std::vector<Matrix> xtx_block(blocks, Matrix(p, p, 0.0));
    std::vector<std::vector<double>> xty_block(blocks, std::vector<double>(p, 0.0));

    for (size_t b = 0; b < blocks; ++b) {
        const int start = out.separators[b];
        const int end = out.separators[b + 1];
        for (int ii = start; ii < end; ++ii) {
            const size_t i = static_cast<size_t>(ii);
            for (size_t j = 0; j < p; ++j) {
                xty_block[b][j] += x(i, j) * y[i];
                for (size_t k = 0; k < p; ++k) xtx_block[b](j, k) += x(i, j) * x(i, k);
            }
        }
    }

    Matrix xtx_tot(p, p, 0.0);
    std::vector<double> xty_tot(p, 0.0);
    for (size_t b = 0; b < blocks; ++b) {
        for (size_t j = 0; j < p; ++j) {
            xty_tot[j] += xty_block[b][j];
            for (size_t k = 0; k < p; ++k) xtx_tot(j, k) += xtx_block[b](j, k);
        }
    }
    out.est = solve_linear(xtx_tot, xty_tot);
    out.delete_values = Matrix(blocks, p, 0.0);
    for (size_t b = 0; b < blocks; ++b) {
        Matrix delete_xtx = xtx_tot;
        std::vector<double> delete_xty = xty_tot;
        for (size_t j = 0; j < p; ++j) {
            delete_xty[j] -= xty_block[b][j];
            for (size_t k = 0; k < p; ++k) delete_xtx(j, k) -= xtx_block[b](j, k);
        }
        const std::vector<double> del = solve_linear(std::move(delete_xtx), std::move(delete_xty));
        for (size_t j = 0; j < p; ++j) out.delete_values(b, j) = del[j];
    }

    Matrix pseudovalues = delete_values_to_pseudovalues(out.delete_values, out.est);
    jackknife_from_pseudovalues(pseudovalues, out.jknife_est, out.jknife_var, out.jknife_se, out.jknife_cov);
    return out;
}

RatioJackknifeResult ratio_jackknife(const std::vector<double>& est, const Matrix& numer, const Matrix& denom) {
    if (numer.rows != denom.rows || numer.cols != denom.cols || est.size() != numer.cols) {
        throw Error("Ratio jackknife dimension mismatch.");
    }
    Matrix pseudovalues(numer.rows, numer.cols);
    const double n_blocks = static_cast<double>(numer.rows);
    for (size_t i = 0; i < numer.rows; ++i) {
        for (size_t j = 0; j < numer.cols; ++j) {
            pseudovalues(i, j) = n_blocks * est[j] - (n_blocks - 1.0) * numer(i, j) / denom(i, j);
        }
    }
    RatioJackknifeResult out;
    out.est = est;
    jackknife_from_pseudovalues(pseudovalues, out.jknife_est, out.jknife_var, out.jknife_se, out.jknife_cov);
    return out;
}

JackknifeResult irwls(const Matrix& x, const std::vector<double>& y, const UpdateFunc& update_func, int n_blocks,
                      const std::vector<double>& initial_inverse_cvf,
                      std::optional<std::vector<int>> separators,
                      const Matrix* ref_ld_for_update,
                      const std::vector<double>* w_ld_for_update,
                      const std::vector<double>* n_core_for_update,
                      double m_tot, double nbar,
                      std::optional<double> intercept,
                      const std::vector<int>* source_rows) {
    if (x.rows != y.size() || y.size() != initial_inverse_cvf.size()) throw Error("IRWLS dimension mismatch.");
    std::vector<double> w(initial_inverse_cvf.size());
    for (size_t i = 0; i < w.size(); ++i) w[i] = std::sqrt(initial_inverse_cvf[i]);

    const Matrix& ref = ref_ld_for_update ? *ref_ld_for_update : x;
    const std::vector<double>& wld = w_ld_for_update ? *w_ld_for_update : initial_inverse_cvf;
    const std::vector<double>& ncore = n_core_for_update ? *n_core_for_update : initial_inverse_cvf;
    for (int iter = 0; iter < 2; ++iter) {
        Matrix xw = weighted_matrix(x, w);
        std::vector<double> yw = weighted_vec(y, w);
        std::vector<double> coef = least_squares(xw, yw);
        std::vector<double> new_inverse = update_func(coef, ref, wld, ncore, m_tot, nbar, intercept, source_rows);
        if (new_inverse.size() != w.size()) throw Error("IRWLS updated weights have wrong length.");
        for (size_t i = 0; i < w.size(); ++i) w[i] = std::sqrt(new_inverse[i]);
    }
    Matrix xw = weighted_matrix(x, w);
    std::vector<double> yw = weighted_vec(y, w);
    return lstsq_jackknife_fast(xw, yw, n_blocks, separators);
}

std::vector<int> update_separators(const std::vector<int>& masked_separators, const std::vector<bool>& mask) {
    std::vector<int> kept;
    for (size_t i = 0; i < mask.size(); ++i) {
        if (mask[i]) kept.push_back(static_cast<int>(i));
    }
    std::vector<int> out;
    out.push_back(0);
    for (size_t i = 1; i + 1 < masked_separators.size(); ++i) {
        const int idx = masked_separators[i];
        if (idx < 0 || static_cast<size_t>(idx) >= kept.size()) throw Error("Invalid masked separator.");
        out.push_back(kept[static_cast<size_t>(idx)]);
    }
    out.push_back(static_cast<int>(mask.size()));
    return out;
}

Table read_table(const std::string& path) {
    Table table;
    std::vector<std::string> lines = read_text_lines(path);
    size_t line_idx = 0;
    for (; line_idx < lines.size(); ++line_idx) {
        std::string line = lines[line_idx];
        line = trim(line);
        if (!line.empty()) {
            table.columns = split_table_line(line);
            ++line_idx;
            break;
        }
    }
    if (table.columns.empty()) throw Error("Empty table: " + path);
    for (size_t i = 0; i < table.columns.size(); ++i) table.index[table.columns[i]] = i;
    for (; line_idx < lines.size(); ++line_idx) {
        std::string line = lines[line_idx];
        line = trim(line);
        if (line.empty()) continue;
        std::vector<std::string> row = split_table_line(line);
        if (row.size() < table.columns.size()) {
            throw Error("Malformed row in " + path + ": " + line);
        }
        table.rows.push_back(std::move(row));
    }
    return table;
}

std::string resolve_existing(const std::vector<std::string>& candidates, const std::string& label) {
    for (const auto& c : candidates) {
        if (fs::exists(c)) return c;
    }
    std::ostringstream msg;
    msg << "Could not find " << label << ". Tried:";
    for (const auto& c : candidates) msg << "\n  " << c;
    throw Error(msg.str());
}

std::string resolve_ldscore_file(const std::string& prefix) {
    std::vector<std::string> candidates;
    if (ends_with(prefix, ".l2.ldscore") || ends_with(prefix, ".l2.ldscore.gz") || ends_with(prefix, ".l2.ldscore.bz2")) {
        candidates.push_back(prefix);
    }
    candidates.push_back(prefix + ".l2.ldscore");
    candidates.push_back(prefix + ".l2.ldscore.gz");
    candidates.push_back(prefix + ".l2.ldscore.bz2");
    return resolve_existing(candidates, "LD Score file for prefix " + prefix);
}

std::string resolve_m_file(const std::string& prefix, bool common) {
    const std::string suffix = common ? ".l2.M_5_50" : ".l2.M";
    return resolve_existing({prefix + suffix}, "M file for prefix " + prefix);
}

std::string sub_chr(std::string pattern, int chr) {
    const std::string replacement = std::to_string(chr);
    const size_t pos = pattern.find('@');
    if (pos == std::string::npos) return pattern + replacement;
    size_t start = 0;
    while ((start = pattern.find('@', start)) != std::string::npos) {
        pattern.replace(start, 1, replacement);
        start += replacement.size();
    }
    return pattern;
}

std::string sub_chr_display(std::string pattern) {
    const std::string replacement = "[1-22]";
    if (pattern.find('@') == std::string::npos) pattern += '@';
    size_t start = 0;
    while ((start = pattern.find('@', start)) != std::string::npos) {
        pattern.replace(start, 1, replacement);
        start += replacement.size();
    }
    return pattern;
}

std::vector<std::string> present_chr_ldscore_files(const std::string& prefix) {
    std::vector<std::string> out;
    for (int chr = 1; chr <= 22; ++chr) {
        const std::string chr_prefix = sub_chr(prefix, chr);
        std::vector<std::string> candidates{chr_prefix + ".l2.ldscore", chr_prefix + ".l2.ldscore.gz",
                                            chr_prefix + ".l2.ldscore.bz2"};
        for (const auto& c : candidates) {
            if (fs::exists(c)) {
                out.push_back(c);
                break;
            }
        }
    }
    if (out.empty()) throw Error("No chromosome-split LD Score files found for prefix: " + prefix);
    return out;
}

std::vector<std::string> present_chr_m_files(const std::string& prefix, bool common) {
    const std::string suffix = common ? ".l2.M_5_50" : ".l2.M";
    std::vector<std::string> out;
    for (int chr = 1; chr <= 22; ++chr) {
        const std::string path = sub_chr(prefix, chr) + suffix;
        if (fs::exists(path)) out.push_back(path);
    }
    if (out.empty()) throw Error("No chromosome-split M files found for prefix: " + prefix);
    return out;
}

LdFrame read_ldscore_paths(const std::vector<std::string>& paths) {
    std::vector<LdTempRow> temp;
    std::vector<std::string> value_names;
    bool have_sort_cols = false;
    size_t order = 0;
    for (const auto& path : paths) {
        Table t = read_table(path);
        const int snp_i = t.col("SNP");
        if (snp_i < 0) throw Error("LD Score file lacks SNP column: " + path);
        const int chr_i = t.col("CHR");
        const int bp_i = t.col("BP");
        have_sort_cols = have_sort_cols || (chr_i >= 0 && bp_i >= 0);
        std::vector<size_t> value_idx;
        if (value_names.empty()) {
            for (size_t i = 0; i < t.columns.size(); ++i) {
                const std::string& c = t.columns[i];
                if (c == "SNP" || c == "CHR" || c == "BP" || c == "CM" || c == "MAF") continue;
                value_idx.push_back(i);
                value_names.push_back(c);
            }
        } else {
            for (size_t i = 0; i < t.columns.size(); ++i) {
                const std::string& c = t.columns[i];
                if (c == "SNP" || c == "CHR" || c == "BP" || c == "CM" || c == "MAF") continue;
                value_idx.push_back(i);
            }
            if (value_idx.size() != value_names.size()) {
                throw Error("Chromosome-split LD Score files have different numbers of LD columns.");
            }
        }
        for (const auto& row : t.rows) {
            LdTempRow r;
            r.order = order++;
            r.snp = row[static_cast<size_t>(snp_i)];
            if (chr_i >= 0) r.chr = parse_double(row[static_cast<size_t>(chr_i)]);
            if (bp_i >= 0) r.bp = parse_double(row[static_cast<size_t>(bp_i)]);
            for (size_t idx : value_idx) r.values.push_back(parse_double(row[idx]));
            temp.push_back(std::move(r));
        }
    }
    if (have_sort_cols) {
        std::sort(temp.begin(), temp.end(), [](const LdTempRow& a, const LdTempRow& b) {
            if (a.chr != b.chr) return a.chr < b.chr;
            if (a.bp != b.bp) return a.bp < b.bp;
            return a.order < b.order;
        });
    }
    std::unordered_set<std::string> seen;
    std::vector<LdTempRow> unique;
    for (auto& r : temp) {
        if (seen.insert(r.snp).second) unique.push_back(std::move(r));
    }

    LdFrame out;
    out.names = std::move(value_names);
    out.snp.reserve(unique.size());
    out.values = Matrix(unique.size(), out.names.size());
    for (size_t i = 0; i < unique.size(); ++i) {
        out.snp.push_back(unique[i].snp);
        for (size_t j = 0; j < out.names.size(); ++j) out.values(i, j) = unique[i].values[j];
    }
    return out;
}

LdFrame read_ldscore_prefix(const std::string& prefix, bool chr_split) {
    if (chr_split) return read_ldscore_paths(present_chr_ldscore_files(prefix));
    return read_ldscore_paths({resolve_ldscore_file(prefix)});
}

LdFrame read_ldscore_list(const std::string& arg, bool chr_split) {
    const std::vector<std::string> prefixes = split_char(arg, ',');
    LdFrame combined;
    for (size_t i = 0; i < prefixes.size(); ++i) {
        LdFrame current = read_ldscore_prefix(prefixes[i], chr_split);
        if (i == 0) {
            combined.snp = current.snp;
            combined.values = Matrix(current.values.rows, 0);
        } else if (current.snp != combined.snp) {
            throw Error("LD Scores for concatenation must have identical SNP columns.");
        }
        const size_t old_cols = combined.names.size();
        Matrix next(combined.snp.size(), old_cols + current.names.size());
        for (size_t r = 0; r < combined.snp.size(); ++r) {
            for (size_t c = 0; c < old_cols; ++c) next(r, c) = combined.values(r, c);
            for (size_t c = 0; c < current.names.size(); ++c) next(r, old_cols + c) = current.values(r, c);
        }
        combined.values = std::move(next);
        for (const auto& name : current.names) combined.names.push_back(name + "_" + std::to_string(i));
    }
    return combined;
}

std::vector<double> read_m_line(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw Error("Could not open M file: " + path);
    std::string line;
    std::getline(in, line);
    std::vector<double> out;
    for (const auto& tok : split_ws(line)) out.push_back(parse_double(tok));
    if (out.empty()) throw Error("Empty M file: " + path);
    return out;
}

std::vector<double> read_m_prefix(const std::string& prefix, bool chr_split, bool common) {
    if (!chr_split) return read_m_line(resolve_m_file(prefix, common));
    std::vector<double> acc;
    for (const auto& path : present_chr_m_files(prefix, common)) {
        std::vector<double> row = read_m_line(path);
        if (acc.empty()) acc.assign(row.size(), 0.0);
        if (row.size() != acc.size()) throw Error("Chromosome-split M files have inconsistent widths.");
        for (size_t i = 0; i < row.size(); ++i) acc[i] += row[i];
    }
    return acc;
}

std::vector<double> read_m(const std::optional<std::string>& manual_m, const std::string& ref_arg, bool chr_split,
                           bool common) {
    if (manual_m) {
        std::vector<double> out;
        for (const auto& tok : split_char(*manual_m, ',')) out.push_back(parse_double(tok));
        return out;
    }
    std::vector<double> out;
    for (const auto& prefix : split_char(ref_arg, ',')) {
        std::vector<double> part = read_m_prefix(prefix, chr_split, common);
        out.insert(out.end(), part.begin(), part.end());
    }
    return out;
}

void remove_zero_variance(LdFrame& ref_ld, std::vector<double>& m_annot) {
    if (ref_ld.values.cols != m_annot.size()) throw Error("# terms in M must match # of LD Scores in --ref-ld.");
    std::vector<bool> keep(ref_ld.values.cols, true);
    size_t kept = 0;
    for (size_t c = 0; c < ref_ld.values.cols; ++c) {
        bool zero = false;
        if (ref_ld.values.rows >= 2) {
            double mu = 0.0;
            for (size_t r = 0; r < ref_ld.values.rows; ++r) mu += ref_ld.values(r, c);
            mu /= static_cast<double>(ref_ld.values.rows);
            double ss = 0.0;
            for (size_t r = 0; r < ref_ld.values.rows; ++r) ss += (ref_ld.values(r, c) - mu) * (ref_ld.values(r, c) - mu);
            zero = ss == 0.0;
        }
        keep[c] = !zero;
        if (keep[c]) ++kept;
    }
    if (kept == 0) throw Error("All LD Scores have zero variance.");
    if (kept == ref_ld.values.cols) return;

    Matrix values(ref_ld.values.rows, kept);
    std::vector<std::string> names;
    std::vector<double> m;
    size_t dst = 0;
    for (size_t c = 0; c < keep.size(); ++c) {
        if (!keep[c]) continue;
        for (size_t r = 0; r < ref_ld.values.rows; ++r) values(r, dst) = ref_ld.values(r, c);
        names.push_back(ref_ld.names[c]);
        m.push_back(m_annot[c]);
        ++dst;
    }
    ref_ld.values = std::move(values);
    ref_ld.names = std::move(names);
    m_annot = std::move(m);
}

std::vector<SumRecord> read_sumstats(const std::string& path, bool alleles, bool dropna) {
    Table t = read_table(path);
    const int snp_i = t.col("SNP");
    const int z_i = t.col("Z");
    const int n_i = t.col("N");
    const int a1_i = t.col("A1");
    const int a2_i = t.col("A2");
    if (snp_i < 0 || z_i < 0 || n_i < 0 || (alleles && (a1_i < 0 || a2_i < 0))) {
        throw Error("Improperly formatted sumstats file: " + path);
    }

    std::vector<SumRecord> out;
    std::unordered_set<std::string> seen;
    for (const auto& row : t.rows) {
        SumRecord r;
        r.snp = row[static_cast<size_t>(snp_i)];
        r.z = parse_double(row[static_cast<size_t>(z_i)]);
        r.n = parse_double(row[static_cast<size_t>(n_i)]);
        if (alleles) {
            r.a1 = row[static_cast<size_t>(a1_i)];
            r.a2 = row[static_cast<size_t>(a2_i)];
        }
        if (dropna && (std::isnan(r.z) || std::isnan(r.n) || (alleles && (r.a1 == "." || r.a2 == ".")))) continue;
        if (seen.insert(r.snp).second) out.push_back(std::move(r));
    }
    return out;
}

std::unordered_map<std::string, size_t> index_by_snp(const std::vector<std::string>& snps) {
    std::unordered_map<std::string, size_t> out;
    for (size_t i = 0; i < snps.size(); ++i) out.emplace(snps[i], i);
    return out;
}

BaseData read_ld_sumstats(const std::string& sumstats_path, const std::string& ref_arg, bool ref_chr,
                          const std::string& w_arg, bool w_chr, const std::optional<std::string>& manual_m,
                          bool common_m, std::vector<double>& m_annot, Logger& log,
                          bool alleles, bool dropna) {
    log.log("Reading summary statistics from " + sumstats_path + " ...");
    std::vector<SumRecord> sumstats = read_sumstats(sumstats_path, alleles, dropna);
    log.log("Read summary statistics for " + std::to_string(sumstats.size()) + " SNPs.");

    log.log("Reading reference panel LD Score from " + (ref_chr ? sub_chr_display(ref_arg) : ref_arg) +
            " ... (ldscore_fromlist)");
    LdFrame ref_ld = read_ldscore_list(ref_arg, ref_chr);
    log.log("Read reference panel LD Scores for " + std::to_string(ref_ld.snp.size()) + " SNPs.");
    m_annot = read_m(manual_m, ref_arg, ref_chr, common_m);
    log.log("Removing partitioned LD Scores with zero variance.");
    remove_zero_variance(ref_ld, m_annot);

    log.log("Reading regression weight LD Score from " + (w_chr ? sub_chr_display(w_arg) : w_arg) +
            " ... (ldscore_fromlist)");
    LdFrame w_ld = read_ldscore_list(w_arg, w_chr);
    if (w_ld.values.cols != 1) throw Error("--w-ld may only have one LD Score column.");
    log.log("Read regression weight LD Scores for " + std::to_string(w_ld.snp.size()) + " SNPs.");

    std::unordered_map<std::string, size_t> sum_i;
    for (size_t i = 0; i < sumstats.size(); ++i) sum_i.emplace(sumstats[i].snp, i);
    std::unordered_map<std::string, size_t> w_i = index_by_snp(w_ld.snp);

    BaseData out;
    out.ld_names = ref_ld.names;
    std::vector<std::vector<double>> ld_rows;
    size_t after_ref = 0;
    for (size_t i = 0; i < ref_ld.snp.size(); ++i) {
        auto si = sum_i.find(ref_ld.snp[i]);
        if (si == sum_i.end()) continue;
        ++after_ref;
        auto wi = w_i.find(ref_ld.snp[i]);
        if (wi == w_i.end()) continue;
        const SumRecord& sr = sumstats[si->second];
        out.snp.push_back(ref_ld.snp[i]);
        out.z.push_back(sr.z);
        out.n.push_back(sr.n);
        out.w_ld.push_back(w_ld.values(wi->second, 0));
        if (alleles) {
            out.a1.push_back(sr.a1);
            out.a2.push_back(sr.a2);
        }
        std::vector<double> row(ref_ld.values.cols);
        for (size_t c = 0; c < ref_ld.values.cols; ++c) row[c] = ref_ld.values(i, c);
        ld_rows.push_back(std::move(row));
    }
    log.log("After merging with reference panel LD, " + std::to_string(after_ref) + " SNPs remain.");
    log.log("After merging with regression SNP LD, " + std::to_string(out.snp.size()) + " SNPs remain.");
    if (out.snp.empty()) throw Error("No SNPs retained for analysis.");
    out.ld = Matrix(out.snp.size(), ref_ld.values.cols);
    for (size_t r = 0; r < ld_rows.size(); ++r) {
        for (size_t c = 0; c < ld_rows[r].size(); ++c) out.ld(r, c) = ld_rows[r][c];
    }
    return out;
}

std::vector<double> hsq_weights(const std::vector<double>& ld, const std::vector<double>& w_ld,
                                const std::vector<double>& n, double m, double hsq,
                                std::optional<double> intercept) {
    if (ld.size() != w_ld.size() || ld.size() != n.size()) throw Error("Hsq weight dimension mismatch.");
    const double intc = intercept.value_or(1.0);
    hsq = std::max(0.0, std::min(1.0, hsq));
    std::vector<double> out(ld.size());
    for (size_t i = 0; i < ld.size(); ++i) {
        const double l = std::max(ld[i], 1.0);
        const double wl = std::max(w_ld[i], 1.0);
        const double c = hsq * n[i] / m;
        const double het_w = 1.0 / (2.0 * std::pow(intc + c * l, 2.0));
        out[i] = het_w / wl;
    }
    return out;
}

std::vector<double> gencov_weights(const std::vector<double>& ld, const std::vector<double>& w_ld,
                                   const std::vector<double>& n1, const std::vector<double>& n2,
                                   double m, double h1, double h2, double rho_g,
                                   std::optional<double> intercept_gencov,
                                   std::optional<double> intercept_hsq1,
                                   std::optional<double> intercept_hsq2) {
    if (ld.size() != w_ld.size() || ld.size() != n1.size() || ld.size() != n2.size()) {
        throw Error("Gencov weight dimension mismatch.");
    }
    h1 = std::max(0.0, std::min(1.0, h1));
    h2 = std::max(0.0, std::min(1.0, h2));
    rho_g = std::max(-1.0, std::min(1.0, rho_g));
    const double int_g = intercept_gencov.value_or(0.0);
    const double int_h1 = intercept_hsq1.value_or(1.0);
    const double int_h2 = intercept_hsq2.value_or(1.0);
    std::vector<double> out(ld.size());
    for (size_t i = 0; i < ld.size(); ++i) {
        const double l = std::max(ld[i], 1.0);
        const double wl = std::max(w_ld[i], 1.0);
        const double a = n1[i] * h1 * l / m + int_h1;
        const double b = n2[i] * h2 * l / m + int_h2;
        const double c = std::sqrt(n1[i] * n2[i]) * rho_g * l / m + int_g;
        out[i] = 1.0 / (a * b + c * c) / wl;
    }
    return out;
}

void save_matrix(const std::string& path, const Matrix& m) {
    std::ofstream out(path);
    if (!out) throw Error("Could not open output file: " + path);
    out << std::scientific << std::setprecision(18);
    for (size_t r = 0; r < m.rows; ++r) {
        for (size_t c = 0; c < m.cols; ++c) {
            if (c) out << ' ';
            out << m(r, c);
        }
        out << '\n';
    }
}

void ensure_parent_dir(const std::string& out_prefix) {
    fs::path p(out_prefix);
    if (p.has_parent_path()) fs::create_directories(p.parent_path());
}

char complement_base(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return '?';
    }
}

bool valid_base(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

bool alleles_match(const std::string& a) {
    if (a.size() != 4) return false;
    for (char c : a) {
        if (!valid_base(c)) return false;
    }
    if (a[0] == a[1] || a[2] == a[3]) return false;
    const std::string s1 = a.substr(0, 2);
    const std::string s2 = a.substr(2, 2);
    auto strand_ambig = [](const std::string& s) { return s[0] == complement_base(s[1]); };
    if (strand_ambig(s1) || strand_ambig(s2)) return false;
    const bool strand_ref = a[0] == a[2] && a[1] == a[3];
    const bool strand_flip = a[0] == complement_base(a[2]) && a[1] == complement_base(a[3]);
    const bool ref_flip = a[0] == a[3] && a[1] == a[2];
    const bool both_flip = a[0] == complement_base(a[3]) && a[1] == complement_base(a[2]);
    return strand_ref || strand_flip || ref_flip || both_flip;
}

bool flip_alleles(const std::string& a) {
    const bool ref_flip = a[0] == a[3] && a[1] == a[2];
    const bool both_flip = a[0] == complement_base(a[3]) && a[1] == complement_base(a[2]);
    return ref_flip || both_flip;
}

}  // namespace ldsc_cpp
