#include "continuous.hpp"

namespace ldsc_cpp {

std::vector<double> read_numeric_row_tail(const std::vector<std::string>& row) {
    std::vector<double> out;
    for (size_t i = 1; i < row.size(); ++i) out.push_back(parse_double(row[i]));
    return out;
}

double jackknife_se(double theta, const std::vector<double>& theta_j) {
    const double n = static_cast<double>(theta_j.size());
    double theta_hat = n * theta;
    for (double v : theta_j) theta_hat -= (n - 1.0) * v / n;
    double acc = 0.0;
    for (double v : theta_j) {
        const double tau = n * theta - (n - 1.0) * v;
        acc += (tau - theta_hat) * (tau - theta_hat) / (n - 1.0);
    }
    return std::sqrt(acc / n);
}

void quantile_h2g_cpp(const Args& args) {
    std::string annotfile;
    std::string resultfile;
    std::string outfile;
    if (args.positionals.size() >= 3) {
        annotfile = args.positionals[0];
        resultfile = args.positionals[1];
        outfile = args.positionals[2];
    } else {
        if (!args.quantile_annotfile || !args.quantile_resultfile || !args.quantile_outfile) {
            throw Error("quantile-h2g requires annotfile, resultfile prefix, and outfile.");
        }
        annotfile = *args.quantile_annotfile;
        resultfile = *args.quantile_resultfile;
        outfile = *args.quantile_outfile;
    }
    std::vector<std::vector<std::string>> annot_rows = read_rows_no_header(annotfile);
    if (annot_rows.size() < 2) throw Error("Quantile annotation file is empty.");
    std::vector<double> chunk_size = read_numeric_row_tail(annot_rows[0]);
    const size_t q = chunk_size.size();
    Matrix baseline(annot_rows.size() - 1, q);
    for (size_t i = 1; i < annot_rows.size(); ++i) {
        std::vector<double> vals = read_numeric_row_tail(annot_rows[i]);
        if (vals.size() != q) throw Error("Quantile annotation row has inconsistent width.");
        for (size_t j = 0; j < q; ++j) baseline(i - 1, j) = vals[j];
    }
    Table results = read_table(resultfile + ".results");
    int coef_i = results.col("Coefficient");
    if (coef_i < 0) throw Error(".results file lacks Coefficient column.");
    std::vector<double> tau;
    for (const auto& row : results.rows) tau.push_back(parse_double(row[static_cast<size_t>(coef_i)]));
    if (tau.size() != baseline.rows) throw Error("Coefficient count does not match quantile annotation rows.");
    auto compute_h2 = [&](const std::vector<double>& coef) {
        std::vector<double> h(q, 0.0);
        for (size_t a = 0; a < baseline.rows; ++a) {
            for (size_t j = 0; j < q; ++j) h[j] += coef[a] * baseline(a, j);
        }
        return h;
    };
    std::vector<double> h2g = compute_h2(tau);
    const double h2_sum = std::accumulate(h2g.begin(), h2g.end(), 0.0);
    const double m_sum = std::accumulate(chunk_size.begin(), chunk_size.end(), 0.0);
    std::vector<double> h2g_prop(q), enr(q), enr_z(q);
    for (size_t j = 0; j < q; ++j) {
        h2g_prop[j] = h2g[j] / h2_sum;
        enr[j] = h2g_prop[j] / (chunk_size[j] / m_sum);
        enr_z[j] = (h2g[j] / chunk_size[j]) - ((h2_sum - h2g[j]) / (m_sum - chunk_size[j]));
    }
    std::vector<std::vector<std::string>> jack_rows = read_rows_no_header(resultfile + ".part_delete");
    std::vector<std::vector<double>> h2_j(q), prop_j(q), enr_z_j(q);
    for (const auto& row : jack_rows) {
        std::vector<double> coef;
        for (const auto& x : row) coef.push_back(parse_double(x));
        std::vector<double> h = compute_h2(coef);
        double hs = std::accumulate(h.begin(), h.end(), 0.0);
        for (size_t j = 0; j < q; ++j) {
            h2_j[j].push_back(h[j]);
            prop_j[j].push_back(h[j] / hs);
            enr_z_j[j].push_back((h[j] / chunk_size[j]) - ((hs - h[j]) / (m_sum - chunk_size[j])));
        }
    }
    std::ofstream out(outfile);
    if (!out) throw Error("Could not open quantile-h2g output: " + outfile);
    out << "h2g\th2g_se\tprop_h2g\tprop_h2g_se\tenr\tenr_se\tenr_pval\n";
    out << std::setprecision(10);
    for (size_t j = 0; j < q; ++j) {
        double h2_se = jackknife_se(h2g[j], h2_j[j]);
        double prop_se = jackknife_se(h2g_prop[j], prop_j[j]);
        double enr_z_se = jackknife_se(enr_z[j], enr_z_j[j]);
        double enr_se = prop_se / (chunk_size[j] / m_sum);
        double p = std::erfc(std::abs(enr_z[j] / enr_z_se) / std::sqrt(2.0));
        out << h2g[j] << '\t' << h2_se << '\t' << h2g_prop[j] << '\t' << prop_se << '\t'
            << enr[j] << '\t' << enr_se << '\t' << p << "\n";
    }
}

std::string chr_file_existing(const std::string& prefix, int chr, const std::string& suffix) {
    std::string path = prefix + std::to_string(chr) + suffix;
    if (fs::exists(path)) return path;
    throw Error("Could not open chromosome-split file: " + path);
}

void quantile_m_cpp(const Args& args) {
    if (!args.frqfile_chr || !args.ref_annot_chr || !args.annot_header || !args.out_set) {
        throw Error("quantile-m requires --frqfile-chr, --ref-annot-chr, --annot-header, and --out.");
    }
    const double maf_threshold = args.maf_set ? args.maf : 0.05;
    if (maf_threshold < 0.0 || maf_threshold > 0.5) throw Error("--maf should be between 0 and 0.5.");
    std::vector<bool> common;
    for (int chr = 1; chr <= 22; ++chr) {
        Table frq = read_table(chr_file_existing(*args.frqfile_chr, chr, ".frq"));
        int maf_i = frq.col("MAF");
        if (maf_i < 0) maf_i = frq.col("FRQ");
        if (maf_i < 0 && frq.columns.size() > 4) maf_i = 4;
        if (maf_i < 0) throw Error("Frequency file lacks MAF/FRQ column.");
        for (const auto& row : frq.rows) {
            double maf = parse_double(row[static_cast<size_t>(maf_i)]);
            common.push_back(maf >= maf_threshold && maf <= 1.0 - maf_threshold);
        }
    }
    std::vector<std::string> prefixes = split_char(*args.ref_annot_chr, ',');
    std::string annot_prefix;
    int annot_col = -1;
    size_t total_annots = 0;
    for (const auto& pref : prefixes) {
        Table a = read_table(chr_file_existing(pref, 22, ".annot"));
        for (size_t i = 0; i < a.columns.size(); ++i) {
            if (a.columns[i] == *args.annot_header) {
                annot_prefix = pref;
                annot_col = static_cast<int>(i);
            }
        }
        if (a.columns.size() > 4) total_annots += a.columns.size() - 4;
    }
    if (annot_col < 0) throw Error("Annotation " + *args.annot_header + " not found in --ref-annot-chr files.");
    std::vector<double> annot_all;
    std::vector<double> annot_for_quantile;
    size_t snp_i = 0;
    for (int chr = 1; chr <= 22; ++chr) {
        Table a = read_table(chr_file_existing(annot_prefix, chr, ".annot"));
        for (const auto& row : a.rows) {
            double v = parse_double(row[static_cast<size_t>(annot_col)]);
            annot_all.push_back(v);
            if (snp_i < common.size() && common[snp_i] && (!args.exclude0 || v != 0.0)) annot_for_quantile.push_back(v);
            ++snp_i;
        }
    }
    if (annot_for_quantile.empty()) throw Error("No annotation values available for quantiles.");
    std::sort(annot_for_quantile.begin(), annot_for_quantile.end());
    std::vector<double> qv;
    for (int i = 0; i <= args.nb_quantile; ++i) {
        size_t idx = static_cast<size_t>(std::llround(static_cast<double>(i) *
                                                      static_cast<double>(annot_for_quantile.size() - 1) /
                                                      static_cast<double>(args.nb_quantile)));
        qv.push_back(annot_for_quantile[idx]);
    }
    qv.front() -= 1e-10;
    Matrix matrix(total_annots, static_cast<size_t>(args.nb_quantile), 0.0);
    std::vector<double> nb_per_q(static_cast<size_t>(args.nb_quantile), 0.0);
    std::vector<std::string> annot_names;
    size_t annot_offset = 0;
    for (const auto& pref : prefixes) {
        bool names_added = false;
        snp_i = 0;
        for (int chr = 1; chr <= 22; ++chr) {
            Table a = read_table(chr_file_existing(pref, chr, ".annot"));
            if (!names_added) {
                for (size_t c = 4; c < a.columns.size(); ++c) annot_names.push_back(a.columns[c]);
                names_added = true;
            }
            for (const auto& row : a.rows) {
                if (snp_i < common.size() && common[snp_i]) {
                    double av = annot_all[snp_i];
                    if (!args.exclude0 || av != 0.0) {
                        int bin = -1;
                        while (bin + 1 < args.nb_quantile && av > qv[static_cast<size_t>(bin + 1)]) ++bin;
                        if (bin < 0) bin = 0;
                        if (bin >= args.nb_quantile) bin = args.nb_quantile - 1;
                        for (size_t c = 4; c < a.columns.size(); ++c) {
                            matrix(annot_offset + c - 4, static_cast<size_t>(bin)) += parse_double(row[c]);
                        }
                        nb_per_q[static_cast<size_t>(bin)] += 1.0;
                    }
                }
                ++snp_i;
            }
        }
        Table first = read_table(chr_file_existing(pref, 22, ".annot"));
        annot_offset += first.columns.size() > 4 ? first.columns.size() - 4 : 0;
    }
    std::ofstream out(args.out);
    if (!out) throw Error("Could not open quantile-m output: " + args.out);
    out << "N";
    for (double n : nb_per_q) out << '\t' << n;
    out << "\n";
    for (size_t i = 0; i < matrix.rows; ++i) {
        out << (i < annot_names.size() ? annot_names[i] : ("ANNOT" + std::to_string(i)));
        for (size_t j = 0; j < matrix.cols; ++j) out << '\t' << matrix(i, j);
        out << "\n";
    }
}

}  // namespace ldsc_cpp
