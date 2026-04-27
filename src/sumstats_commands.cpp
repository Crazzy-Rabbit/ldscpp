#include "sumstats_commands.hpp"

namespace ldsc_cpp {

struct OverlapAnnotData {
    Matrix overlap;
    double m_tot = 0.0;
};

std::string py2_section(const std::string& title) {
    return "\n" + title + "\n" + std::string(title.size(), '-');
}

std::string format_rg_p_summary(double p) {
    if (p == 0.0) return "0.";
    return format_num(p);
}

std::string format_rg_p_table(double p) {
    if (p == 0.0) return "0.0";
    return format_num(p);
}

std::string rg_table_row(const std::string& p1, const std::string& p2, const RGResult& r) {
    const std::vector<std::string> values = {
        format_num(r.rg_ratio),
        format_num(r.rg_se),
        format_num(r.z_score),
        format_rg_p_table(r.p_value),
        format_num(r.hsq2.tot),
        format_num(r.hsq2.tot_se),
        format_num(r.hsq2.intercept),
        format_num(r.hsq2.intercept_se),
        format_num(r.gencov.intercept),
        format_num(r.gencov.intercept_se),
    };
    const std::vector<int> widths = {6, 6, 8, 3, 6, 9, 6, 9, 8, 11};
    std::ostringstream out;
    out << p1 << "  " << p2;
    for (size_t i = 0; i < values.size(); ++i) out << "  " << std::setw(widths[i]) << values[i];
    return out.str();
}

double beta_continued_fraction(double a, double b, double x) {
    constexpr int max_iter = 200;
    constexpr double eps = 3e-14;
    constexpr double fpmin = 1e-300;
    const double qab = a + b;
    const double qap = a + 1.0;
    const double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::abs(d) < fpmin) d = fpmin;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= max_iter; ++m) {
        const int m2 = 2 * m;
        double aa = static_cast<double>(m) * (b - static_cast<double>(m)) * x /
                    ((qam + static_cast<double>(m2)) * (a + static_cast<double>(m2)));
        d = 1.0 + aa * d;
        if (std::abs(d) < fpmin) d = fpmin;
        c = 1.0 + aa / c;
        if (std::abs(c) < fpmin) c = fpmin;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + static_cast<double>(m)) * (qab + static_cast<double>(m)) * x /
             ((a + static_cast<double>(m2)) * (qap + static_cast<double>(m2)));
        d = 1.0 + aa * d;
        if (std::abs(d) < fpmin) d = fpmin;
        c = 1.0 + aa / c;
        if (std::abs(c) < fpmin) c = fpmin;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (std::abs(del - 1.0) <= eps) break;
    }
    return h;
}

double regularized_beta(double x, double a, double b) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    const double bt = std::exp(std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
                               a * std::log(x) + b * std::log1p(-x));
    if (x < (a + 1.0) / (a + b + 2.0)) return bt * beta_continued_fraction(a, b, x) / a;
    return 1.0 - bt * beta_continued_fraction(b, a, 1.0 - x) / b;
}

double student_t_two_sided_p(double t, int df) {
    if (!(df > 0)) return kNan;
    if (std::isnan(t)) return kNan;
    if (std::isinf(t)) return 0.0;
    const double v = static_cast<double>(df);
    const double x = v / (v + t * t);
    return regularized_beta(x, 0.5 * v, 0.5);
}

std::string resolve_annot_file(const std::string& prefix) {
    return resolve_existing({prefix + ".annot", prefix}, "annotation file");
}

std::string resolve_frq_file(const std::string& prefix) {
    return resolve_existing({prefix + ".frq", prefix}, "frequency file");
}

std::vector<int> present_annot_chrs(const std::string& prefix) {
    std::vector<int> out;
    for (int chr = 1; chr <= 22; ++chr) {
        const std::string p = sub_chr(prefix, chr);
        if (fs::exists(p + ".annot") || fs::exists(p)) out.push_back(chr);
    }
    if (out.empty()) throw Error("No chromosome-split .annot files found for prefix: " + prefix);
    return out;
}

Matrix read_annot_numeric_matrix(const std::string& annot_path, const std::optional<std::string>& frq_path) {
    Table annot = read_table(annot_path);
    std::vector<size_t> annot_cols;
    for (size_t i = 0; i < annot.columns.size(); ++i) {
        const std::string& c = annot.columns[i];
        if (c == "SNP" || c == "CHR" || c == "BP" || c == "CM") continue;
        annot_cols.push_back(i);
    }
    if (annot_cols.empty()) throw Error("Annotation file has no annotation columns: " + annot_path);

    std::vector<bool> keep(annot.rows.size(), true);
    if (frq_path) {
        Table frq = read_table(*frq_path);
        int frq_i = frq.col("FRQ");
        if (frq_i < 0) frq_i = frq.col("MAF");
        if (frq_i < 0) throw Error("Frequency file must contain FRQ or MAF: " + *frq_path);
        if (frq.rows.size() != annot.rows.size()) {
            throw Error("Annotation and frequency files have different row counts: " + annot_path + " / " + *frq_path);
        }
        for (size_t r = 0; r < frq.rows.size(); ++r) {
            const double f = parse_double(frq.rows[r][static_cast<size_t>(frq_i)]);
            keep[r] = f > 0.05 && f < 0.95;
        }
    }

    size_t kept = 0;
    for (bool k : keep) {
        if (k) ++kept;
    }
    Matrix out(kept, annot_cols.size());
    size_t dst = 0;
    for (size_t r = 0; r < annot.rows.size(); ++r) {
        if (!keep[r]) continue;
        for (size_t c = 0; c < annot_cols.size(); ++c) {
            out(dst, c) = parse_double(annot.rows[r][annot_cols[c]]);
        }
        ++dst;
    }
    return out;
}

Matrix hstack_matrices(const std::vector<Matrix>& parts) {
    if (parts.empty()) throw Error("No annotation matrices to concatenate.");
    const size_t rows = parts.front().rows;
    size_t cols = 0;
    for (const Matrix& m : parts) {
        if (m.rows != rows) throw Error("Annotation files have different row counts after frequency filtering.");
        cols += m.cols;
    }
    Matrix out(rows, cols);
    size_t offset = 0;
    for (const Matrix& m : parts) {
        for (size_t r = 0; r < rows; ++r) {
            for (size_t c = 0; c < m.cols; ++c) out(r, offset + c) = m(r, c);
        }
        offset += m.cols;
    }
    return out;
}

void add_overlap(Matrix& overlap, const Matrix& annot) {
    if (overlap.rows == 0) overlap = Matrix(annot.cols, annot.cols, 0.0);
    if (overlap.rows != annot.cols || overlap.cols != annot.cols) throw Error("Annotation column count mismatch.");
    for (size_t r = 0; r < annot.rows; ++r) {
        for (size_t i = 0; i < annot.cols; ++i) {
            for (size_t j = 0; j < annot.cols; ++j) overlap(i, j) += annot(r, i) * annot(r, j);
        }
    }
}

OverlapAnnotData read_overlap_annotations(const Args& args, Logger& log) {
    const bool chr_split = args.ref_ld_chr.has_value();
    const std::string ref_arg = chr_split ? *args.ref_ld_chr : args.ref_ld.value_or("");
    std::vector<std::string> prefixes = split_char(ref_arg, ',');
    OverlapAnnotData out;
    log.log("Reading annot matrix for overlapping categories ...");
    if (chr_split) {
        const std::vector<int> chrs = present_annot_chrs(prefixes.front());
        for (int chr : chrs) {
            std::optional<std::string> frq_path;
            if (args.frqfile_chr) frq_path = resolve_frq_file(sub_chr(*args.frqfile_chr, chr));
            std::vector<Matrix> parts;
            for (const std::string& prefix : prefixes) {
                parts.push_back(read_annot_numeric_matrix(resolve_annot_file(sub_chr(prefix, chr)), frq_path));
            }
            Matrix combined = hstack_matrices(parts);
            add_overlap(out.overlap, combined);
            out.m_tot += static_cast<double>(combined.rows);
        }
    } else {
        std::optional<std::string> frq_path;
        if (args.frqfile) frq_path = resolve_frq_file(*args.frqfile);
        std::vector<Matrix> parts;
        for (const std::string& prefix : prefixes) {
            parts.push_back(read_annot_numeric_matrix(resolve_annot_file(prefix), frq_path));
        }
        Matrix combined = hstack_matrices(parts);
        add_overlap(out.overlap, combined);
        out.m_tot = static_cast<double>(combined.rows);
    }
    return out;
}

std::string overlap_cell(double x) {
    return format_num(x, 17);
}

void write_overlap_results(const std::string& path, const std::vector<std::string>& names, const HsqResult& hsq,
                           const Matrix& overlap, const std::vector<double>& m_annot, double m_tot,
                           bool print_coefficients, int n_blocks) {
    const size_t p = hsq.n_annot;
    if (overlap.rows != p || overlap.cols != p || m_annot.size() != p || names.size() != p) {
        throw Error("Overlap annotation dimensions do not match regression annotations.");
    }
    std::ofstream out(path);
    if (!out) throw Error("Could not open overlap results file: " + path);
    out << std::setprecision(17);
    out << "Category\tProp._SNPs\tProp._h2\tProp._h2_std_error\tEnrichment\tEnrichment_std_error\tEnrichment_p";
    if (print_coefficients) {
        out << "\tCoefficient\tCoefficient_std_error\tCoefficient_z-score";
    }
    out << "\n";

    Matrix overlap_prop(p, p);
    for (size_t i = 0; i < p; ++i) {
        for (size_t j = 0; j < p; ++j) overlap_prop(i, j) = overlap(i, j) / m_annot[j];
    }

    Matrix overlap_diff(p, p, 0.0);
    for (size_t i = 0; i < p; ++i) {
        if (m_tot == m_annot[i]) continue;
        for (size_t j = 0; j < p; ++j) {
            overlap_diff(i, j) = overlap(i, j) / m_annot[i] -
                                 (m_annot[j] - overlap(i, j)) / (m_tot - m_annot[i]);
        }
    }

    for (size_t i = 0; i < p; ++i) {
        double prop_h2 = 0.0;
        for (size_t j = 0; j < p; ++j) prop_h2 += overlap_prop(i, j) * hsq.prop[j];
        double prop_h2_var = 0.0;
        for (size_t j = 0; j < p; ++j) {
            for (size_t k = 0; k < p; ++k) {
                prop_h2_var += overlap_prop(i, j) * hsq.prop_cov(j, k) * overlap_prop(i, k);
            }
        }
        const double prop_h2_se = std::sqrt(std::max(0.0, prop_h2_var));
        const double prop_m = m_annot[i] / m_tot;
        const double enrichment = prop_h2 / prop_m;
        const double enrichment_se = prop_h2_se / prop_m;

        double diff_est = 0.0;
        for (size_t j = 0; j < p; ++j) diff_est += overlap_diff(i, j) * hsq.coef[j];
        double diff_var = 0.0;
        for (size_t j = 0; j < p; ++j) {
            for (size_t k = 0; k < p; ++k) {
                diff_var += overlap_diff(i, j) * hsq.coef_cov(j, k) * overlap_diff(i, k);
            }
        }
        const double diff_se = std::sqrt(std::max(0.0, diff_var));

        out << names[i] << '\t' << overlap_cell(prop_m) << '\t' << overlap_cell(prop_h2) << '\t'
            << overlap_cell(prop_h2_se) << '\t' << overlap_cell(enrichment) << '\t'
            << overlap_cell(enrichment_se) << '\t';
        if (diff_se == 0.0) {
            out << "NA";
        } else {
            out << overlap_cell(student_t_two_sided_p(std::abs(diff_est / diff_se), n_blocks));
        }
        if (print_coefficients) {
            out << '\t' << overlap_cell(hsq.coef[i]) << '\t' << overlap_cell(hsq.coef_se[i]) << '\t'
                << overlap_cell(hsq.coef[i] / hsq.coef_se[i]);
        }
        out << "\n";
    }
}

void estimate_h2(const Args& args, Logger& log) {
    if (!args.h2) throw Error("--h2 is required.");
    if (args.overlap_annot && !args.not_m_5_50) {
        if (!((args.frqfile && args.ref_ld) || (args.frqfile_chr && args.ref_ld_chr))) {
            throw Error("Must set either --frqfile and --ref-ld or --frqfile-chr and --ref-ld-chr with --overlap-annot.");
        }
    }
    const bool ref_chr = args.ref_ld_chr.has_value();
    const bool w_chr = args.w_ld_chr.has_value();
    const std::string ref_arg = ref_chr ? *args.ref_ld_chr : args.ref_ld.value_or("");
    const std::string w_arg = w_chr ? *args.w_ld_chr : args.w_ld.value_or("");
    if (ref_arg.empty() || w_arg.empty()) throw Error("--ref-ld/--w-ld are required.");
    std::vector<double> m;
    BaseData data = read_ld_sumstats(*args.h2, ref_arg, ref_chr, w_arg, w_chr, args.manual_m,
                                     !args.not_m_5_50, m, log, false, true);
    const size_t n0 = data.snp.size();
    if (n0 < 200000) log.log("WARNING: number of SNPs less than 200k; this is almost always bad.");
    int n_blocks = std::min<int>(args.n_blocks, static_cast<int>(n0));
    std::optional<double> intercept = args.no_intercept ? std::optional<double>(1.0) : args.intercept_h2_scalar;
    std::optional<double> twostep = args.twostep;
    double chisq_max = args.chisq_max.value_or(kNan);
    bool has_chisq_max = args.chisq_max.has_value();
    bool old_weights = false;
    if (data.ld.cols == 1) {
        if (!twostep && !intercept) twostep = 30.0;
    } else {
        old_weights = true;
        if (!has_chisq_max) {
            const double max_n = *std::max_element(data.n.begin(), data.n.end());
            chisq_max = std::max(0.001 * max_n, 80.0);
            has_chisq_max = true;
        }
    }
    std::vector<double> chisq(data.z.size());
    for (size_t i = 0; i < data.z.size(); ++i) chisq[i] = data.z[i] * data.z[i];
    if (has_chisq_max) {
        std::vector<int> keep;
        for (size_t i = 0; i < chisq.size(); ++i) {
            if (chisq[i] < chisq_max) keep.push_back(static_cast<int>(i));
        }
        log.log("Removed " + std::to_string(chisq.size() - keep.size()) + " SNPs with chi^2 > " +
                format_py_float(chisq_max) + " (" + std::to_string(keep.size()) + " SNPs remain)");
        data.ld = data.ld.select_rows(keep);
        data.w_ld = select_vec(data.w_ld, keep);
        data.n = select_vec(data.n, keep);
        chisq = select_vec(chisq, keep);
    }
    if (twostep) log.log("Using two-step estimator with cutoff at " + format_num(*twostep) + ".");
    HsqResult hsq = run_hsq(chisq, data.ld, data.w_ld, data.n, m, n_blocks, intercept, twostep, old_weights);
    if (args.print_cov) {
        const std::string path = args.out + ".cov";
        log.log("Printing covariance matrix of the estimates to " + path + ".");
        save_matrix(path, hsq.coef_cov);
    }
    if (args.print_delete_vals) {
        const std::string delete_path = args.out + ".delete";
        const std::string part_delete_path = args.out + ".part_delete";
        log.log("Printing block jackknife delete values to " + delete_path + ".");
        save_matrix(delete_path, hsq.tot_delete_values);
        log.log("Printing partitioned block jackknife delete values to " + part_delete_path + ".");
        save_matrix(part_delete_path, hsq.part_delete_values);
    }
    std::optional<double> sp, pp;
    if (args.samp_prev && args.pop_prev) {
        sp = parse_double(*args.samp_prev);
        pp = parse_double(*args.pop_prev);
    }
    log.log(hsq_summary(hsq, data.ld_names, sp, pp, args.overlap_annot));
    if (args.overlap_annot) {
        OverlapAnnotData overlap = read_overlap_annotations(args, log);
        write_overlap_results(args.out + ".results", data.ld_names, hsq, overlap.overlap, m, overlap.m_tot,
                              args.print_coefficients, n_blocks);
        log.log("Results printed to " + args.out + ".results");
    }
}

void estimate_h2_cts(const Args& args, Logger& log) {
    if (!args.h2_cts) throw Error("--h2-cts is required.");
    if (!args.ref_ld_chr || !args.ref_ld_chr_cts || !(args.w_ld || args.w_ld_chr)) {
        throw Error("--h2-cts requires --ref-ld-chr, --w-ld/--w-ld-chr, and --ref-ld-chr-cts.");
    }
    const bool w_chr = args.w_ld_chr.has_value();
    const std::string w_arg = w_chr ? *args.w_ld_chr : args.w_ld.value_or("");
    std::vector<double> base_m;
    BaseData base = read_ld_sumstats(*args.h2_cts, *args.ref_ld_chr, true, w_arg, w_chr, args.manual_m,
                                     !args.not_m_5_50, base_m, log, false, true);
    if (base.snp.size() < 200000) log.log("WARNING: number of SNPs less than 200k; this is almost always bad.");
    double chisq_max = args.chisq_max.value_or(std::max(0.001 * (*std::max_element(base.n.begin(), base.n.end())), 80.0));
    std::vector<double> chisq(base.z.size());
    for (size_t i = 0; i < base.z.size(); ++i) chisq[i] = base.z[i] * base.z[i];
    std::vector<int> keep;
    for (size_t i = 0; i < chisq.size(); ++i) {
        if (chisq[i] < chisq_max) keep.push_back(static_cast<int>(i));
    }
    log.log("Removed " + std::to_string(chisq.size() - keep.size()) + " SNPs with chi^2 > " +
            format_py_float(chisq_max) + " (" + std::to_string(keep.size()) + " SNPs remain)");
    base.ld = base.ld.select_rows(keep);
    base.w_ld = select_vec(base.w_ld, keep);
    base.n = select_vec(base.n, keep);
    chisq = select_vec(chisq, keep);
    std::vector<std::string> keep_snps;
    for (int idx : keep) keep_snps.push_back(base.snp[static_cast<size_t>(idx)]);

    struct CtsResult { std::string name; double coef; double se; double p; };
    std::vector<CtsResult> results;
    std::ifstream list(*args.ref_ld_chr_cts);
    if (!list) throw Error("Could not open --ref-ld-chr-cts: " + *args.ref_ld_chr_cts);
    std::string line;
    while (std::getline(list, line)) {
        line = trim(line);
        if (line.empty()) continue;
        std::vector<std::string> parts = split_ws(line);
        if (parts.size() < 2) throw Error("Each --ref-ld-chr-cts line must contain name and LD prefix.");
        const std::string name = parts[0];
        const std::string ct_prefix = parts[1];
        LdFrame cts_all = read_ldscore_list(ct_prefix, true);
        std::unordered_map<std::string, size_t> cts_idx = index_by_snp(cts_all.snp);
        Matrix cts_ld(keep_snps.size(), cts_all.values.cols);
        for (size_t i = 0; i < keep_snps.size(); ++i) {
            auto it = cts_idx.find(keep_snps[i]);
            if (it == cts_idx.end()) throw Error("Missing some LD scores from cts files.");
            for (size_t c = 0; c < cts_all.values.cols; ++c) cts_ld(i, c) = cts_all.values(it->second, c);
        }
        std::vector<double> cts_m = read_m(std::nullopt, ct_prefix, true, !args.not_m_5_50);
        Matrix all_ld(keep_snps.size(), cts_ld.cols + base.ld.cols);
        for (size_t r = 0; r < all_ld.rows; ++r) {
            for (size_t c = 0; c < cts_ld.cols; ++c) all_ld(r, c) = cts_ld(r, c);
            for (size_t c = 0; c < base.ld.cols; ++c) all_ld(r, cts_ld.cols + c) = base.ld(r, c);
        }
        std::vector<double> all_m = cts_m;
        all_m.insert(all_m.end(), base_m.begin(), base_m.end());
        std::optional<double> intercept = args.no_intercept ? std::optional<double>(1.0) : args.intercept_h2_scalar;
        HsqResult h = run_hsq(chisq, all_ld, base.w_ld, base.n, all_m,
                              std::min<int>(args.n_blocks, static_cast<int>(keep_snps.size())),
                              intercept, std::nullopt, true);
        auto add_result = [&](size_t idx, const std::string& label) {
            double z = h.coef[idx] / h.coef_se[idx];
            double p = 0.5 * std::erfc(z / std::sqrt(2.0));
            results.push_back({label, h.coef[idx], h.coef_se[idx], p});
        };
        add_result(0, name);
        if (args.print_all_cts) {
            for (size_t i = 1; i < cts_ld.cols; ++i) add_result(i, name + "_" + std::to_string(i));
        }
    }
    std::sort(results.begin(), results.end(), [](const CtsResult& a, const CtsResult& b) { return a.p < b.p; });
    std::ofstream out(args.out + ".cell_type_results.txt");
    if (!out) throw Error("Could not open cell type output.");
    out << "Name\tCoefficient\tCoefficient_std_error\tCoefficient_P_value\n";
    out << std::setprecision(10);
    for (const auto& r : results) out << r.name << '\t' << r.coef << '\t' << r.se << '\t' << r.p << "\n";
    log.log("Results printed to " + args.out + ".cell_type_results.txt");
}

RGMergedData merge_other_sumstats(const BaseData& base, const std::string& path, bool no_check_alleles, Logger& log) {
    log.log("Reading summary statistics from " + path + " ...");
    std::vector<SumRecord> other = read_sumstats(path, true, false);
    log.log("Read summary statistics for " + std::to_string(other.size()) + " SNPs.");
    std::unordered_map<std::string, size_t> other_i;
    for (size_t i = 0; i < other.size(); ++i) other_i.emplace(other[i].snp, i);

    std::vector<int> keep;
    std::vector<double> z2_aligned;
    size_t merged = 0;
    for (size_t i = 0; i < base.snp.size(); ++i) {
        auto it = other_i.find(base.snp[i]);
        if (it == other_i.end()) continue;
        ++merged;
        const SumRecord& o = other[it->second];
        if (std::isnan(o.z) || std::isnan(o.n) || o.a1 == "." || o.a2 == ".") continue;
        double z2 = o.z;
        if (!no_check_alleles) {
            const std::string allele = base.a1[i] + base.a2[i] + o.a1 + o.a2;
            if (!alleles_match(allele)) continue;
            if (flip_alleles(allele)) z2 = -z2;
        }
        keep.push_back(static_cast<int>(i));
        z2_aligned.push_back(z2);
    }
    log.log("After merging with summary statistics, " + std::to_string(merged) + " SNPs remain.");
    if (!no_check_alleles) log.log(std::to_string(keep.size()) + " SNPs with valid alleles.");
    RGMergedData out;
    out.ld = base.ld.select_rows(keep);
    out.z1 = select_vec(base.z, keep);
    out.n1 = select_vec(base.n, keep);
    out.w_ld = select_vec(base.w_ld, keep);
    out.z2 = std::move(z2_aligned);
    out.n2.reserve(keep.size());
    for (int idx : keep) {
        const auto it = other_i.find(base.snp[static_cast<size_t>(idx)]);
        out.n2.push_back(other[it->second].n);
    }
    return out;
}

void estimate_rg(const Args& args, Logger& log) {
    if (!args.rg) throw Error("--rg is required.");
    std::vector<std::string> paths = split_char(*args.rg, ',');
    if (paths.size() < 2) throw Error("Must specify at least two phenotypes for rg estimation.");
    const bool ref_chr = args.ref_ld_chr.has_value();
    const bool w_chr = args.w_ld_chr.has_value();
    const std::string ref_arg = ref_chr ? *args.ref_ld_chr : args.ref_ld.value_or("");
    const std::string w_arg = w_chr ? *args.w_ld_chr : args.w_ld.value_or("");
    if (ref_arg.empty() || w_arg.empty()) throw Error("--ref-ld/--w-ld are required.");

    std::vector<std::optional<double>> int_h2 = split_optional_doubles(args.intercept_h2_list, paths.size());
    std::vector<std::optional<double>> int_g = split_optional_doubles(args.intercept_gencov_list, paths.size());
    std::vector<std::optional<double>> samp_prev = split_optional_doubles(args.samp_prev, paths.size());
    std::vector<std::optional<double>> pop_prev = split_optional_doubles(args.pop_prev, paths.size());
    if (args.no_intercept) {
        for (auto& v : int_h2) v = 1.0;
        for (auto& v : int_g) v = 0.0;
    }

    std::vector<double> m;
    BaseData base = read_ld_sumstats(paths[0], ref_arg, ref_chr, w_arg, w_chr, args.manual_m,
                                     !args.not_m_5_50, m, log, true, true);
    std::optional<double> twostep = args.twostep;
    if (twostep) log.log("Using two-step estimator with cutoff at " + format_num(*twostep) + ".");

    std::vector<RGResult> results;
    for (size_t i = 1; i < paths.size(); ++i) {
        log.log("Computing rg for phenotype " + std::to_string(i + 1) + "/" + std::to_string(paths.size()));
        RGMergedData loop = merge_other_sumstats(base, paths[i], args.no_check_alleles, log);
        if (loop.z1.size() < 200000) log.log("WARNING: number of SNPs less than 200k; this is almost always bad.");
        if (args.chisq_max) {
            std::vector<int> keep;
            for (size_t r = 0; r < loop.z1.size(); ++r) {
                if (loop.z1[r] * loop.z1[r] * loop.z2[r] * loop.z2[r] < (*args.chisq_max) * (*args.chisq_max)) {
                    keep.push_back(static_cast<int>(r));
                }
            }
            loop.ld = loop.ld.select_rows(keep);
            loop.z1 = select_vec(loop.z1, keep);
            loop.z2 = select_vec(loop.z2, keep);
            loop.n1 = select_vec(loop.n1, keep);
            loop.n2 = select_vec(loop.n2, keep);
            loop.w_ld = select_vec(loop.w_ld, keep);
        }
        int n_blocks = std::min<int>(args.n_blocks, static_cast<int>(loop.z1.size()));
        RGResult rg = run_rg_model(loop.z1, loop.z2, loop.ld, loop.w_ld, loop.n1, loop.n2, m, int_h2[0], int_h2[i],
                                   int_g[i], n_blocks, twostep);
        results.push_back(rg);

        if (i == 1) {
            log.log(py2_section("Heritability of phenotype 1"));
            log.log(hsq_summary(rg.hsq1, base.ld_names, samp_prev[0], pop_prev[0]));
        }
        log.log(py2_section("Heritability of phenotype " + std::to_string(i + 1) + "/" + std::to_string(paths.size())));
        log.log(hsq_summary(rg.hsq2, base.ld_names, samp_prev[i], pop_prev[i]));
        log.log(py2_section("Genetic Covariance"));
        log.log(gencov_summary(rg.gencov, base.ld_names, {samp_prev[0], samp_prev[i]}, {pop_prev[0], pop_prev[i]}));
        log.log(py2_section("Genetic Correlation"));
        log.log(rg_summary(rg, args.return_silly_things) + "\n");

        std::string out_loop = args.out + fs::path(paths[0]).filename().string() + "_" + fs::path(paths[i]).filename().string();
        if (args.print_cov) {
            log.log("Printing covariance matrix of the estimates to " + out_loop + ".hsq1.cov.");
            save_matrix(out_loop + ".hsq1.cov", rg.hsq1.coef_cov);
            log.log("Printing covariance matrix of the estimates to " + out_loop + ".hsq2.cov.");
            save_matrix(out_loop + ".hsq2.cov", rg.hsq2.coef_cov);
            log.log("Printing covariance matrix of the estimates to " + out_loop + ".gencov.cov.");
            save_matrix(out_loop + ".gencov.cov", rg.gencov.coef_cov);
        }
        if (args.print_delete_vals) {
            log.log("Printing block jackknife delete values to " + out_loop + ".hsq1.delete.");
            save_matrix(out_loop + ".hsq1.delete", rg.hsq1.tot_delete_values);
            log.log("Printing block jackknife delete values to " + out_loop + ".hsq2.delete.");
            save_matrix(out_loop + ".hsq2.delete", rg.hsq2.tot_delete_values);
            log.log("Printing block jackknife delete values to " + out_loop + ".gencov.delete.");
            save_matrix(out_loop + ".gencov.delete", rg.gencov.tot_delete_values);
        }
    }

    log.log("\nSummary of Genetic Correlation Results");
    log.log("p1                   p2      rg      se         z    p  h2_obs  h2_obs_se  h2_int  h2_int_se  gcov_int  gcov_int_se");
    for (size_t i = 0; i < results.size(); ++i) {
        log.log(rg_table_row(paths[0], paths[i + 1], results[i]));
    }
    log.log("");
}

}  // namespace ldsc_cpp
