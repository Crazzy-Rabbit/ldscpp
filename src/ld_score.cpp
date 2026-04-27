#include "ld_score.hpp"

namespace ldsc_cpp {

std::vector<BimRecord> read_bim_records(const std::string& path) {
    std::vector<std::vector<std::string>> rows = read_rows_no_header(path);
    std::vector<BimRecord> out;
    for (const auto& row : rows) {
        if (row.size() < 6) throw Error("BIM file must have 6 columns.");
        BimRecord r;
        r.chr = row[0];
        r.snp = row[1];
        r.cm = parse_double(row[2]);
        r.bp = static_cast<int>(parse_double(row[3]));
        r.a1 = row[4];
        r.a2 = row[5];
        out.push_back(r);
    }
    return out;
}

std::vector<std::string> read_fam_iids(const std::string& path) {
    std::vector<std::vector<std::string>> rows = read_rows_no_header(path);
    std::vector<std::string> out;
    for (const auto& row : rows) {
        if (row.size() < 2) throw Error("FAM file must have at least 2 columns.");
        out.push_back(row[1]);
    }
    return out;
}

std::set<std::string> read_one_col_set(const std::string& path) {
    std::vector<std::vector<std::string>> rows = read_rows_no_header(path);
    std::set<std::string> out;
    for (const auto& row : rows) {
        if (!row.empty()) out.insert(row[0]);
    }
    return out;
}

GenotypeData read_plink_bed(const std::string& prefix, const Args& args, Logger& log, bool log_messages = true) {
    std::vector<BimRecord> bim = read_bim_records(prefix + ".bim");
    std::vector<std::string> fam = read_fam_iids(prefix + ".fam");
    if (log_messages) {
        log.log("Read list of " + std::to_string(bim.size()) + " SNPs from " + prefix + ".bim");
        log.log("Read list of " + std::to_string(fam.size()) + " individuals from " + prefix + ".fam");
        log.log("Reading genotypes from " + prefix + ".bed");
    }
    std::vector<int> keep_indivs;
    if (args.keep) {
        std::set<std::string> keep_ids = read_one_col_set(*args.keep);
        for (size_t i = 0; i < fam.size(); ++i) {
            if (keep_ids.count(fam[i])) keep_indivs.push_back(static_cast<int>(i));
        }
    } else {
        for (size_t i = 0; i < fam.size(); ++i) keep_indivs.push_back(static_cast<int>(i));
    }
    std::set<std::string> keep_snps;
    if (args.extract) keep_snps = read_one_col_set(*args.extract);

    std::ifstream in(prefix + ".bed", std::ios::binary);
    if (!in) throw Error("Could not open PLINK BED: " + prefix + ".bed");
    unsigned char magic[3];
    in.read(reinterpret_cast<char*>(magic), 3);
    if (magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01) {
        throw Error("Plink .bed file must be SNP-major with recognized magic number.");
    }
    const size_t n_raw = fam.size();
    const size_t bytes_per_snp = (n_raw + 3) / 4;
    std::vector<std::vector<double>> kept_cols;
    GenotypeData out;
    for (size_t s = 0; s < bim.size(); ++s) {
        std::vector<unsigned char> bytes(bytes_per_snp);
        in.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
        if (!in) throw Error("Unexpected end of BED file.");
        if (!keep_snps.empty() && !keep_snps.count(bim[s].snp)) continue;
        std::vector<double> geno;
        geno.reserve(keep_indivs.size());
        for (int idx : keep_indivs) {
            const unsigned char byte = bytes[static_cast<size_t>(idx) / 4];
            const int code = (byte >> ((idx % 4) * 2)) & 0x3;
            double g = 9.0;
            if (code == 0) g = 0.0;
            else if (code == 2) g = 1.0;
            else if (code == 3) g = 2.0;
            geno.push_back(g);
        }
        double sum = 0.0;
        int nonmiss = 0;
        int het_or_miss = 0;
        for (double g : geno) {
            if (g == 9.0) {
                ++het_or_miss;
                continue;
            }
            if (g == 1.0) ++het_or_miss;
            sum += g;
            ++nonmiss;
        }
        double f = nonmiss > 0 ? sum / (2.0 * nonmiss) : 0.0;
        double maf = std::min(f, 1.0 - f);
        if (!(maf > args.maf) || het_or_miss == static_cast<int>(geno.size())) continue;
        double avg = nonmiss > 0 ? sum / nonmiss : 0.0;
        for (double& g : geno) {
            if (g == 9.0) g = avg;
        }
        double var = 0.0;
        for (double g : geno) var += (g - avg) * (g - avg);
        var /= static_cast<double>(geno.size());
        double sd = std::sqrt(var);
        if (sd == 0.0) sd = 1.0;
        for (double& g : geno) g = (g - avg) / sd;
        kept_cols.push_back(std::move(geno));
        out.records.push_back(bim[s]);
        out.kept_indices.push_back(s);
        out.freq.push_back(f);
        out.maf.push_back(maf);
    }
    out.x = Matrix(keep_indivs.size(), kept_cols.size());
    for (size_t j = 0; j < kept_cols.size(); ++j) {
        for (size_t i = 0; i < kept_cols[j].size(); ++i) out.x(i, j) = kept_cols[j][i];
    }
    std::cout << "After filtering, " << out.records.size() << " SNPs remain\n";
    return out;
}

std::vector<std::string> read_one_col_list(const std::string& path) {
    std::vector<std::vector<std::string>> rows = read_rows_no_header(path);
    std::vector<std::string> out;
    for (const auto& row : rows) {
        if (!row.empty()) out.push_back(row[0]);
    }
    return out;
}

double ldscore_quantile(std::vector<double> x, double q) {
    if (x.empty()) return kNan;
    std::sort(x.begin(), x.end());
    const double pos = (static_cast<double>(x.size()) - 1.0) * q;
    const size_t lo = static_cast<size_t>(std::floor(pos));
    const size_t hi = static_cast<size_t>(std::ceil(pos));
    if (lo == hi) return x[lo];
    return x[lo] + (x[hi] - x[lo]) * (pos - static_cast<double>(lo));
}

double sample_std(const std::vector<double>& x) {
    if (x.size() < 2) return kNan;
    const double mu = mean(x);
    double ss = 0.0;
    for (double v : x) ss += (v - mu) * (v - mu);
    if (ss < 1e-24) return 0.0;
    return std::sqrt(ss / static_cast<double>(x.size() - 1));
}

double corr(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size() || a.size() < 2) return kNan;
    const double ma = mean(a);
    const double mb = mean(b);
    double va = 0.0, vb = 0.0, cov = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        const double da = a[i] - ma;
        const double db = b[i] - mb;
        va += da * da;
        vb += db * db;
        cov += da * db;
    }
    if (va < 1e-24 || vb < 1e-24) return kNan;
    return cov / std::sqrt(va * vb);
}

std::string pandas_summary_num(double x) {
    if (std::isnan(x)) return "NaN";
    std::ostringstream out;
    out << std::fixed << std::setprecision(4) << x;
    return out.str();
}

std::string pandas_corr_num(double x) {
    if (std::isnan(x)) return "NaN";
    if (std::abs(x - 1.0) < 5e-13) return "1.0";
    return format_num(x, 4);
}

void log_ldscore_summary(Logger& log, const std::string& path, const std::vector<double>& maf,
                         const Matrix& ld, const std::vector<std::string>& names) {
    std::vector<std::string> colnames{"MAF"};
    for (const auto& n : names) colnames.push_back(n);
    std::vector<std::vector<double>> cols;
    cols.push_back(maf);
    for (size_t c = 0; c < ld.cols; ++c) {
        std::vector<double> v(ld.rows);
        for (size_t r = 0; r < ld.rows; ++r) v[r] = ld(r, c);
        cols.push_back(std::move(v));
    }

    log.log("\nSummary of LD Scores in " + path);
    const std::vector<std::string> rows{"mean", "std", "min", "25%", "50%", "75%", "max"};
    std::ostringstream out;
    out << std::setw(4) << "";
    for (size_t c = 0; c < colnames.size(); ++c) {
        out << (c == 0 ? "  " : "  ") << std::setw(6) << colnames[c];
    }
    out << "\n";
    for (const auto& row : rows) {
        out << std::left << std::setw(4) << row << std::right;
        for (const auto& values : cols) {
            double v = kNan;
            if (row == "mean") v = mean(values);
            else if (row == "std") v = sample_std(values);
            else if (row == "min") v = *std::min_element(values.begin(), values.end());
            else if (row == "25%") v = ldscore_quantile(values, 0.25);
            else if (row == "50%") v = ldscore_quantile(values, 0.50);
            else if (row == "75%") v = ldscore_quantile(values, 0.75);
            else if (row == "max") v = *std::max_element(values.begin(), values.end());
            out << "  " << std::setw(6) << pandas_summary_num(v);
        }
        if (row != rows.back()) out << "\n";
    }
    log.log(out.str());

    log.log("\nMAF/LD Score Correlation Matrix");
    std::ostringstream cor;
    cor << std::setw(3) << "";
    for (size_t c = 0; c < colnames.size(); ++c) {
        cor << (c == 0 ? "  " : " ") << std::setw(3) << colnames[c];
    }
    cor << "\n";
    for (size_t i = 0; i < cols.size(); ++i) {
        cor << std::left << std::setw(3) << colnames[i] << std::right;
        for (size_t j = 0; j < cols.size(); ++j) {
            cor << (j == 0 ? "  " : " ") << std::setw(3) << pandas_corr_num(corr(cols[i], cols[j]));
        }
        if (i + 1 < cols.size()) cor << "\n";
    }
    log.log(cor.str());
}

std::vector<int> get_block_lefts(const std::vector<double>& coords, double max_dist) {
    std::vector<int> out(coords.size(), 0);
    size_t j = 0;
    for (size_t i = 0; i < coords.size(); ++i) {
        while (j < coords.size() && std::abs(coords[j] - coords[i]) > max_dist) ++j;
        out[i] = static_cast<int>(j);
    }
    return out;
}

double unbiased_ld_corr(const GenotypeData& geno, size_t j, size_t k) {
    double r = 0.0;
    for (size_t i = 0; i < geno.x.rows; ++i) r += geno.x(i, j) * geno.x(i, k);
    r /= static_cast<double>(geno.x.rows);
    const double r2 = r * r;
    const double denom = geno.x.rows > 2 ? static_cast<double>(geno.x.rows - 2) : static_cast<double>(geno.x.rows);
    return r2 - (1.0 - r2) / denom;
}

void add_ld_block(Matrix& ld, const Matrix& annot, const GenotypeData& geno,
                  size_t row_begin, size_t row_count, size_t col_begin, size_t col_count) {
    for (size_t rr = 0; rr < row_count; ++rr) {
        const size_t j = row_begin + rr;
        for (size_t cc = 0; cc < col_count; ++cc) {
            const size_t k = col_begin + cc;
            const double value = unbiased_ld_corr(geno, j, k);
            for (size_t a = 0; a < annot.cols; ++a) ld(j, a) += value * annot(k, a);
        }
    }
}

void add_ld_cross_block(Matrix& ld, const Matrix& annot, const GenotypeData& geno,
                        size_t left_begin, size_t left_count, size_t right_begin, size_t right_count) {
    for (size_t rr = 0; rr < left_count; ++rr) {
        const size_t j = left_begin + rr;
        for (size_t cc = 0; cc < right_count; ++cc) {
            const size_t k = right_begin + cc;
            const double value = unbiased_ld_corr(geno, j, k);
            for (size_t a = 0; a < annot.cols; ++a) {
                ld(j, a) += value * annot(k, a);
                ld(k, a) += value * annot(j, a);
            }
        }
    }
}

Matrix ldscore_var_blocks_py2(const GenotypeData& geno, const std::vector<int>& block_left,
                              int chunk_size, const Matrix& annot) {
    if (chunk_size <= 0) throw Error("--chunk-size must be a positive integer.");
    const size_t m = geno.records.size();
    Matrix ld(m, annot.cols, 0.0);
    if (m == 0) return ld;

    size_t c = static_cast<size_t>(chunk_size);
    std::vector<size_t> block_sizes(m, 0);
    for (size_t i = 0; i < m; ++i) {
        const size_t left = static_cast<size_t>(block_left[i]);
        const size_t width = i >= left ? i - left : 0;
        block_sizes[i] = static_cast<size_t>(std::ceil(static_cast<double>(width) / static_cast<double>(c))) * c;
    }

    size_t b = m;
    for (size_t i = 0; i < m; ++i) {
        if (block_left[i] > 0) {
            b = i;
            break;
        }
    }
    b = static_cast<size_t>(std::ceil(static_cast<double>(b) / static_cast<double>(c))) * c;
    if (b > m) {
        c = 1;
        b = m;
    }

    size_t l_a = 0;
    for (size_t l_b = 0; l_b < b; l_b += c) {
        add_ld_block(ld, annot, geno, l_a, b, l_b, c);
    }

    const size_t b0 = b;
    const size_t md = c * (m / c);
    const size_t end = md != m ? md + 1 : md;
    for (size_t l_b = b0; l_b < end; l_b += c) {
        const size_t old_b = b;
        b = block_sizes[l_b];
        if (l_b > b0 && b > 0) {
            l_a += old_b - b + c;
        } else if (l_b == b0 && b > 0) {
            l_a = b0 - b;
        } else if (b == 0) {
            l_a = l_b;
        }

        size_t current_c = c;
        if (l_b == md) current_c = m - md;
        if (current_c == 0) continue;

        const bool left_zero = b == 0;
        bool annot_left_zero = true;
        if (!left_zero) {
            for (size_t r = l_a; r < l_a + b && annot_left_zero; ++r) {
                for (size_t a = 0; a < annot.cols; ++a) {
                    if (annot(r, a) != 0.0) {
                        annot_left_zero = false;
                        break;
                    }
                }
            }
        }
        bool annot_right_zero = true;
        for (size_t r = l_b; r < l_b + current_c && annot_right_zero; ++r) {
            for (size_t a = 0; a < annot.cols; ++a) {
                if (annot(r, a) != 0.0) {
                    annot_right_zero = false;
                    break;
                }
            }
        }
        if (!(annot_left_zero && annot_right_zero)) {
            add_ld_cross_block(ld, annot, geno, l_a, b, l_b, current_c);
            add_ld_block(ld, annot, geno, l_b, current_c, l_b, current_c);
        }
    }
    return ld;
}

Matrix read_annot_matrix_for_ldscore(const Args& args, size_t expected_rows, std::vector<std::string>& names,
                                     Logger& log) {
    if (!args.annot) {
        names = {"L2"};
        return Matrix(expected_rows, 1, 1.0);
    }
    Table a = read_table(*args.annot);
    size_t offset = args.thin_annot ? 0 : 4;
    if (a.columns.size() <= offset) throw Error("Annotation file has no annotation columns.");
    names.assign(a.columns.begin() + static_cast<long long>(offset), a.columns.end());
    log.log("Read " + std::to_string(names.size()) + " annotations for " + std::to_string(a.rows.size()) +
            " SNPs from " + *args.annot);
    if (a.rows.size() != expected_rows) throw Error("Annotation file must contain same number of SNPs as filtered BIM.");
    Matrix out(a.rows.size(), names.size());
    for (size_t i = 0; i < a.rows.size(); ++i) {
        for (size_t j = 0; j < names.size(); ++j) out(i, j) = parse_double(a.rows[i][offset + j]);
    }
    return out;
}

Matrix make_cts_annot_matrix(const Args& args, const std::vector<BimRecord>& bim, std::vector<std::string>& names) {
    if (!args.cts_bin || !args.cts_breaks) throw Error("Must set both --cts-bin and --cts-breaks.");
    std::vector<std::string> files = split_char(*args.cts_bin, ',');
    std::string breaks_arg = *args.cts_breaks;
    std::replace(breaks_arg.begin(), breaks_arg.end(), 'N', '-');
    std::vector<std::string> break_groups = split_char(breaks_arg, 'x');
    if (files.size() != break_groups.size()) throw Error("Need one set of --cts-breaks for each --cts-bin file.");
    std::vector<std::string> cts_names;
    if (args.cts_names) cts_names = split_char(*args.cts_names, ',');
    else {
        for (size_t i = 0; i < files.size(); ++i) cts_names.push_back("ANNOT" + std::to_string(i));
    }
    if (cts_names.size() != files.size()) throw Error("--cts-names must have one value per --cts-bin file.");

    std::unordered_map<std::string, size_t> snp_i;
    for (size_t i = 0; i < bim.size(); ++i) snp_i[bim[i].snp] = i;
    std::vector<std::vector<int>> bin_by_var(files.size(), std::vector<int>(bim.size(), -1));
    std::vector<std::vector<std::string>> labels_by_var;
    std::vector<int> radix;
    for (size_t f = 0; f < files.size(); ++f) {
        std::vector<std::vector<std::string>> rows = read_rows_no_header(files[f]);
        std::vector<double> values(bim.size(), kNan);
        for (const auto& row : rows) {
            if (row.size() < 2) continue;
            auto it = snp_i.find(row[0]);
            if (it != snp_i.end()) values[it->second] = parse_double(row[1]);
        }
        std::vector<double> present;
        for (double v : values) {
            if (!std::isnan(v)) present.push_back(v);
        }
        if (present.empty()) throw Error("--cts-bin file has no SNPs matching the genotype data.");
        const double vmin = *std::min_element(present.begin(), present.end());
        const double vmax = *std::max_element(present.begin(), present.end());
        std::vector<double> cut;
        for (const auto& x : split_char(break_groups[f], ',')) cut.push_back(parse_double(x));
        bool all_ge_max = std::all_of(cut.begin(), cut.end(), [&](double x) { return x >= vmax; });
        bool all_le_min = std::all_of(cut.begin(), cut.end(), [&](double x) { return x <= vmin; });
        if (all_ge_max || all_le_min) throw Error("All breaks lie outside the range of the cts variable.");
        std::vector<std::string> name_breaks;
        std::vector<double> name_vals = cut;
        if (std::all_of(cut.begin(), cut.end(), [&](double x) { return x <= vmax; })) {
            cut.push_back(vmax + 1.0);
            name_vals.push_back(vmax);
        }
        if (std::all_of(cut.begin(), cut.end(), [&](double x) { return x >= vmin; })) {
            cut.push_back(vmin - 1.0);
            name_vals.push_back(vmin);
        }
        std::sort(cut.begin(), cut.end());
        std::sort(name_vals.begin(), name_vals.end());
        for (double x : name_vals) name_breaks.push_back(format_num(x));
        if (!name_breaks.empty()) {
            name_breaks.front() = "min";
            name_breaks.back() = "max";
        }
        std::vector<std::string> labels;
        for (size_t i = 0; i + 1 < name_breaks.size(); ++i) labels.push_back(name_breaks[i] + "_" + name_breaks[i + 1]);
        if (labels.empty()) throw Error("No cts bins were generated.");
        for (size_t i = 0; i < values.size(); ++i) {
            double v = values[i];
            if (std::isnan(v)) throw Error("--cts-bin is missing SNP " + bim[i].snp);
            for (size_t b = 0; b + 1 < cut.size(); ++b) {
                if (v > cut[b] && v <= cut[b + 1]) {
                    bin_by_var[f][i] = static_cast<int>(b);
                    break;
                }
            }
            if (bin_by_var[f][i] < 0) throw Error("Could not assign cts value to a bin.");
        }
        radix.push_back(static_cast<int>(labels.size()));
        labels_by_var.push_back(std::move(labels));
    }
    size_t n_cols = 1;
    for (int r : radix) n_cols *= static_cast<size_t>(r);
    Matrix out(bim.size(), n_cols, 0.0);
    names.assign(n_cols, "");
    for (size_t col = 0; col < n_cols; ++col) {
        size_t tmp = col;
        std::vector<int> parts(files.size(), 0);
        for (int f = static_cast<int>(files.size()) - 1; f >= 0; --f) {
            parts[static_cast<size_t>(f)] = static_cast<int>(tmp % static_cast<size_t>(radix[static_cast<size_t>(f)]));
            tmp /= static_cast<size_t>(radix[static_cast<size_t>(f)]);
        }
        std::ostringstream label;
        for (size_t f = 0; f < files.size(); ++f) {
            if (f) label << '_';
            label << cts_names[f] << '_' << labels_by_var[f][static_cast<size_t>(parts[f])];
        }
        names[col] = label.str();
    }
    for (size_t i = 0; i < bim.size(); ++i) {
        size_t col = 0;
        for (size_t f = 0; f < files.size(); ++f) {
            col = col * static_cast<size_t>(radix[f]) + static_cast<size_t>(bin_by_var[f][i]);
        }
        out(i, col) = 1.0;
    }
    return out;
}

std::vector<std::string> ldscore_colnames(const std::vector<std::string>& annot_names, size_t n_annot,
                                          const std::string& scale_suffix) {
    if (n_annot == 1) return {"L2" + scale_suffix};
    std::vector<std::string> out;
    out.reserve(annot_names.size());
    for (const auto& name : annot_names) out.push_back(name + "L2" + scale_suffix);
    return out;
}

void estimate_ldscore_cpp(const Args& args, Logger& log) {
    if (!args.bfile || !args.l2) throw Error("--bfile and --l2 are required for LD Score estimation.");
    int wind_count = (args.ld_wind_snps ? 1 : 0) + (args.ld_wind_kb ? 1 : 0) + (args.ld_wind_cm ? 1 : 0);
    if (wind_count != 1) throw Error("Must specify exactly one --ld-wind option.");
    if (args.annot && args.extract) throw Error("--annot and --extract are currently incompatible.");
    if (args.cts_bin && args.extract) throw Error("--cts-bin and --extract are currently incompatible.");
    if (args.annot && args.cts_bin) throw Error("--annot and --cts-bin are currently incompatible.");
    if (args.cts_bin.has_value() != args.cts_breaks.has_value()) throw Error("Must set both or neither of --cts-bin and --cts-breaks.");

    std::vector<BimRecord> bim = read_bim_records(*args.bfile + ".bim");
    log.log("Read list of " + std::to_string(bim.size()) + " SNPs from " + *args.bfile + ".bim");
    std::vector<std::string> annot_names;
    Matrix annot_full;
    bool have_full_annot = false;
    if (args.cts_bin) log.log("Reading numbers with which to bin SNPs from " + *args.cts_bin);
    if (args.cts_bin) {
        annot_full = make_cts_annot_matrix(args, bim, annot_names);
        have_full_annot = true;
    } else if (args.annot) {
        annot_full = read_annot_matrix_for_ldscore(args, bim.size(), annot_names, log);
        have_full_annot = true;
    }

    std::vector<std::string> fam = read_fam_iids(*args.bfile + ".fam");
    log.log("Read list of " + std::to_string(fam.size()) + " individuals from " + *args.bfile + ".fam");
    log.log("Reading genotypes from " + *args.bfile + ".bed");
    GenotypeData geno = read_plink_bed(*args.bfile, args, log, false);
    Matrix annot;
    if (have_full_annot) {
        std::vector<int> keep;
        keep.reserve(geno.kept_indices.size());
        for (size_t idx : geno.kept_indices) keep.push_back(static_cast<int>(idx));
        annot = annot_full.select_rows(keep);
    } else {
        annot_names = {"L2"};
        annot = Matrix(geno.records.size(), 1, 1.0);
    }
    if (args.per_allele && args.pq_exp) throw Error("Cannot set both --per-allele and --pq-exp.");
    std::optional<double> pq_exp = args.pq_exp;
    if (args.per_allele) pq_exp = 1.0;
    std::string scale_suffix;
    if (pq_exp) {
        log.log("Computing LD with pq ^ " + format_py_float(*pq_exp) + ".");
        log.log("Note that LD Scores with pq raised to a nonzero power arenot directly comparable to normal LD Scores.");
        scale_suffix = "_S" + format_py_float(*pq_exp);
    }
    if (pq_exp) {
        for (size_t i = 0; i < annot.rows; ++i) {
            double pq = std::pow(geno.maf[i] * (1.0 - geno.maf[i]), *pq_exp);
            for (size_t j = 0; j < annot.cols; ++j) annot(i, j) *= pq;
        }
    }
    std::vector<double> coords(geno.records.size());
    double max_dist = 0.0;
    if (args.ld_wind_snps) {
        max_dist = *args.ld_wind_snps;
        for (size_t i = 0; i < coords.size(); ++i) coords[i] = static_cast<double>(i);
    } else if (args.ld_wind_kb) {
        max_dist = *args.ld_wind_kb * 1000.0;
        for (size_t i = 0; i < coords.size(); ++i) coords[i] = geno.records[i].bp;
    } else {
        max_dist = *args.ld_wind_cm;
        for (size_t i = 0; i < coords.size(); ++i) coords[i] = geno.records[i].cm;
    }
    std::vector<int> left = get_block_lefts(coords, max_dist);
    if (!left.empty() && left.back() == 0 && !args.yes_really) {
        throw Error("Do you really want to compute whole-chromosome LD Score? If so, set --yes-really.");
    }
    log.log("Estimating LD Score.");
    Matrix ld = ldscore_var_blocks_py2(geno, left, args.chunk_size, annot);
    std::vector<std::string> out_names = ldscore_colnames(annot_names, annot.cols, scale_suffix);
    std::set<std::string> print_keep;
    if (args.print_snps) {
        std::vector<std::string> print_list = read_one_col_list(*args.print_snps);
        log.log("Reading list of " + std::to_string(print_list.size()) +
                " SNPs for which to print LD Scores from " + *args.print_snps);
        print_keep.insert(print_list.begin(), print_list.end());
    }
    std::string out_path = args.out + ".l2.ldscore";
    std::ostringstream out;
    out << "CHR\tSNP\tBP";
    for (const auto& n : out_names) out << '\t' << n;
    out << "\n";
    out << std::fixed << std::setprecision(3);
    size_t printed = 0;
    std::vector<size_t> printed_indices;
    for (size_t i = 0; i < geno.records.size(); ++i) {
        if (!print_keep.empty() && !print_keep.count(geno.records[i].snp)) continue;
        ++printed;
        printed_indices.push_back(i);
        out << geno.records[i].chr << '\t' << geno.records[i].snp << '\t' << geno.records[i].bp;
        for (size_t a = 0; a < ld.cols; ++a) out << '\t' << ld(i, a);
        out << "\n";
    }
    if (args.print_snps) {
        if (printed == 0) throw Error("After merging with --print-snps, no SNPs remain.");
        log.log("After merging with --print-snps, LD Scores for " + std::to_string(printed) + " SNPs will be printed.");
    }
    std::vector<double> m(annot.cols, 0.0), m_common(annot.cols, 0.0);
    for (size_t i = 0; i < annot.rows; ++i) {
        for (size_t a = 0; a < annot.cols; ++a) {
            m[a] += annot(i, a);
            if (geno.maf[i] > 0.05) m_common[a] += annot(i, a);
        }
    }
    auto write_m = [&](const std::string& path, const std::vector<double>& vals) {
        std::ofstream f(path);
        if (!f) throw Error("Could not open M output: " + path);
        for (size_t i = 0; i < vals.size(); ++i) {
            if (i) f << '\t';
            f << std::setprecision(12) << vals[i];
        }
        f << "\n";
    };
    write_m(args.out + ".l2.M", m);
    write_m(args.out + ".l2.M_5_50", m_common);
    const std::string gz_path = out_path + ".gz";
    write_gzip_text(gz_path, out.str());
    log.log("Writing LD Scores for " + std::to_string(printed) + " SNPs to " + gz_path);
    if (args.cts_bin && !args.no_print_annot) {
        std::ostringstream annot_out;
        annot_out << "CHR\tSNP\tBP\tCM";
        for (const auto& n : out_names) annot_out << '\t' << n;
        annot_out << "\n";
        for (size_t i = 0; i < geno.records.size(); ++i) {
            annot_out << geno.records[i].chr << '\t' << geno.records[i].snp << '\t' << geno.records[i].bp << '\t'
                      << format_py_float(geno.records[i].cm);
            for (size_t a = 0; a < annot.cols; ++a) annot_out << '\t' << annot(i, a);
            annot_out << "\n";
        }
        log.log("Writing annot matrix produced by --cts-bin to " + gz_path);
        write_gzip_text(args.out + ".annot.gz", annot_out.str());
    }
    std::vector<double> summary_maf;
    Matrix summary_ld(printed_indices.size(), ld.cols);
    summary_maf.reserve(printed_indices.size());
    for (size_t r = 0; r < printed_indices.size(); ++r) {
        const size_t src = printed_indices[r];
        summary_maf.push_back(geno.maf[src]);
        for (size_t c = 0; c < ld.cols; ++c) summary_ld(r, c) = ld(src, c);
    }
    log_ldscore_summary(log, gz_path, summary_maf, summary_ld, out_names);
}

}  // namespace ldsc_cpp
