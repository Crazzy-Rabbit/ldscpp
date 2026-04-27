#include "munge.hpp"

namespace ldsc_cpp {

std::string clean_header_cpp(std::string s) {
    for (char& c : s) {
        if (c == '-' || c == '.') c = '_';
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return s;
}

std::map<std::string, std::string> default_cname_map() {
    return {
        {"SNP", "SNP"}, {"MARKERNAME", "SNP"}, {"SNPID", "SNP"}, {"RS", "SNP"}, {"RSID", "SNP"},
        {"RS_NUMBER", "SNP"}, {"RS_NUMBERS", "SNP"},
        {"NSTUDY", "NSTUDY"}, {"N_STUDY", "NSTUDY"}, {"NSTUDIES", "NSTUDY"}, {"N_STUDIES", "NSTUDY"},
        {"P", "P"}, {"PVALUE", "P"}, {"P_VALUE", "P"}, {"PVAL", "P"}, {"P_VAL", "P"}, {"GC_PVALUE", "P"},
        {"A1", "A1"}, {"ALLELE1", "A1"}, {"ALLELE_1", "A1"}, {"EFFECT_ALLELE", "A1"},
        {"REFERENCE_ALLELE", "A1"}, {"INC_ALLELE", "A1"}, {"EA", "A1"},
        {"A2", "A2"}, {"ALLELE2", "A2"}, {"ALLELE_2", "A2"}, {"OTHER_ALLELE", "A2"},
        {"NON_EFFECT_ALLELE", "A2"}, {"DEC_ALLELE", "A2"}, {"NEA", "A2"},
        {"N", "N"}, {"WEIGHT", "N"},
        {"NCASE", "N_CAS"}, {"CASES_N", "N_CAS"}, {"N_CASE", "N_CAS"}, {"N_CASES", "N_CAS"}, {"N_CAS", "N_CAS"},
        {"NCONTROL", "N_CON"}, {"CONTROLS_N", "N_CON"}, {"N_CONTROL", "N_CON"}, {"N_CONTROLS", "N_CON"},
        {"N_CON", "N_CON"},
        {"ZSCORE", "Z"}, {"Z_SCORE", "Z"}, {"GC_ZSCORE", "Z"}, {"Z", "Z"},
        {"OR", "OR"}, {"B", "BETA"}, {"BETA", "BETA"}, {"LOG_ODDS", "LOG_ODDS"},
        {"EFFECTS", "BETA"}, {"EFFECT", "BETA"}, {"SIGNED_SUMSTAT", "SIGNED_SUMSTAT"},
        {"INFO", "INFO"},
        {"EAF", "FRQ"}, {"FRQ", "FRQ"}, {"MAF", "FRQ"}, {"FRQ_U", "FRQ"}, {"F_U", "FRQ"}
    };
}

std::map<std::string, double> null_values_map() {
    return {{"LOG_ODDS", 0.0}, {"BETA", 0.0}, {"OR", 1.0}, {"Z", 0.0}};
}

std::map<std::string, std::string> cname_descriptions() {
    return {
        {"SNP", "Variant ID (e.g., rs number)"},
        {"P", "p-Value"},
        {"A1", "Allele 1, interpreted as ref allele for signed sumstat."},
        {"A2", "Allele 2, interpreted as non-ref allele for signed sumstat."},
        {"N", "Sample size"},
        {"N_CAS", "Number of cases"},
        {"N_CON", "Number of controls"},
        {"Z", "Z-score (0 --> no effect; above 0 --> A1 is trait/risk increasing)"},
        {"OR", "Odds ratio (1 --> no effect; above 1 --> A1 is risk increasing)"},
        {"BETA", "[linear/logistic] regression coefficient (0 --> no effect; above 0 --> A1 is trait/risk increasing)"},
        {"LOG_ODDS", "Log odds ratio (0 --> no effect; above 0 --> A1 is risk increasing)"},
        {"INFO", "INFO score (imputation quality; higher --> better imputation)"},
        {"FRQ", "Allele frequency"},
        {"SIGNED_SUMSTAT", "Directional summary statistic as specified by --signed-sumstats."},
        {"NSTUDY", "Number of studies in which the SNP was genotyped."}
    };
}

double quantile(std::vector<double> x, double q) {
    if (x.empty()) throw Error("Cannot take quantile of empty vector.");
    std::sort(x.begin(), x.end());
    const double pos = (static_cast<double>(x.size()) - 1.0) * q;
    const size_t lo = static_cast<size_t>(std::floor(pos));
    const size_t hi = static_cast<size_t>(std::ceil(pos));
    if (lo == hi) return x[lo];
    return x[lo] + (x[hi] - x[lo]) * (pos - static_cast<double>(lo));
}

double p_to_z(double p) {
    if (!(p > 0.0 && p <= 1.0)) return kNan;
    return norm_ppf(1.0 - p / 2.0);
}

std::string py2_str_float(double x) {
    std::ostringstream out;
    out << std::setprecision(12) << x;
    std::string s = out.str();
    if (s.find('.') == std::string::npos && s.find('e') == std::string::npos && s.find('E') == std::string::npos) {
        s += ".0";
    }
    return s;
}

std::string fixed3(double x) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(3) << x;
    return out.str();
}

std::map<std::string, std::string> flag_cname_map(const Args& args, double& signed_null) {
    std::map<std::string, std::string> out;
    auto add = [&](const std::optional<std::string>& v, const std::string& internal) {
        if (v) out[clean_header_cpp(*v)] = internal;
    };
    add(args.nstudy_col, "NSTUDY");
    add(args.snp_col, "SNP");
    add(args.n_col, "N");
    add(args.n_cas_col, "N_CAS");
    add(args.n_con_col, "N_CON");
    add(args.a1_col, "A1");
    add(args.a2_col, "A2");
    add(args.p_col, "P");
    add(args.frq_col, "FRQ");
    add(args.info_col, "INFO");
    if (args.info_list) {
        for (const auto& x : split_char(*args.info_list, ',')) out[clean_header_cpp(x)] = "INFO";
    }
    if (args.signed_sumstats) {
        std::vector<std::string> parts = split_char(*args.signed_sumstats, ',');
        if (parts.size() != 2) throw Error("The argument to --signed-sumstats should be column header comma number.");
        signed_null = parse_double(parts[1]);
        out[clean_header_cpp(parts[0])] = "SIGNED_SUMSTAT";
    }
    return out;
}

std::map<std::string, std::string> build_cname_translation(const Table& t, const Args& args, double& signed_null) {
    std::set<std::string> ignore;
    if (args.ignore_cols) {
        for (const auto& x : split_char(*args.ignore_cols, ',')) ignore.insert(clean_header_cpp(x));
    }
    std::map<std::string, std::string> flag = flag_cname_map(args, signed_null);
    std::map<std::string, std::string> defs = default_cname_map();
    if (args.signed_sumstats || args.a1_inc) {
        std::map<std::string, double> nulls = null_values_map();
        for (auto it = defs.begin(); it != defs.end();) {
            if (nulls.count(it->second)) it = defs.erase(it);
            else ++it;
        }
    }
    std::map<std::string, std::string> clean_to_internal;
    for (const auto& kv : flag) {
        if (!ignore.count(kv.first)) clean_to_internal[kv.first] = kv.second;
    }
    for (const auto& kv : defs) {
        if (!ignore.count(kv.first) && !clean_to_internal.count(kv.first)) clean_to_internal[kv.first] = kv.second;
    }

    std::map<std::string, std::string> translation;
    for (const auto& c : t.columns) {
        std::string clean = clean_header_cpp(c);
        auto it = clean_to_internal.find(clean);
        if (it != clean_to_internal.end()) translation[c] = it->second;
    }
    std::map<std::string, double> nulls = null_values_map();
    if (!args.signed_sumstats && !args.a1_inc) {
        std::vector<std::string> sign_cols;
        for (const auto& kv : translation) {
            if (nulls.count(kv.second)) sign_cols.push_back(kv.first);
        }
        if (sign_cols.size() > 1) throw Error("Too many signed sumstat columns. Specify which to ignore with --ignore.");
        if (sign_cols.empty()) throw Error("Could not find a signed summary statistic column.");
        signed_null = nulls[translation[sign_cols[0]]];
        translation[sign_cols[0]] = "SIGNED_SUMSTAT";
    }
    return translation;
}

void require_internal(const std::map<std::string, std::string>& tr, const std::string& name) {
    for (const auto& kv : tr) {
        if (kv.second == name) return;
    }
    throw Error("Could not find " + name + " column.");
}

void munge_sumstats_cpp(const Args& args) {
    const auto start = std::chrono::steady_clock::now();
    if (!args.sumstats) throw Error("The --sumstats flag is required.");
    if (!args.out_set) throw Error("The --out flag is required.");
    if (args.no_alleles && args.merge_alleles) throw Error("--no-alleles and --merge-alleles are not compatible.");
    if (args.daner && args.daner_n) throw Error("--daner and --daner-n are not compatible.");
    Logger log(args.out + ".log");
    log.log(command_header_from_args(args.raw_args, args.program_name));
    Table t = read_table(*args.sumstats);
    std::optional<double> forced_n = args.n_value;
    std::optional<double> forced_n_cas = args.n_cas;
    std::optional<double> forced_n_con = args.n_con;
    double signed_null = kNan;
    std::map<std::string, std::string> tr = build_cname_translation(t, args, signed_null);
    if (args.daner || args.daner_n) {
        std::optional<std::string> frq_u;
        std::optional<std::string> frq_a;
        for (const auto& c : t.columns) {
            if (c.rfind("FRQ_U_", 0) == 0) frq_u = c;
            if (c.rfind("FRQ_A_", 0) == 0) frq_a = c;
        }
        if (!frq_u) throw Error("Could not find FRQ_U_* column expected for daner format.");
        tr[*frq_u] = "FRQ";
        if (args.daner) {
            if (!frq_a) throw Error("Could not find FRQ_A_* column expected for daner format.");
            forced_n_con = parse_double(frq_u->substr(6));
            forced_n_cas = parse_double(frq_a->substr(6));
            log.log("Inferred that N_cas = " + format_py_float(*forced_n_cas) + ", N_con = " +
                    format_py_float(*forced_n_con) + " from the FRQ_[A/U] columns.");
        } else {
            for (const auto& c : t.columns) {
                if (c == "Nca") tr[c] = "N_CAS";
                if (c == "Nco") tr[c] = "N_CON";
            }
        }
    }
    if (args.a1_inc) {
        require_internal(tr, "SNP");
        require_internal(tr, "P");
    } else {
        require_internal(tr, "SNP");
        require_internal(tr, "P");
        require_internal(tr, "SIGNED_SUMSTAT");
    }
    if (!args.no_alleles) {
        require_internal(tr, "A1");
        require_internal(tr, "A2");
    }

    std::map<std::string, std::string> by_internal;
    for (const auto& kv : tr) {
        if (by_internal.count(kv.second) && kv.second != "INFO") {
            throw Error("Found multiple columns mapping to " + kv.second);
        }
        by_internal[kv.second] = kv.first;
    }

    if (!forced_n && !(forced_n_cas && forced_n_con) && !by_internal.count("N") &&
        !(by_internal.count("N_CAS") && by_internal.count("N_CON"))) {
        throw Error("Could not determine N.");
    }
    std::string sign_cname = "signed sumstat";
    for (const auto& kv : tr) {
        if (kv.second == "SIGNED_SUMSTAT") {
            sign_cname = kv.first;
            break;
        }
    }
    if (args.signed_sumstats) sign_cname = "SIGNED_SUMSTATS";
    std::map<std::string, std::string> descriptions = cname_descriptions();
    struct DescEntry {
        std::string column;
        std::string internal;
    };
    std::vector<DescEntry> desc_entries;
    for (const auto& kv : tr) {
        std::string internal = kv.second;
        if (internal == "SIGNED_SUMSTAT" && !args.signed_sumstats) {
            std::map<std::string, std::string> defs = default_cname_map();
            auto original = defs.find(clean_header_cpp(kv.first));
            if (original != defs.end() && null_values_map().count(original->second)) internal = original->second;
        }
        desc_entries.push_back({kv.first, internal});
    }
    bool desc_has_beta = false;
    bool desc_has_n = false;
    for (const auto& entry : desc_entries) {
        desc_has_beta = desc_has_beta || entry.internal == "BETA";
        desc_has_n = desc_has_n || entry.internal == "N";
    }
    std::vector<std::string> py2_desc_order = {
        "INFO", "A1", "P", "FRQ", "A2", "SNP", "OR", "BETA", "LOG_ODDS", "Z",
        "SIGNED_SUMSTAT", "N", "N_CAS", "N_CON", "NSTUDY"
    };
    if (args.daner) {
        py2_desc_order = {"INFO", "A1", "P", "A2", "SNP", "FRQ", "OR", "BETA", "LOG_ODDS", "Z",
                          "SIGNED_SUMSTAT", "N", "N_CAS", "N_CON", "NSTUDY"};
    } else if (args.signed_sumstats) {
        py2_desc_order = {"P", "SIGNED_SUMSTAT", "SNP", "N", "INFO", "A1", "FRQ", "A2", "OR", "BETA",
                          "LOG_ODDS", "Z", "N_CAS", "N_CON", "NSTUDY"};
    } else if (desc_has_beta && desc_has_n) {
        py2_desc_order = {"INFO", "BETA", "N", "A1", "P", "FRQ", "A2", "SNP", "OR", "LOG_ODDS", "Z",
                          "SIGNED_SUMSTAT", "N_CAS", "N_CON", "NSTUDY"};
    }
    log.log("Interpreting column names as follows:");
    std::ostringstream desc;
    bool first_desc = true;
    for (const auto& key : py2_desc_order) {
        for (const auto& entry : desc_entries) {
            if (entry.internal != key) continue;
            auto di = descriptions.find(entry.internal);
            if (di == descriptions.end()) continue;
            if (!first_desc) desc << "\n";
            desc << entry.column << ":\t" << di->second;
            first_desc = false;
        }
    }
    desc << "\n";
    log.log(desc.str());

    std::unordered_set<std::string> merge_keep;
    std::unordered_map<std::string, std::string> merge_ma;
    std::vector<std::string> merge_order;
    if (args.merge_alleles) {
        log.log("Reading list of SNPs for allele merge from " + *args.merge_alleles);
        Table mt = read_table(*args.merge_alleles);
        int snp_i = mt.col("SNP");
        int a1_i = mt.col("A1");
        int a2_i = mt.col("A2");
        if (snp_i < 0 || a1_i < 0 || a2_i < 0) throw Error("--merge-alleles must have columns SNP, A1, A2.");
        for (const auto& row : mt.rows) {
            std::string snp = row[static_cast<size_t>(snp_i)];
            merge_keep.insert(snp);
            merge_order.push_back(snp);
            merge_ma[snp] = row[static_cast<size_t>(a1_i)] + row[static_cast<size_t>(a2_i)];
            for (char& c : merge_ma[snp]) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        }
        log.log("Read " + std::to_string(merge_keep.size()) + " SNPs for allele merge.");
    }
    log.log("Reading sumstats from " + *args.sumstats + " into memory " + std::to_string(args.chunksize) +
            " SNPs at a time.");

    std::vector<MungeRow> rows;
    size_t drop_na = 0, drop_p = 0, drop_info = 0, drop_frq = 0, drop_a = 0, drop_merge = 0;
    std::vector<size_t> info_cols;
    for (const auto& kv : tr) {
        if (kv.second == "INFO") info_cols.push_back(t.index.at(kv.first));
    }

    for (const auto& row : t.rows) {
        auto get = [&](const std::string& internal) -> std::optional<std::string> {
            auto it = by_internal.find(internal);
            if (it == by_internal.end()) return std::nullopt;
            return row[t.index.at(it->second)];
        };
        MungeRow r;
        r.snp = get("SNP").value_or("");
        if (r.snp.empty() || r.snp == "." || r.snp == "NA") {
            ++drop_na;
            continue;
        }
        if (args.merge_alleles && !merge_keep.count(r.snp)) {
            ++drop_merge;
            continue;
        }
        try {
            r.p = parse_double(get("P").value_or("."));
            if (!args.a1_inc) r.signed_sumstat = parse_double(get("SIGNED_SUMSTAT").value_or("."));
            if (auto v = get("N")) r.n = parse_double(*v);
            if (auto v = get("N_CAS")) r.n_cas = parse_double(*v);
            if (auto v = get("N_CON")) r.n_con = parse_double(*v);
            if (auto v = get("NSTUDY")) r.nstudy = parse_double(*v);
            if (auto v = get("FRQ")) r.frq = parse_double(*v);
        } catch (const Error&) {
            ++drop_na;
            continue;
        }
        if (std::isnan(r.p) || (!args.a1_inc && std::isnan(r.signed_sumstat))) {
            ++drop_na;
            continue;
        }
        if (!args.no_alleles) {
            r.a1 = get("A1").value_or("");
            r.a2 = get("A2").value_or("");
            for (char& c : r.a1) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            for (char& c : r.a2) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            if (r.a1.empty() || r.a2.empty() || r.a1 == "." || r.a2 == ".") {
                ++drop_na;
                continue;
            }
        }
        bool ok = true;
        if (!info_cols.empty()) {
            double info_sum = 0.0;
            int info_n = 0;
            for (size_t c : info_cols) {
                double info = parse_double(row[c]);
                if (!std::isnan(info)) {
                    info_sum += info;
                    ++info_n;
                }
            }
            if (info_n > 0 && info_sum < args.info_min * static_cast<double>(info_n)) ok = false;
            if (!ok) ++drop_info;
        }
        if (ok && !std::isnan(r.frq)) {
            double maf = std::min(r.frq, 1.0 - r.frq);
            if (r.frq < 0.0 || r.frq > 1.0 || maf <= args.maf_min) {
                ok = false;
                ++drop_frq;
            }
        }
        if (ok && !(r.p > 0.0 && r.p <= 1.0)) {
            ok = false;
            ++drop_p;
        }
        if (ok && !args.no_alleles && !alleles_match(r.a1 + r.a2 + r.a1 + r.a2)) {
            ok = false;
            ++drop_a;
        }
        if (!ok) continue;
        rows.push_back(std::move(r));
    }

    log.log("Read " + std::to_string(t.rows.size()) + " SNPs from --sumstats file.");
    if (args.merge_alleles) log.log("Removed " + std::to_string(drop_merge) + " SNPs not in --merge-alleles.");
    log.log("Removed " + std::to_string(drop_na) + " SNPs with missing values.");
    log.log("Removed " + std::to_string(drop_info) + " SNPs with INFO <= " + format_num(args.info_min) + ".");
    log.log("Removed " + std::to_string(drop_frq) + " SNPs with MAF <= " + format_num(args.maf_min) + ".");
    log.log("Removed " + std::to_string(drop_p) + " SNPs with out-of-bounds p-values.");
    log.log("Removed " + std::to_string(drop_a) + " variants that were not SNPs or were strand-ambiguous.");
    log.log(std::to_string(rows.size()) + " SNPs remain.");

    std::unordered_set<std::string> seen;
    std::vector<MungeRow> dedup;
    for (auto& r : rows) {
        if (seen.insert(r.snp).second) dedup.push_back(std::move(r));
    }
    log.log("Removed " + std::to_string(rows.size() - dedup.size()) + " SNPs with duplicated rs numbers (" +
            std::to_string(dedup.size()) + " SNPs remain).");
    rows = std::move(dedup);

    const bool n_from_columns = by_internal.count("N") || (by_internal.count("N_CAS") && by_internal.count("N_CON"));
    if (!std::isnan(rows.empty() ? kNan : rows[0].n_cas)) {
        for (auto& r : rows) {
            if (!std::isnan(r.n_cas) && !std::isnan(r.n_con)) r.n = r.n_cas + r.n_con;
        }
    }
    bool has_n = n_from_columns;
    if (n_from_columns) {
        for (const auto& r : rows) has_n = has_n && !std::isnan(r.n);
    }
    if (has_n) {
        std::vector<double> ns;
        for (const auto& r : rows) ns.push_back(r.n);
        double nmin = args.n_min.value_or(quantile(ns, 0.9) / 1.5);
        std::vector<MungeRow> kept;
        for (auto& r : rows) {
            if (r.n >= nmin) kept.push_back(std::move(r));
        }
        log.log("Removed " + std::to_string(rows.size() - kept.size()) + " SNPs with N < " + py2_str_float(nmin) +
                " (" + std::to_string(kept.size()) + " SNPs remain).");
        rows = std::move(kept);
    } else if (forced_n) {
        for (auto& r : rows) r.n = *forced_n;
        log.log("Using N = " + format_py_float(*forced_n));
    } else if (forced_n_cas && forced_n_con) {
        for (auto& r : rows) r.n = *forced_n_cas + *forced_n_con;
        if (!args.daner) {
            log.log("Using N_cas = " + format_py_float(*forced_n_cas) + "; N_con = " +
                    format_py_float(*forced_n_con));
        }
    } else {
        throw Error("Cannot determine N.");
    }

    std::vector<double> signed_values;
    for (auto& r : rows) {
        r.z = p_to_z(r.p);
        if (!args.a1_inc) {
            signed_values.push_back(r.signed_sumstat);
            if (r.signed_sumstat < signed_null) r.z *= -1.0;
        }
    }
    if (!args.a1_inc) {
        double med = median(signed_values);
        if (std::abs(med - signed_null) > 0.1) {
            throw Error("WARNING: median value of " + sign_cname + " is " + format_num(med, 2) +
                        " (should be close to " + format_py_float(signed_null) + "). This column may be mislabeled.");
        }
        log.log("Median value of " + sign_cname + " was " + format_py_float(med) + ", which seems sensible.");
    }

    if (args.merge_alleles) {
        std::unordered_map<std::string, size_t> row_by_snp;
        for (size_t i = 0; i < rows.size(); ++i) row_by_snp.emplace(rows[i].snp, i);
        std::vector<MungeRow> merged;
        merged.reserve(merge_order.size());
        size_t old = 0;
        size_t mismatch = 0;
        for (const auto& snp : merge_order) {
            MungeRow out_row;
            out_row.snp = snp;
            auto ri = row_by_snp.find(snp);
            if (ri != row_by_snp.end()) {
                ++old;
                const MungeRow& r = rows[ri->second];
                auto ai = merge_ma.find(snp);
                if (ai != merge_ma.end() && alleles_match(r.a1 + r.a2 + ai->second)) {
                    out_row = r;
                } else {
                    ++mismatch;
                }
            }
            merged.push_back(std::move(out_row));
        }
        if (mismatch < old) {
            log.log("Removed " + std::to_string(mismatch) +
                    " SNPs whose alleles did not match --merge-alleles (" +
                    std::to_string(old - mismatch) + " SNPs remain).");
        } else {
            throw Error("All SNPs have alleles that do not match --merge-alleles.");
        }
        rows = std::move(merged);
    }

    std::string out_path = args.out + ".sumstats";
    std::ostringstream out;
    out << "SNP";
    if (!args.no_alleles) out << "\tA1\tA2";
    out << "\tZ\tN";
    if (args.keep_maf && by_internal.count("FRQ")) out << "\tFRQ";
    out << "\n";
    size_t nonmissing = 0;
    for (const auto& r : rows) {
        out << r.snp;
        if (!args.no_alleles) out << '\t' << r.a1 << '\t' << r.a2;
        out << '\t';
        if (!std::isnan(r.z)) out << fixed3(r.z);
        out << '\t';
        if (!std::isnan(r.n)) {
            if (n_from_columns && std::abs(r.n - std::round(r.n)) < 1e-9) out << static_cast<long long>(std::llround(r.n));
            else out << fixed3(r.n);
            ++nonmissing;
        }
        if (args.keep_maf && by_internal.count("FRQ")) out << '\t' << fixed3(r.frq);
        out << "\n";
    }
    const std::string gz_path = out_path + ".gz";
    log.log("Writing summary statistics for " + std::to_string(rows.size()) + " SNPs (" + std::to_string(nonmissing) +
            " with nonmissing beta) to " + gz_path + ".");
    write_gzip_text(gz_path, out.str());
    std::vector<double> chisq;
    for (const auto& r : rows) {
        if (!std::isnan(r.z)) chisq.push_back(r.z * r.z);
    }
    if (!chisq.empty()) {
        log.log("\nMetadata:");
        log.log("Mean chi^2 = " + format_num(mean(chisq), 3));
        if (mean(chisq) < 1.02) log.log("WARNING: mean chi^2 may be too small.");
        log.log("Lambda GC = " + format_num(median(chisq) / 0.4549, 3));
        log.log("Max chi^2 = " + format_num(*std::max_element(chisq.begin(), chisq.end()), 3));
        size_t gws = std::count_if(chisq.begin(), chisq.end(), [](double x) { return x > 29.0; });
        log.log(std::to_string(gws) + " Genome-wide significant SNPs (some may have been removed by filtering).");
    }
    const auto end = std::chrono::steady_clock::now();
    const double elapsed = std::round(std::chrono::duration<double>(end - start).count() * 100.0) / 100.0;
    log.log("\nConversion finished at " + time_now_string());
    log.log("Total time elapsed: " + sec_to_str(elapsed));
}

}  // namespace ldsc_cpp
