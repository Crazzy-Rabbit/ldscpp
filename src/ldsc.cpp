#include "ldsc.hpp"

namespace ldsc_cpp {

void print_help() {
    std::cout << "ldsc: C++ LDSC summary-statistics regression\n"
              << "Subcommands: h2, rg, munge, make-annot, ld-score, quantile-m, quantile-h2g\n"
              << "Required for h2: --h2 FILE --ref-ld PREFIX --w-ld PREFIX\n"
              << "Required for rg: --rg FILE1,FILE2[,FILE3...] --ref-ld PREFIX --w-ld PREFIX\n"
              << "Required for munge: munge --sumstats FILE --out PREFIX [--N N]\n"
              << "Required for make-annot: make-annot --bed-file BED --bimfile BIM --annot-file OUT\n"
              << "Required for ld-score: ld-score --bfile PREFIX --l2 --ld-wind-snps N --out PREFIX\n";
}

Args parse_args(int argc, char** argv) {
    Args args;
    if (argc > 0) args.program_name = argv[0];
    for (int i = 1; i < argc; ++i) args.raw_args.emplace_back(argv[i]);
    int start = 1;
    if (argc > 1 && std::string(argv[1]).rfind("--", 0) != 0) {
        args.command = argv[1];
        start = 2;
    }
    for (int i = start; i < argc; ++i) {
        std::string key = argv[i];
        auto need_value = [&](const std::string& k) -> std::string {
            if (i + 1 >= argc) throw Error("Missing value for " + k);
            return argv[++i];
        };
        if (key == "-h" || key == "--help") {
            print_help();
            std::exit(0);
        } else if (key == "--h2") {
            args.h2 = need_value(key);
        } else if (key == "--h2-cts") {
            args.h2_cts = need_value(key);
        } else if (key == "--rg") {
            args.rg = need_value(key);
        } else if (key == "--ref-ld") {
            args.ref_ld = need_value(key);
        } else if (key == "--ref-ld-chr") {
            args.ref_ld_chr = need_value(key);
        } else if (key == "--w-ld") {
            args.w_ld = need_value(key);
        } else if (key == "--w-ld-chr") {
            args.w_ld_chr = need_value(key);
        } else if (key == "--M") {
            args.manual_m = need_value(key);
        } else if (key == "--intercept-h2") {
            std::string v = need_value(key);
            args.intercept_h2_list = v;
            if (v.find(',') == std::string::npos) args.intercept_h2_scalar = parse_double(v);
        } else if (key == "--intercept-gencov") {
            args.intercept_gencov_list = need_value(key);
        } else if (key == "--two-step") {
            args.twostep = parse_double(need_value(key));
        } else if (key == "--chisq-max") {
            args.chisq_max = parse_double(need_value(key));
        } else if (key == "--n-blocks") {
            args.n_blocks = std::stoi(need_value(key));
        } else if (key == "--out") {
            args.out = need_value(key);
            args.out_set = true;
        } else if (key == "--samp-prev") {
            args.samp_prev = need_value(key);
        } else if (key == "--pop-prev") {
            args.pop_prev = need_value(key);
        } else if (key == "--no-intercept") {
            args.no_intercept = true;
        } else if (key == "--not-M-5-50") {
            args.not_m_5_50 = true;
        } else if (key == "--print-cov") {
            args.print_cov = true;
        } else if (key == "--print-delete-vals") {
            args.print_delete_vals = true;
        } else if (key == "--no-check-alleles") {
            args.no_check_alleles = true;
        } else if (key == "--return-silly-things") {
            args.return_silly_things = true;
        } else if (key == "--ref-ld-chr-cts") {
            args.ref_ld_chr_cts = need_value(key);
        } else if (key == "--print-all-cts") {
            args.print_all_cts = true;
        } else if (key == "--overlap-annot") {
            args.overlap_annot = true;
        } else if (key == "--print-coefficients") {
            args.print_coefficients = true;
        } else if (key == "--frqfile") {
            args.frqfile = need_value(key);
        } else if (key == "--frqfile-chr") {
            args.frqfile_chr = need_value(key);
        } else if (key == "--sumstats") {
            args.sumstats = need_value(key);
        } else if (key == "--N") {
            args.n_value = parse_double(need_value(key));
        } else if (key == "--N-cas") {
            args.n_cas = parse_double(need_value(key));
        } else if (key == "--N-con") {
            args.n_con = parse_double(need_value(key));
        } else if (key == "--info-min") {
            args.info_min = parse_double(need_value(key));
        } else if (key == "--maf-min") {
            args.maf_min = parse_double(need_value(key));
        } else if (key == "--daner") {
            args.daner = true;
        } else if (key == "--daner-n") {
            args.daner_n = true;
        } else if (key == "--no-alleles") {
            args.no_alleles = true;
        } else if (key == "--merge-alleles") {
            args.merge_alleles = need_value(key);
        } else if (key == "--n-min") {
            args.n_min = parse_double(need_value(key));
        } else if (key == "--snp") {
            args.snp_col = need_value(key);
        } else if (key == "--N-col") {
            args.n_col = need_value(key);
        } else if (key == "--N-cas-col") {
            args.n_cas_col = need_value(key);
        } else if (key == "--N-con-col") {
            args.n_con_col = need_value(key);
        } else if (key == "--a1") {
            args.a1_col = need_value(key);
        } else if (key == "--a2") {
            args.a2_col = need_value(key);
        } else if (key == "--p") {
            args.p_col = need_value(key);
        } else if (key == "--frq") {
            args.frq_col = need_value(key);
        } else if (key == "--signed-sumstats") {
            args.signed_sumstats = need_value(key);
        } else if (key == "--info") {
            args.info_col = need_value(key);
        } else if (key == "--info-list") {
            args.info_list = need_value(key);
        } else if (key == "--nstudy") {
            args.nstudy_col = need_value(key);
        } else if (key == "--nstudy-min") {
            args.nstudy_min = parse_double(need_value(key));
        } else if (key == "--ignore") {
            args.ignore_cols = need_value(key);
        } else if (key == "--a1-inc") {
            args.a1_inc = true;
        } else if (key == "--keep-maf") {
            args.keep_maf = true;
        } else if (key == "--gene-set-file") {
            args.gene_set_file = need_value(key);
        } else if (key == "--gene-coord-file") {
            args.gene_coord_file = need_value(key);
        } else if (key == "--windowsize") {
            args.windowsize = std::stoi(need_value(key));
        } else if (key == "--bed-file") {
            args.bed_file = need_value(key);
        } else if (key == "--nomerge") {
            args.nomerge = true;
        } else if (key == "--bimfile") {
            args.bimfile = need_value(key);
        } else if (key == "--annot-file") {
            args.annot_file = need_value(key);
        } else if (key == "--bfile") {
            args.bfile = need_value(key);
        } else if (key == "--l2") {
            args.l2 = true;
        } else if (key == "--extract") {
            args.extract = need_value(key);
        } else if (key == "--keep") {
            args.keep = need_value(key);
        } else if (key == "--ld-wind-snps") {
            args.ld_wind_snps = std::stoi(need_value(key));
        } else if (key == "--ld-wind-kb") {
            args.ld_wind_kb = parse_double(need_value(key));
        } else if (key == "--ld-wind-cm") {
            args.ld_wind_cm = parse_double(need_value(key));
        } else if (key == "--print-snps") {
            args.print_snps = need_value(key);
        } else if (key == "--annot") {
            args.annot = need_value(key);
        } else if (key == "--thin-annot") {
            args.thin_annot = true;
        } else if (key == "--cts-bin") {
            args.cts_bin = need_value(key);
        } else if (key == "--cts-breaks") {
            args.cts_breaks = need_value(key);
        } else if (key == "--cts-names") {
            args.cts_names = need_value(key);
        } else if (key == "--per-allele") {
            args.per_allele = true;
        } else if (key == "--pq-exp") {
            args.pq_exp = parse_double(need_value(key));
        } else if (key == "--no-print-annot") {
            args.no_print_annot = true;
        } else if (key == "--maf") {
            args.maf = parse_double(need_value(key));
            args.maf_set = true;
        } else if (key == "--chunk-size") {
            args.chunk_size = std::stoi(need_value(key));
        } else if (key == "--chunksize") {
            args.chunksize = std::stoi(need_value(key));
        } else if (key == "--pickle") {
            args.pickle = true;
        } else if (key == "--invert-anyway") {
            args.invert_anyway = true;
        } else if (key == "--yes-really") {
            args.yes_really = true;
        } else if (key == "--annotfile") {
            args.quantile_annotfile = need_value(key);
        } else if (key == "--resultfile") {
            args.quantile_resultfile = need_value(key);
        } else if (key == "--outfile") {
            args.quantile_outfile = need_value(key);
        } else if (key == "--ref-annot-chr") {
            args.ref_annot_chr = need_value(key);
        } else if (key == "--annot-header") {
            args.annot_header = need_value(key);
        } else if (key == "--nb-quantile") {
            args.nb_quantile = std::stoi(need_value(key));
        } else if (key == "--exclude0") {
            args.exclude0 = true;
        } else if (key.rfind("--", 0) != 0) {
            args.positionals.push_back(key);
        } else {
            throw Error("Unknown argument: " + key);
        }
    }
    return args;
}

}  // namespace ldsc_cpp

using namespace ldsc_cpp;

int main(int argc, char** argv) {
    const auto start = std::chrono::steady_clock::now();
    try {
        ldsc_cpp::Args args = ldsc_cpp::parse_args(argc, argv);
        if (args.n_blocks <= 1) throw Error("--n-blocks must be an integer > 1.");
        if ((args.h2 ? 1 : 0) + (args.h2_cts ? 1 : 0) + (args.rg ? 1 : 0) > 1) {
            throw Error("Cannot set more than one of --h2, --h2-cts, and --rg.");
        }
        if (args.ref_ld && args.ref_ld_chr) throw Error("Cannot set both --ref-ld and --ref-ld-chr.");
        if (args.w_ld && args.w_ld_chr) throw Error("Cannot set both --w-ld and --w-ld-chr.");
        if (args.samp_prev.has_value() != args.pop_prev.has_value()) throw Error("Must set both or neither of --samp-prev and --pop-prev.");
        ensure_parent_dir(args.out);
        if (args.command == "munge" || args.sumstats) {
            munge_sumstats_cpp(args);
            return 0;
        }
        if (args.command == "make-annot" || args.annot_file) {
            make_annot_cpp(args);
            return 0;
        }
        if (args.command == "quantile-h2g") {
            quantile_h2g_cpp(args);
            return 0;
        }
        if (args.command == "quantile-m") {
            quantile_m_cpp(args);
            return 0;
        }
        Logger log(args.out + ".log");
        log.log(command_header_from_args(args.raw_args, args.program_name));
        log.log("Beginning analysis at " + time_now_string());
        if (args.command == "ld-score" || (args.bfile && args.l2)) {
            estimate_ldscore_cpp(args, log);
        } else if (args.h2_cts) {
            estimate_h2_cts(args, log);
        } else if (args.h2) {
            estimate_h2(args, log);
        } else if (args.rg) {
            estimate_rg(args, log);
        } else {
            print_help();
            return 2;
        }
        const auto end = std::chrono::steady_clock::now();
        const double elapsed = std::round(std::chrono::duration<double>(end - start).count() * 100.0) / 100.0;
        log.log("Analysis finished at " + time_now_string());
        log.log("Total time elapsed: " + sec_to_str(elapsed));
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
