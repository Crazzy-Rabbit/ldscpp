#include "annotation.hpp"

namespace ldsc_cpp {

std::string norm_chr(std::string c) {
    if (c.rfind("chr", 0) == 0) return c;
    return "chr" + c;
}

std::vector<Interval> merge_intervals(std::vector<Interval> intervals) {
    std::sort(intervals.begin(), intervals.end(), [](const Interval& a, const Interval& b) {
        if (a.chr != b.chr) return a.chr < b.chr;
        if (a.start != b.start) return a.start < b.start;
        return a.end < b.end;
    });
    std::vector<Interval> out;
    for (const auto& iv : intervals) {
        if (out.empty() || out.back().chr != iv.chr || iv.start > out.back().end) {
            out.push_back(iv);
        } else {
            out.back().end = std::max(out.back().end, iv.end);
        }
    }
    return out;
}

std::vector<Interval> read_bed_intervals(const std::string& path) {
    std::vector<std::vector<std::string>> rows = read_rows_no_header(path);
    std::vector<Interval> out;
    for (const auto& row : rows) {
        if (row.size() < 3) continue;
        Interval iv;
        iv.chr = norm_chr(row[0]);
        iv.start = static_cast<int>(parse_double(row[1]));
        iv.end = static_cast<int>(parse_double(row[2]));
        out.push_back(iv);
    }
    return out;
}

std::vector<Interval> gene_set_to_intervals(const Args& args) {
    if (!args.gene_set_file) throw Error("--gene-set-file is required.");
    std::vector<std::vector<std::string>> genes_raw = read_rows_no_header(*args.gene_set_file);
    std::set<std::string> genes;
    for (const auto& row : genes_raw) {
        if (!row.empty()) genes.insert(row[0]);
    }
    Table coords = read_table(args.gene_coord_file);
    int gene_i = coords.col("GENE");
    int chr_i = coords.col("CHR");
    int start_i = coords.col("START");
    int end_i = coords.col("END");
    if (gene_i < 0 || chr_i < 0 || start_i < 0 || end_i < 0) {
        throw Error("--gene-coord-file must contain GENE, CHR, START, END columns.");
    }
    std::vector<Interval> intervals;
    for (const auto& row : coords.rows) {
        if (!genes.count(row[static_cast<size_t>(gene_i)])) continue;
        Interval iv;
        iv.chr = norm_chr(row[static_cast<size_t>(chr_i)]);
        iv.start = std::max(1, static_cast<int>(parse_double(row[static_cast<size_t>(start_i)])) - args.windowsize) - 1;
        iv.end = static_cast<int>(parse_double(row[static_cast<size_t>(end_i)])) + args.windowsize;
        intervals.push_back(iv);
    }
    return merge_intervals(intervals);
}

void make_annot_cpp(const Args& args) {
    if (!args.bimfile || !args.annot_file) throw Error("--bimfile and --annot-file are required.");
    if (!args.gene_set_file && !args.bed_file) throw Error("Must specify either --gene-set-file or --bed-file.");
    std::vector<Interval> intervals = args.gene_set_file ? gene_set_to_intervals(args) : read_bed_intervals(args.bed_file.value_or(""));
    if (!args.nomerge) intervals = merge_intervals(intervals);
    std::vector<std::vector<std::string>> bim = read_rows_no_header(*args.bimfile);
    std::ofstream out(*args.annot_file);
    if (!out) throw Error("Could not open annot output: " + *args.annot_file);
    out << "ANNOT\n";
    for (const auto& row : bim) {
        if (row.size() < 4) throw Error("BIM file must have at least 4 columns.");
        std::string chr = norm_chr(row[0]);
        int bp = static_cast<int>(parse_double(row[3]));
        int count = 0;
        for (const auto& iv : intervals) {
            if (iv.chr == chr && bp - 1 >= iv.start && bp - 1 < iv.end) ++count;
        }
        out << (args.nomerge ? count : (count > 0 ? 1 : 0)) << "\n";
    }
}

}  // namespace ldsc_cpp
