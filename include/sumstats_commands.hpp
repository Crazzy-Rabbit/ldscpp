#pragma once

#include "regression.hpp"

namespace ldsc_cpp {

void estimate_h2(const Args& args, Logger& log);
void estimate_h2_cts(const Args& args, Logger& log);
RGMergedData merge_other_sumstats(const BaseData& base, const std::string& path, bool no_check_alleles, Logger& log);
void estimate_rg(const Args& args, Logger& log);

}  // namespace ldsc_cpp
