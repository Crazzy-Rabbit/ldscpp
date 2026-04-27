#pragma once

#include "annotation.hpp"
#include "continuous.hpp"
#include "ld_score.hpp"
#include "munge.hpp"
#include "sumstats_commands.hpp"

namespace ldsc_cpp {

void print_help();
Args parse_args(int argc, char** argv);

}  // namespace ldsc_cpp
