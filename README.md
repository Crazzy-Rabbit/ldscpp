# ldscpp

`ldscpp` is a C++17 rewrite of the Python2 [LDSC](https://github.com/bulik/ldsc) command-line tools. It kept LDSC inputs, options, output files, and numerical results compatible with the original Python implementation.

The C++ version is invoked as:

```bash
build/ldsc [same LDSC options as Python LDSC]
```

For compatibility testing, the first command token is allowed to differ (`./ldsc.py` or `./munge_sumstats.py` on the Python side, `build/ldsc` on the C++ side). All LDSC parameters after the program name are checked for exact equality.

## Implemented Features

The current C++ implementation includes:

- Summary-statistics munging equivalent to `munge_sumstats.py`.
- Heritability regression with `--h2`.
- Cell-type-specific regression with `--h2-cts`.
- Genetic correlation with `--rg`.
- LD Score estimation from PLINK `.bed/.bim/.fam` files with `--bfile ... --l2`.
- Annotation-aware LD Scores with `--annot`, `--thin-annot`, and `--cts-bin`.
- `--print-snps`, including bzip2 input.
- `--pq-exp` and `--per-allele` LD Score scaling.
- `--print-cov`, `--print-delete-vals`, block jackknife output, and partitioned delete values.
- Annotation generation equivalent to `make_annot.py` using `--bed-file` or `--gene-set-file`.
- Continuous-annotation helpers equivalent to `quantile_h2g.r` and `quantile_M.pl`.
- gzip and bzip2 text input reading for LDSC data files.

## Build

Build on WSL/Linux with CMake:

```bash
cd ldsc_cpp
cmake -S . -B build
cmake --build build -j2
```

The executable is:

```text
ldsc_cpp/build/ldsc
```

The repository also contains a simple Windows Makefile:

```powershell
cd ldsc_cpp
mingw32-make
```

## Python2 Parity Checks

The test suite generates small LDSC fixtures, runs the original Python2 LDSC and this C++ rewrite, and compares:

- Command parameters, ignoring only the program name.
- Text and gzip output files.
- LDSC log content after normalizing timestamps, author/institution banner lines, executable path, and known pandas diagnostic table formatting.
- Numeric matrices such as `.cov`, `.delete`, and `.part_delete` with strict floating-point tolerance.

Run:

```bash
cd ldsc_cpp
bash test/run_wiki_flow_parity_wsl.sh
bash test/run_full_parity_wsl.sh
```

The latest check passed:

```text
OK command parameters
OK out/scz.sumstats.gz
OK out/bip.sumstats.gz
OK out/22.l2.ldscore.gz
OK out/22.l2.M
OK out/22.l2.M_5_50
OK out/scz.log
OK out/bip.log
OK out/rg.log
OK out/h2.log
OK out/22.log
```

The full parity check also passed all generated outputs, including:

- munge outputs for fixed `N`, per-SNP `N`, `--no-alleles`, and `--daner`.
- h2 and rg logs.
- h2 and rg covariance/delete matrices.
- PLINK LD Score outputs.
- annotation, thin annotation, cts-bin annotation, print-snps, and pq-exp outputs.
- make-annot output.

Example C++ log `Call:` block:

```text
../../../../build/ldsc \
--h2 out/scz.sumstats.gz \
--ref-ld-chr data/eur_w_ld_chr/ \
--w-ld-chr data/eur_w_ld_chr/ \
--out out/h2
```

That is intentional: C++ logs show the real compiled executable, not `./ldsc.py`.

## Speed Benchmark

Run the speed benchmark:

```bash
cd ldsc_cpp
bash test/benchmark_speed_wsl.sh 5
```

The benchmark uses the generated `test/full_flow` fixture, runs each Python2 and C++ command 5 times, and reports medians. Results are saved to:

- `test/speed/results.tsv`
- `test/speed/summary.tsv`

Latest benchmark:

| Analysis | Python2 median (s) | C++ median (s) | Speedup | Time reduction |
|---|---:|---:|---:|---:|
| Full flow | 10.680 | 0.737 | 14.50x | 93.1% |
| make annot | 0.856 | 0.029 | 29.73x | 96.6% |
| munge SCZ | 0.911 | 0.029 | 31.96x | 96.9% |
| munge BIP | 0.862 | 0.034 | 25.58x | 96.1% |
| munge BETA | 0.910 | 0.032 | 28.82x | 96.5% |
| munge no-alleles | 0.878 | 0.029 | 30.53x | 96.7% |
| munge daner | 0.875 | 0.028 | 31.32x | 96.8% |
| h2 regression | 0.821 | 0.139 | 5.92x | 83.1% |
| rg regression | 0.840 | 0.158 | 5.31x | 81.2% |
| LD Score with annot | 0.764 | 0.040 | 19.29x | 94.8% |
| LD Score thin annot | 0.754 | 0.037 | 20.13x | 95.0% |
| LD Score cts-bin | 0.760 | 0.040 | 18.87x | 94.7% |
| LD Score print-snps | 0.726 | 0.038 | 18.99x | 94.7% |
| LD Score pq-exp | 0.734 | 0.034 | 21.43x | 95.3% |

These numbers are from small generated fixtures. They are useful for regression testing and for showing the startup and command-flow overhead reduction. On large real LDSC datasets, speedup will depend on file I/O, matrix size, annotation count, SNP count, and regression settings.

## Command Migration

Python LDSC:

```bash
./munge_sumstats.py --sumstats input.txt --N 1000 --out out/trait
./ldsc.py --h2 out/trait.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out out/h2
./make_annot.py --bed-file regions.bed --bimfile 22.bim --annot-file 22.annot
```

C++ LDSC:

```bash
build/ldsc --sumstats input.txt --N 1000 --out out/trait
build/ldsc --h2 out/trait.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out out/h2
build/ldsc --bed-file regions.bed --bimfile 22.bim --annot-file 22.annot
```

Only the executable changes. LDSC options remain the same.

## Usage Examples

### 1. Munge Summary Statistics

Fixed sample size:

```bash
build/ldsc \
  --sumstats data/pgc.cross.SCZ17.2013-05.txt \
  --N 1000 \
  --out out/scz \
  --merge-alleles data/w_hm3.snplist
```

Per-SNP sample size from the input file:

```bash
build/ldsc \
  --sumstats data/beta.sumstats \
  --out out/beta \
  --keep-maf
```

No allele columns, using a signed statistic:

```bash
build/ldsc \
  --sumstats data/noalleles.sumstats \
  --out out/noalleles \
  --no-alleles \
  --signed-sumstats BETA,0
```

DANER input:

```bash
build/ldsc \
  --sumstats data/daner.sumstats \
  --out out/daner \
  --daner
```

### 2. Heritability Regression

Basic h2:

```bash
build/ldsc \
  --h2 out/scz.sumstats.gz \
  --ref-ld-chr data/eur_w_ld_chr/ \
  --w-ld-chr data/eur_w_ld_chr/ \
  --out out/h2
```

h2 with covariance and jackknife delete values:

```bash
build/ldsc \
  --h2 out/scz.sumstats.gz \
  --ref-ld-chr data/eur_w_ld_chr/ \
  --w-ld-chr data/eur_w_ld_chr/ \
  --out out/h2_extra \
  --print-cov \
  --print-delete-vals \
  --M 80 \
  --n-blocks 10 \
  --chisq-max 999.0
```

### 3. Genetic Correlation

Basic rg:

```bash
build/ldsc \
  --rg out/scz.sumstats.gz,out/bip.sumstats.gz \
  --ref-ld-chr data/eur_w_ld_chr/ \
  --w-ld-chr data/eur_w_ld_chr/ \
  --out out/rg
```

rg with covariance and delete values:

```bash
build/ldsc \
  --rg out/scz.sumstats.gz,out/bip.sumstats.gz \
  --ref-ld-chr data/eur_w_ld_chr/ \
  --w-ld-chr data/eur_w_ld_chr/ \
  --out out/rg_extra \
  --print-cov \
  --print-delete-vals \
  --n-blocks 10
```

### 4. LD Score Estimation

Basic LD Score:

```bash
build/ldsc \
  --bfile data/plink/22 \
  --l2 \
  --ld-wind-cm 1 \
  --out out/22
```

LD Score with full annotation:

```bash
build/ldsc \
  --bfile data/plink/22 \
  --l2 \
  --ld-wind-snps 1 \
  --out out/annot \
  --annot data/plink/22.annot
```

LD Score with thin annotation:

```bash
build/ldsc \
  --bfile data/plink/22 \
  --l2 \
  --ld-wind-snps 1 \
  --out out/thin \
  --annot data/plink/22.thin.annot \
  --thin-annot
```

LD Score with continuous-variable bins:

```bash
build/ldsc \
  --bfile data/plink/22 \
  --l2 \
  --ld-wind-snps 1 \
  --out out/cts \
  --cts-bin data/cts.txt \
  --cts-breaks 0.25,0.55 \
  --cts-names C
```

Print LD Scores only for selected SNPs:

```bash
build/ldsc \
  --bfile data/plink/22 \
  --l2 \
  --ld-wind-snps 1 \
  --out out/print \
  --print-snps data/print_snps.txt.bz2
```

LD Score with `pq-exp`:

```bash
build/ldsc \
  --bfile data/plink/22 \
  --l2 \
  --ld-wind-snps 1 \
  --out out/pq \
  --pq-exp 1.0
```

### 5. Make Annotation

From BED intervals:

```bash
build/ldsc \
  --bed-file data/regions.bed \
  --bimfile data/plink/22.bim \
  --annot-file out/make.annot
```

From a gene set:

```bash
build/ldsc \
  --gene-set-file genes.txt \
  --gene-coord-file gene_coords.txt \
  --bimfile data/plink/22.bim \
  --windowsize 100000 \
  --annot-file out/gene.annot
```

### 6. Cell-Type-Specific h2

```bash
build/ldsc \
  --h2-cts out/scz.sumstats.gz \
  --ref-ld-chr data/baselineLD. \
  --ref-ld-chr-cts data/cell_type_ldcts.txt \
  --w-ld-chr data/weights. \
  --out out/h2_cts
```

### 7. Continuous Annotation Helpers

Quantile h2g:

```bash
build/ldsc quantile-h2g \
  --annotfile annot.txt \
  --resultfile results.txt \
  --outfile out/quantile_h2g.txt \
  --nb-quantile 10
```

Quantile M:

```bash
build/ldsc quantile-m \
  --frqfile-chr data/frq/ \
  --ref-annot-chr data/annot/ \
  --annot-header ANNOT \
  --outfile out/quantile_m.txt \
  --nb-quantile 10
```

## Test Directory

- `test/generate_wiki_example.py`: creates the small wiki-style fixture.
- `test/generate_full_parity.py`: creates the broader full-flow fixture.
- `test/compare_wiki_flow.py`: compares wiki-flow parameters and outputs.
- `test/compare_full_parity.py`: compares full-flow parameters and outputs.
- `test/run_wiki_flow_parity_wsl.sh`: rebuilds and runs wiki parity.
- `test/run_full_parity_wsl.sh`: rebuilds and runs full parity.
- `test/benchmark_speed_wsl.sh`: runs the Python2/C++ speed benchmark.
- `test/wiki_flow/run`: latest wiki parity outputs.
- `test/full_flow/run`: latest full parity outputs.
- `test/speed`: latest benchmark outputs.

## Source Layout

- `src/ldsc.cpp`: main CLI entry, argument parsing, and command dispatch.
- `src/common.cpp`: shared I/O, table parsing, matrix utilities, jackknife, and IRWLS helpers.
- `src/regression.cpp`: h2, genetic covariance, and rg regression math.
- `src/sumstats_commands.cpp`: top-level `--h2`, `--h2-cts`, and `--rg` command flow.
- `src/munge.cpp`: C++ rewrite of summary-statistics munging.
- `src/annotation.cpp`: annotation generation.
- `src/ld_score.cpp`: PLINK reader and LD Score estimation.
- `src/continuous.cpp`: continuous annotation helper commands.
- `include/*.hpp`: headers for the C++ modules.
