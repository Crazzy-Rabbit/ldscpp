#!/usr/bin/env python3
"""Generate broader py2-vs-C++ LDSC parity fixtures."""

from __future__ import annotations

import bz2
import gzip
import math
import os
import shutil
from pathlib import Path

from generate_wiki_example import ALLELES, two_sided_p, write_ldscores, write_merge_alleles, write_sumstats


def write_py_wrappers(run_dir: Path, repo_root: Path) -> None:
    ldsc_root = repo_root / "ldsc"
    for script in ["ldsc.py", "munge_sumstats.py", "make_annot.py"]:
        py_wrapper = run_dir / "py" / script
        py_wrapper.write_text(
            "#!/usr/bin/env python\n"
            "import runpy, sys\n"
            f"sys.path.insert(0, {str(ldsc_root)!r})\n"
            f"runpy.run_path({str(ldsc_root / script)!r}, run_name='__main__')\n",
            encoding="utf-8",
        )
        py_wrapper.chmod(0o755)


def write_plink_fixture(data_dir: Path, repo_root: Path) -> None:
    plink_dir = data_dir / "plink"
    plink_dir.mkdir(parents=True, exist_ok=True)
    src = repo_root / "ldsc" / "test" / "plink_test" / "plink"
    shutil.copyfile(src.with_suffix(".bed"), plink_dir / "22.bed")
    shutil.copyfile(src.with_suffix(".fam"), plink_dir / "22.fam")
    lines = []
    for idx, line in enumerate(src.with_suffix(".bim").read_text(encoding="utf-8").splitlines()):
        parts = line.split()
        parts[2] = f"{idx * 0.4:.1f}"
        lines.append("\t".join(parts))
    (plink_dir / "22.bim").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_beta_sumstats(path: Path) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write("SNP\tA1\tA2\tBETA\tP\tN\tINFO\tFRQ\n")
        for i in range(1, 61):
            a1, a2 = ALLELES[i % len(ALLELES)]
            z = 1.1 + 0.025 * (i % 11)
            if i % 2 == 0:
                z = -z
            beta = 0.02 if z >= 0 else -0.02
            p = max(two_sided_p(z), 1e-300)
            info = 0.96 + (i % 3) * 0.01
            frq = 0.12 + (i % 12) * 0.03
            out.write(f"rs{i}\t{a1}\t{a2}\t{beta:.6f}\t{p:.12g}\t1000\t{info:.3f}\t{frq:.3f}\n")


def write_noalleles_sumstats(path: Path) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write("SNP\tBETA\tP\tN\n")
        for i in range(1, 41):
            z = 0.9 + 0.02 * (i % 9)
            if i % 2 == 0:
                z = -z
            beta = 0.01 if z >= 0 else -0.01
            out.write(f"rs{i}\t{beta:.6f}\t{two_sided_p(z):.12g}\t900\n")


def write_ldscore_annotations(data_dir: Path) -> None:
    bim = (data_dir / "plink" / "22.bim").read_text(encoding="utf-8").splitlines()
    with (data_dir / "plink" / "22.annot").open("w", encoding="utf-8") as out:
        out.write("CHR\tSNP\tBP\tCM\tA\tB\n")
        for idx, line in enumerate(bim):
            chrom, snp, cm, bp, *_ = line.split()
            out.write(f"{chrom}\t{snp}\t{bp}\t{cm}\t{1 if idx % 2 == 0 else 0}\t{1 if idx >= 4 else 0}\n")
    with (data_dir / "plink" / "22.thin.annot").open("w", encoding="utf-8") as out:
        out.write("A\tB\n")
        for idx, _line in enumerate(bim):
            out.write(f"{1 if idx % 2 == 0 else 0}\t{1 if idx >= 4 else 0}\n")
    with (data_dir / "cts.txt").open("w", encoding="utf-8") as out:
        for idx, line in enumerate(bim):
            snp = line.split()[1]
            value = 0.1 + idx * 0.12
            out.write(f"{snp}\t{value:.2f}\n")
    with bz2.open(data_dir / "print_snps.txt.bz2", "wt", encoding="utf-8") as out:
        out.write("rs_4\nrs_5\nrs_6\n")


def write_make_annot_inputs(data_dir: Path) -> None:
    (data_dir / "regions.bed").write_text("chr1\t4\t7\nchr1\t6\t9\n", encoding="utf-8")


def main() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    root = repo_root / "ldsc_cpp"
    run_dir = root / "test" / "full_flow" / "run"
    if run_dir.exists():
        shutil.rmtree(run_dir)
    for side in ["py", "cpp"]:
        (run_dir / side / "data").mkdir(parents=True)
        (run_dir / side / "out").mkdir()

    data = run_dir / "py" / "data"
    write_sumstats(data / "pgc.cross.SCZ17.2013-05.txt", 0.0)
    write_sumstats(data / "pgc.cross.BIP11.2013-05.txt", 0.9)
    write_beta_sumstats(data / "beta.sumstats")
    write_noalleles_sumstats(data / "noalleles.sumstats")
    write_merge_alleles(data / "w_hm3.snplist")
    write_ldscores(data)
    write_plink_fixture(data, repo_root)
    write_ldscore_annotations(data)
    write_make_annot_inputs(data)
    shutil.copyfile(repo_root / "ldsc" / "test" / "munge_test" / "sumstats", data / "daner.sumstats")
    shutil.copytree(data, run_dir / "cpp" / "data", dirs_exist_ok=True)
    write_py_wrappers(run_dir, repo_root)

    commands_py = """\
./make_annot.py --bed-file data/regions.bed --bimfile data/plink/22.bim --annot-file out/make.annot
./munge_sumstats.py --sumstats data/pgc.cross.SCZ17.2013-05.txt --N 1000 --out out/scz --merge-alleles data/w_hm3.snplist
./munge_sumstats.py --sumstats data/pgc.cross.BIP11.2013-05.txt --N 1000 --out out/bip --merge-alleles data/w_hm3.snplist
./munge_sumstats.py --sumstats data/beta.sumstats --out out/beta --keep-maf
./munge_sumstats.py --sumstats data/noalleles.sumstats --out out/noalleles --no-alleles --signed-sumstats BETA,0
./munge_sumstats.py --sumstats data/daner.sumstats --out out/daner --daner
./ldsc.py --h2 out/scz.sumstats.gz --ref-ld-chr data/eur_w_ld_chr/ --w-ld-chr data/eur_w_ld_chr/ --out out/h2_extra --print-cov --print-delete-vals --M 80 --n-blocks 10 --chisq-max 999.0
./ldsc.py --rg out/scz.sumstats.gz,out/bip.sumstats.gz --ref-ld-chr data/eur_w_ld_chr/ --w-ld-chr data/eur_w_ld_chr/ --out out/rg_extra --print-cov --print-delete-vals --n-blocks 10
./ldsc.py --bfile data/plink/22 --l2 --ld-wind-snps 1 --out out/annot --annot data/plink/22.annot
./ldsc.py --bfile data/plink/22 --l2 --ld-wind-snps 1 --out out/thin --annot data/plink/22.thin.annot --thin-annot
./ldsc.py --bfile data/plink/22 --l2 --ld-wind-snps 1 --out out/cts --cts-bin data/cts.txt --cts-breaks 0.25,0.55 --cts-names C
./ldsc.py --bfile data/plink/22 --l2 --ld-wind-snps 1 --out out/print --print-snps data/print_snps.txt.bz2
./ldsc.py --bfile data/plink/22 --l2 --ld-wind-snps 1 --out out/pq --pq-exp 1.0
"""
    cpp_ldsc = Path(os.path.relpath(root / "build" / "ldsc", run_dir / "cpp")).as_posix()
    commands_cpp = f"""\
{cpp_ldsc} --bed-file data/regions.bed --bimfile data/plink/22.bim --annot-file out/make.annot
{cpp_ldsc} --sumstats data/pgc.cross.SCZ17.2013-05.txt --N 1000 --out out/scz --merge-alleles data/w_hm3.snplist
{cpp_ldsc} --sumstats data/pgc.cross.BIP11.2013-05.txt --N 1000 --out out/bip --merge-alleles data/w_hm3.snplist
{cpp_ldsc} --sumstats data/beta.sumstats --out out/beta --keep-maf
{cpp_ldsc} --sumstats data/noalleles.sumstats --out out/noalleles --no-alleles --signed-sumstats BETA,0
{cpp_ldsc} --sumstats data/daner.sumstats --out out/daner --daner
{cpp_ldsc} --h2 out/scz.sumstats.gz --ref-ld-chr data/eur_w_ld_chr/ --w-ld-chr data/eur_w_ld_chr/ --out out/h2_extra --print-cov --print-delete-vals --M 80 --n-blocks 10 --chisq-max 999.0
{cpp_ldsc} --rg out/scz.sumstats.gz,out/bip.sumstats.gz --ref-ld-chr data/eur_w_ld_chr/ --w-ld-chr data/eur_w_ld_chr/ --out out/rg_extra --print-cov --print-delete-vals --n-blocks 10
{cpp_ldsc} --bfile data/plink/22 --l2 --ld-wind-snps 1 --out out/annot --annot data/plink/22.annot
{cpp_ldsc} --bfile data/plink/22 --l2 --ld-wind-snps 1 --out out/thin --annot data/plink/22.thin.annot --thin-annot
{cpp_ldsc} --bfile data/plink/22 --l2 --ld-wind-snps 1 --out out/cts --cts-bin data/cts.txt --cts-breaks 0.25,0.55 --cts-names C
{cpp_ldsc} --bfile data/plink/22 --l2 --ld-wind-snps 1 --out out/print --print-snps data/print_snps.txt.bz2
{cpp_ldsc} --bfile data/plink/22 --l2 --ld-wind-snps 1 --out out/pq --pq-exp 1.0
"""
    flow_dir = root / "test" / "full_flow"
    flow_dir.mkdir(parents=True, exist_ok=True)
    (flow_dir / "commands_py.sh").write_text(commands_py, encoding="utf-8")
    (flow_dir / "commands_cpp.sh").write_text(commands_cpp, encoding="utf-8")


if __name__ == "__main__":
    main()
