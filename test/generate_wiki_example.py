#!/usr/bin/env python3
"""Generate a small LDSC wiki-style parity dataset."""

from __future__ import annotations

import gzip
import math
import os
import shutil
from pathlib import Path


ALLELES = [("A", "C"), ("A", "G"), ("C", "T"), ("G", "T")]


def two_sided_p(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


def write_sumstats(path: Path, trait_shift: float) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write("SNP\tA1\tA2\tOR\tP\tINFO\tFRQ\n")
        for i in range(1, 81):
            snp = f"rs{i}"
            a1, a2 = ALLELES[i % len(ALLELES)]
            ld = 5.0 + (i % 13) * 0.35
            z = 0.15 * math.sin(i * 0.37 + trait_shift) + math.sqrt(1.0 + 0.55 * ld)
            if i % 2 == 0:
                z = -z
            odds = math.exp(0.03 if z >= 0 else -0.03)
            p = max(two_sided_p(z), 1e-300)
            frq = 0.08 + (i % 17) * 0.02
            info = 0.95 + (i % 5) * 0.005
            out.write(f"{snp}\t{a1}\t{a2}\t{odds:.8f}\t{p:.12g}\t{info:.3f}\t{frq:.3f}\n")


def write_merge_alleles(path: Path) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write("SNP\tA1\tA2\n")
        for i in range(1, 81):
            a1, a2 = ALLELES[i % len(ALLELES)]
            out.write(f"rs{i}\t{a1}\t{a2}\n")


def write_ldscores(data_dir: Path) -> None:
    ld_dir = data_dir / "eur_w_ld_chr"
    ld_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(ld_dir / "1.l2.ldscore.gz", "wt", encoding="utf-8") as out:
        out.write("CHR\tSNP\tBP\tL2\n")
        for i in range(1, 81):
            ld = 5.0 + (i % 13) * 0.35
            out.write(f"1\trs{i}\t{100000 + i * 1000}\t{ld:.6f}\n")
    (ld_dir / "1.l2.M").write_text("80\n", encoding="utf-8")
    (ld_dir / "1.l2.M_5_50").write_text("80\n", encoding="utf-8")


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


def write_py_wrappers(run_dir: Path, repo_root: Path) -> None:
    ldsc_root = repo_root / "ldsc"
    for script in ["ldsc.py", "munge_sumstats.py"]:
        py_wrapper = run_dir / "py" / script
        py_wrapper.write_text(
            "#!/usr/bin/env python\n"
            "import runpy, sys\n"
            f"sys.path.insert(0, {str(ldsc_root)!r})\n"
            f"runpy.run_path({str(ldsc_root / script)!r}, run_name='__main__')\n",
            encoding="utf-8",
        )
        py_wrapper.chmod(0o755)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    root = repo_root / "ldsc_cpp"
    run_dir = root / "test" / "wiki_flow" / "run"
    if run_dir.exists():
        shutil.rmtree(run_dir)
    for side in ["py", "cpp"]:
        (run_dir / side / "data").mkdir(parents=True)
        (run_dir / side / "out").mkdir()

    data = run_dir / "py" / "data"
    write_sumstats(data / "pgc.cross.SCZ17.2013-05.txt", 0.0)
    write_sumstats(data / "pgc.cross.BIP11.2013-05.txt", 0.9)
    write_merge_alleles(data / "w_hm3.snplist")
    write_ldscores(data)
    write_plink_fixture(data, repo_root)
    shutil.copytree(data, run_dir / "cpp" / "data", dirs_exist_ok=True)
    write_py_wrappers(run_dir, repo_root)

    commands_py = """\
./munge_sumstats.py --sumstats data/pgc.cross.SCZ17.2013-05.txt --N 1000 --out out/scz --merge-alleles data/w_hm3.snplist
./munge_sumstats.py --sumstats data/pgc.cross.BIP11.2013-05.txt --N 1000 --out out/bip --merge-alleles data/w_hm3.snplist
./ldsc.py --rg out/scz.sumstats.gz,out/bip.sumstats.gz --ref-ld-chr data/eur_w_ld_chr/ --w-ld-chr data/eur_w_ld_chr/ --out out/rg
./ldsc.py --h2 out/scz.sumstats.gz --ref-ld-chr data/eur_w_ld_chr/ --w-ld-chr data/eur_w_ld_chr/ --out out/h2
./ldsc.py --bfile data/plink/22 --l2 --ld-wind-cm 1 --out out/22
"""
    cpp_ldsc = Path(os.path.relpath(root / "build" / "ldsc", run_dir / "cpp")).as_posix()
    commands_cpp = f"""\
{cpp_ldsc} --sumstats data/pgc.cross.SCZ17.2013-05.txt --N 1000 --out out/scz --merge-alleles data/w_hm3.snplist
{cpp_ldsc} --sumstats data/pgc.cross.BIP11.2013-05.txt --N 1000 --out out/bip --merge-alleles data/w_hm3.snplist
{cpp_ldsc} --rg out/scz.sumstats.gz,out/bip.sumstats.gz --ref-ld-chr data/eur_w_ld_chr/ --w-ld-chr data/eur_w_ld_chr/ --out out/rg
{cpp_ldsc} --h2 out/scz.sumstats.gz --ref-ld-chr data/eur_w_ld_chr/ --w-ld-chr data/eur_w_ld_chr/ --out out/h2
{cpp_ldsc} --bfile data/plink/22 --l2 --ld-wind-cm 1 --out out/22
"""
    flow_dir = root / "test" / "wiki_flow"
    flow_dir.mkdir(parents=True, exist_ok=True)
    (flow_dir / "commands_py.sh").write_text(commands_py, encoding="utf-8")
    (flow_dir / "commands_cpp.sh").write_text(commands_cpp, encoding="utf-8")


if __name__ == "__main__":
    main()
