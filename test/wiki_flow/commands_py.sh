./munge_sumstats.py --sumstats data/pgc.cross.SCZ17.2013-05.txt --N 1000 --out out/scz --merge-alleles data/w_hm3.snplist
./munge_sumstats.py --sumstats data/pgc.cross.BIP11.2013-05.txt --N 1000 --out out/bip --merge-alleles data/w_hm3.snplist
./ldsc.py --rg out/scz.sumstats.gz,out/bip.sumstats.gz --ref-ld-chr data/eur_w_ld_chr/ --w-ld-chr data/eur_w_ld_chr/ --out out/rg
./ldsc.py --h2 out/scz.sumstats.gz --ref-ld-chr data/eur_w_ld_chr/ --w-ld-chr data/eur_w_ld_chr/ --out out/h2
./ldsc.py --bfile data/plink/22 --l2 --ld-wind-cm 1 --out out/22
