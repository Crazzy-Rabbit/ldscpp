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
