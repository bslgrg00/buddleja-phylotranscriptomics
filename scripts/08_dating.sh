#!/bin/bash
# 08_dating.sh: divergence times with treePL (penalized likelihood) and BEAST 2 (100 SortaDate loci)
# Calibrations (uniform, TimeTree): Paulowniaceae + Scrophulariaceae 64.4-67.6 Ma;
#                                   Scrophularia + Buddleja 43.9-47.1 Ma
# Output (Zenodo 09_dating): best_ml_clean.tre, RAxML_bootstrap.Astral_treepl_M20_bs,
#   buddleja_cv.txt, buddleja_dated_summary_treePL.tre, sortadate_100loci.txt,
#   buddleja_combined.log, buddleja_dated_summary_beast.tre
set -euo pipefail
THREADS=32
SRC=path/to/SortaDate/src
OUT=Rehmannia_glutinosa_OWAS,Lindenbergia_philippensis_WUZV

# 1. Branch lengths from the M20 supermatrix on the fixed ASTRAL-III topology (best tree and
#    1000 bootstrap replicates), rooted on Rehmannia and Lindenbergia
sed -E 's/\)[0-9.]+:/):/g' astral_M20.tre > astral_M20_clean.tre
raxmlHPC-PTHREADS -T $THREADS -f e -t astral_M20_clean.tre -s M20_supermatrix.phy \
                  -m GTRGAMMA -n treepl_M20
pxrr -t RAxML_result.treepl_M20 -g $OUT -o best_ml_clean.tre
raxmlHPC-PTHREADS -T $THREADS -g astral_M20_clean.tre -s M20_supermatrix.phy -m GTRGAMMA \
                  -p 12345 -b 12345 -# 1000 -n Astral_treepl_M20_bs
pxrr -t RAxML_bootstrap.Astral_treepl_M20_bs -g $OUT -o bootstrap_rooted.tre

# 2. treePL with the config files in Zenodo 09_dating: prime, cross-validation (optimal
#    smoothing 1e-8), and dating (thorough); bootstrap dates summarized as 95% CIs
treePL config_prime.txt
treePL config_cv.txt
treePL config_date.txt
treeannotator -burnin 0 -heights mean -target best_dated.tre bootstrap_dated.tre \
              buddleja_dated_summary_treePL.tre

# 3. SortaDate: root-to-tip variance, tree length, and bipartition support; top 100 loci
python $SRC/get_var_length.py genetrees_rooted/ --flend .tre --outf var --outg $OUT
python $SRC/get_bp_genetrees.py genetrees_rooted/ astral_M20_clean.tre --flend .tre --outf bp
python $SRC/combine_results.py var bp --outf comb
python $SRC/get_good_genes.py comb --max 100 --order 3,1,2 --outf sortadate_100loci.txt

# 4. Concatenated alignment of the 100 loci (154,844 bp) for BEAUti v2.6.3
#    (GTR+G4, uncorrelated lognormal clock, Yule, ASTRAL-III topology fixed) -> Buddleja_dating.xml
for G in $(awk 'NR > 1 {print $1}' sortadate_100loci.txt); do
  ls 09_M20/final/${G%%.*}.*aln-cln
done > beast_loci.txt
pxcat -s $(cat beast_loci.txt) -p beast_100loci.model -o beast_100loci.fa

# 5. BEAST: four runs of 80 million generations (sampled every 8000), 10% burn-in,
#    combined log and MCC tree with mean node heights and 95% HPD
for i in 1 2 3 4; do
  mkdir -p run$i
  (cd run$i && beast -beagle -beagle_CPU -beagle_SSE -beagle_double -threads 16 -seed $i \
                     -overwrite ../Buddleja_dating.xml)
done
logcombiner -b 10 $(printf -- '-log %s ' run*/*.log) -o buddleja_combined.log
logcombiner -b 10 $(printf -- '-log %s ' run*/*.trees) -o buddleja_combined.trees
treeannotator -burnin 0 -heights mean buddleja_combined.trees buddleja_dated_summary_beast.tre
