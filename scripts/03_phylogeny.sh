#!/bin/bash
# 03_phylogeny.sh: gene tree sets, concatenation trees (IQ-TREE, RAxML-NG), coalescent trees
#                  (ASTRAL-III, ASTRAL-Pro3), collapsed gene trees, and RF distances
# Output (Zenodo 03_gene_trees): M*_genetrees_BS.tre, M20_genetrees_BS*.tre
# Output (Zenodo 04_species_trees): M*_IQTREE.treefile, M*_RAxML.support, astral_M*.tre,
#                                   astral_pro3_M*.tre, astral_M20_BS*.tre
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
THREADS=32
ASTRAL=path/to/astral.5.7.8.jar

for PCT in 10 20 30; do
  D=09_M$PCT

  # 1. Gene trees: MO ortholog trees (RAxML GTRCAT, 200 BS) of loci that passed the PHI test,
  #    in alignment order, tip names reduced to taxon names
  for ALN in $D/final/*.aln-cln; do
    ID=$(basename $ALN); ID=${ID%%.*}
    cat 08_MO/M$PCT/$ID.*tre
  done | sed -E 's/@[^:,()]+//g' > M${PCT}_genetrees_BS.tre

  # 2. IQ-TREE (ModelFinder per partition, 1000 UFBoot with -bnni, 1000 SH-aLRT)
  iqtree2 -s M${PCT}_supermatrix.phy -p M${PCT}_supermatrix.model -st DNA -m MFP \
          -B 1000 -bnni -alrt 1000 -T AUTO --prefix M${PCT}_IQTREE

  # 3. RAxML-NG (GTR+G per partition, 100 BS)
  sed 's/^DNA,/GTR+G,/' M${PCT}_supermatrix.model > M${PCT}_raxmlng.model
  raxml-ng --all --msa M${PCT}_supermatrix.phy --model M${PCT}_raxmlng.model \
           --bs-trees 100 --seed 12345 --threads $THREADS --workers 1 --prefix M${PCT}_RAxML
  mv M${PCT}_RAxML.raxml.support M${PCT}_RAxML.support

  # 4. ASTRAL-III (local posterior probabilities) and ASTRAL-Pro3
  java -jar $ASTRAL -i M${PCT}_genetrees_BS.tre -o astral_M$PCT.tre 2> astral_M$PCT.log
  astral-pro3 -i M${PCT}_genetrees_BS.tre -o astral_pro3_M$PCT.tre -t $THREADS -u 1 \
              2> astral_pro3_M$PCT.log
done

# 5. M20 gene trees with branches of BS <= 10, 20, 30, and 33 collapsed, then ASTRAL-III
for BS in 10 20 30 33; do
  nw_ed M20_genetrees_BS.tre "i & b<=$BS" o > M20_genetrees_BS$BS.tre
  java -jar $ASTRAL -i M20_genetrees_BS$BS.tre -o astral_M20_BS$BS.tre 2> astral_M20_BS$BS.log
done

# 6. Collapsed vs uncollapsed ASTRAL-III trees (compare_RF.tsv, compare_key_nodes.tsv)
python $REPO/helpers/compare_astral_trees.py BS0=astral_M20.tre BS10=astral_M20_BS10.tre \
        BS20=astral_M20_BS20.tre BS30=astral_M20_BS30.tre BS33=astral_M20_BS33.tre

# 7. RF distances among the 12 species trees
Rscript $REPO/scripts/03_RF_distances.R
