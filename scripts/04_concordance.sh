#!/bin/bash
# 04_concordance.sh: Phyparts (gene tree concordance and conflict) and gene and site
#                    concordance factors (M20 dataset)
# Output (Zenodo 05_concordance): Buddleja70_phyparts.*, phyparts_pies_asian_buddleja70.svg,
#                                 Buddleja_{ASTRAL,IQTREE,RAxML}_CF.cf.stat
set -euo pipefail
PHYPARTS=path/to/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar
PIES=path/to/phypartspiecharts.py
OUT=Rehmannia_glutinosa_OWAS,Lindenbergia_philippensis_WUZV

# 1. M20 MO ortholog trees, tip names reduced to taxon names, rooted on Rehmannia and
#    Lindenbergia; trees that cannot be rooted are dropped
mkdir -p genetrees_rooted
for f in 08_MO/M20/*.tre; do
  sed -E 's/@[^:,()]+//g' $f | pxrr -g $OUT -o genetrees_rooted/$(basename $f) || true
done
find genetrees_rooted -size 0 -delete
NGENES=$(ls genetrees_rooted | wc -l)
cat genetrees_rooted/* > M20_genetrees_rooted.tre

# 2. Rooted ASTRAL-III species tree
pxrr -t astral_M20.tre -g $OUT -o re-astral_M20.tre

# 3. Phyparts (gene tree branches with BS < 70 collapsed) and pie charts (Figure 2A)
java -jar $PHYPARTS -a 1 -d genetrees_rooted -m re-astral_M20.tre -o Buddleja70_phyparts -s 70
python $PIES re-astral_M20.tre Buddleja70_phyparts $NGENES \
       --svg_name phyparts_pies_asian_buddleja70.svg

# 4. Gene and site concordance factors on the ASTRAL-III, IQ-TREE, and RAxML-NG trees
for T in ASTRAL:astral_M20.tre IQTREE:M20_IQTREE.treefile RAxML:M20_RAxML.support; do
  iqtree2 -t ${T#*:} -s M20_supermatrix.phy --gcf M20_genetrees_BS.tre --scf 100 -T 4 \
          --prefix Buddleja_${T%%:*}_CF
done
