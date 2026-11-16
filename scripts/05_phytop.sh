#!/bin/bash
# 05_phytop.sh: Phytop quartet analysis (ILS and hybridization indices) at each node of the
#               ASTRAL-III M20 tree
# Output (Zenodo 06_ILS_tests): phytop.info.tsv
set -euo pipefail
ASTRAL=path/to/astral.5.7.8.jar

# 1. ASTRAL-III with full branch annotation (-t 2, quartet frequencies)
java -jar $ASTRAL -i M20_genetrees_BS.tre -t 2 -o astral_M20_t2.tre 2> astral_M20_t2.log

# 2. Root on Rehmannia and Lindenbergia, then Phytop (Figure 2C)
pxrr -t astral_M20_t2.tre -g Rehmannia_glutinosa_OWAS,Lindenbergia_philippensis_WUZV \
     -o re-astral_M20_t2.tre
phytop re-astral_M20_t2.tre -pre phytop -align -figsize 4 -fontsize 11 \
       -branch_size 36 -leaf_size 40 -polytomy_test
