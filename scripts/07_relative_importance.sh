#!/bin/bash
# 07_relative_importance.sh: contributions of ILS, gene flow, and gene tree estimation error to
#   gene tree variation, following Cai et al. (2021)
# Input trees (rooted on OWAS + WUZV, no node labels, all branch lengths > 0):
#   species.tre      ASTRAL-III M20 tree, branch lengths in coalescent units
#   bs_species.tre   100 ASTRAL-III bootstrap species trees
#   genetrees.tre    M20 gene trees;  BSgenetrees.tre  M20 bootstrap gene trees
# Output (Zenodo 08_relative_importance): theta.csv, Buddleja_reticulation_index.csv,
#   Buddleja_estimation_error.csv, Buddleja_gCF.csv, relative_contribution.csv,
#   Buddleja_relative_contribution_with_CI.csv
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
CAI=path/to/Coalescent_simulation_and_gene_flow_detection

# 1. Gene flow: gene trees simulated under the MSC (Phybase), triplet frequencies, and
#    reticulation index per node (unbalancedTriplet.perc.tre)
sh $CAI/simulator_tripletCounter_tripletMapper.sh species.tre bs_species.tre genetrees.tre

# 2. ILS: branch lengths in substitutions on the fixed ASTRAL-III topology (GTR+G);
#    theta = bl.mut / bl.coal
iqtree2 -s M20_supermatrix.phy -te species.tre -m GTR+G -T AUTO --prefix theta

# 3. Estimation error: 100 alignments of 1135 bp (Seq-Gen, GTR+G parameters from IQ-TREE),
#    ML gene trees, and recovery (%) of each species tree node
RATES=$(grep -E '^  [ACGT]-[ACGT]:' theta.iqtree | awk '{print $2}' | paste -sd,)
FREQS=$(grep -E '^  pi\([ACGT]\)' theta.iqtree | awk '{print $3}' | paste -sd,)
ALPHA=$(grep 'Gamma shape alpha' theta.iqtree | awk '{print $NF}')
seq-gen -m GTR -r $RATES -f $FREQS -a $ALPHA -l 1135 -n 100 -op < theta.treefile > sim.phy
mkdir -p sim && awk '/^ *[0-9]+ +1135 *$/ {n++} {print > sprintf("sim/sim%03d.phy", n)}' sim.phy
for f in sim/*.phy; do
  raxmlHPC -m GTRGAMMA -p 12345 -s $f -n $(basename $f .phy) -w $PWD/sim > /dev/null
done
cat sim/RAxML_bestTree.* > sim_gene.trees
raxmlHPC -f b -t species.tre -z sim_gene.trees -m GTRGAMMA -n ERR

# 4. Gene tree variation: gCF from bootstrap gene trees (Gene_tree_variation = 100 - gCF)
iqtree2 -t species.tre --gcf BSgenetrees.tre --prefix concord

# 5. Node values combined in relative_contribution.csv (zero estimation error set to 0.5),
#    then relative importance (LMG, relaimpo)
Rscript $REPO/scripts/07_relaimpo.R
