#!/bin/bash
# 06_gene_flow.sh: phylogenetic networks (PhyloNet InferNetwork_MPL, h = 0-3) and HyDe
#                  (NANUQ is in 05_MSCquartets.R)
# Output (Zenodo 07_gene_flow): ret*_output.txt, supermatrix_ingroup29.phy, hyde_map29.txt,
#                               buddleja_scan-out.txt
set -euo pipefail
PHYLONET=path/to/PhyloNetv3_8_2.jar
HYDE=path/to/HyDe/scripts/run_hyde.py

# 1. PhyloNet input: 21 taxa of Clades II and III plus Paulownia (UMUL); gene trees with
#    all 22 taxa, pruned to these taxa and rooted on UMUL
printf '%s\n' Buddleja_colvilei_AB1 Buddleja_sessilifolia_AB2 Buddleja_davidii_AB7 \
  Buddleja_forrestii_AB9 Buddleja_myriantha_AB12 Buddleja_albiflora_AB13 \
  Buddleja_fallowiana_AB14 Buddleja_candida_AB15 Buddleja_macrostachya_AB26 \
  Buddleja_calostachya_AB29 Buddleja_delavayi_AB3 Buddleja_yunnanensis_AB4 Buddleja_crispa_AB8 \
  Buddleja_microstachya_AB10 Buddleja_caryopteridifolia_AB19 Buddleja_curviflora_AB20 \
  Buddleja_lindleyana_AB22 Buddleja_officinalis_AB24 Buddleja_paniculata_AB25 \
  Buddleja_atropurpurea_AB28 Buddleja_purdomii_AB31 Paulownia_fargesii_UMUL > phylonet_taxa22.txt
Rscript -e '
library(ape); taxa <- scan("phylonet_taxa22.txt", "", quiet = TRUE)
gt <- Filter(function(t) all(taxa %in% t$tip.label), read.tree("M20_genetrees_rooted.tre"))
gt <- lapply(gt, function(t) { t <- keep.tip(t, taxa); t$node.label <- NULL
  root(t, "Paulownia_fargesii_UMUL", resolve.root = TRUE) })
class(gt) <- "multiPhylo"; write.tree(gt, "phylonet_genetrees22.tre")
cat(length(gt), "gene trees\n")'

# 2. PhyloNet: 0-3 reticulations, 20 runs each, 4 networks returned
for h in 0 1 2 3; do
  { echo -e "#NEXUS\nBEGIN TREES;"
    awk '{print "Tree gt" NR " = " $0}' phylonet_genetrees22.tre
    echo -e "END;\nBEGIN PHYLONET;\nInferNetwork_MPL (all) $h -x 20 -n 4 -pl 8 -di;\nEND;"
  } > ret${h}_buddleja22_phylonet.nex
  java -Xmx200g -jar $PHYLONET ret${h}_buddleja22_phylonet.nex > ret${h}_output.txt
done

# 3. HyDe: 28 Asian taxa plus OUT1 (outgroup), all sites of the M20 supermatrix
{ echo "29 2421192"
  awk 'NR > 1 && $1 !~ /_(OUT2|SNOD|UMUL|OWAS|WUZV)$/ {print $1 "\t" $2}' M20_supermatrix.phy
} > supermatrix_ingroup29.phy
awk 'NR > 1 {print $1 "\t" $1}' supermatrix_ingroup29.phy > hyde_map29.txt
python $HYDE -i supermatrix_ingroup29.phy -m hyde_map29.txt -o Buddleja_spp_OUT1 \
       -n 29 -t 29 -s 2421192 --prefix buddleja_scan --ignore_amb_sites
