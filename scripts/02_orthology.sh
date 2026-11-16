#!/bin/bash
# 02_orthology.sh: homology search, homolog trees, tree cleaning, orthologs (monophyletic
#                  outgroup method), PHI test, and supermatrices (M10, M20, M30)
# Output (Zenodo 02_orthologs): 09_M*/final/*.fa.mafft.aln-cln, M*_supermatrix.phy, M*_supermatrix.model
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
THREADS=32
PDC=path/to/phylogenomic_dataset_construction/scripts
export PYTHONPATH=$PDC:${PYTHONPATH:-}      # helpers import Yang and Smith modules

mkdir -p 06_blast 06_clusters 07_trees 07_treeshrink 07_masked 07_cut 08_MO

# 1. All-by-all BLASTN of non-redundant CDS
cat 05_cdhit/*.cds.cdhit > 06_blast/all.fa
awk '/^>/{f=!seen[$1]++} f' 06_blast/all.fa > 06_blast/all_dedup.fa
makeblastdb -in 06_blast/all_dedup.fa -parse_seqids -dbtype nucl -out 06_blast/all_dedup.fa
for FA in 05_cdhit/*.cds.cdhit; do
  blastn -db 06_blast/all_dedup.fa -query $FA -evalue 10 -max_target_seqs 1000 \
         -num_threads $THREADS \
         -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore' \
         >> 06_blast/all.rawblast
done

# 2. MCL clustering (hit fraction 0.4, inflation 1.4); clusters with >= 30 taxa
python $PDC/blast_to_mcl.py 06_blast/all.rawblast 0.4
mcl 06_blast/all.rawblast.hit-frac0.4.minusLogEvalue --abc -te $THREADS -tf 'gq(5)' \
    -I 1.4 -o 06_blast/hit-frac0.4_I1.4.mcl
ls 05_cdhit/*.cds.cdhit > 06_blast/cds_files.txt
python $PDC/write_fasta_files_from_mcl.py 06_blast/cds_files.txt \
        06_blast/hit-frac0.4_I1.4.mcl 30 06_clusters

# 3. Homolog trees: MAFFT, pxclsq (min column occupancy 0.1), RAxML (GTRCAT, 200 BS)
for FA in 06_clusters/*.fa; do
  N=$(basename $FA .fa)
  mafft --genafpair --maxiterate 1000 --thread $THREADS $FA > 07_trees/$N.fa.mafft.aln
  pxclsq -p 0.1 -s 07_trees/$N.fa.mafft.aln -o 07_trees/$N.fa.mafft.aln-cln
  (cd 07_trees && raxmlHPC-PTHREADS-AVX2 -T $THREADS -f a -x 12345 -p 12345 -# 200 \
                  -m GTRCAT -s $N.fa.mafft.aln-cln -n $N)
  mv 07_trees/RAxML_bipartitions.$N 07_trees/$N.raxml_bs.tre
done

# 4. Tree cleaning: TreeShrink (quantile 0.05), mask mono- and paraphyletic tips,
#    cut deep paralogs (internal branches > 0.5; keep subtrees with >= 17 taxa)
python $PDC/tree_shrink_wrapper_update.py 07_trees raxml_bs.tre 0.05 07_treeshrink
python $REPO/helpers/restore_bootstrap_support.py 07_trees 07_treeshrink raxml_bs.tre tt tt_bs.tre
for f in 07_treeshrink/*.tt_bs.tre; do cp $f 07_masked/$(basename $f .tt_bs.tre).tt; done
python $PDC/mask_tips_by_taxonID_transcripts.py 07_masked 07_trees y
python $PDC/cut_long_internal_branches.py 07_masked tt.mm 0.5 17 07_cut

for TIER in 10:27 20:24 30:21; do       # max missing taxa (%) : min ingroup taxa (of 30)
  PCT=${TIER%%:*}; MIN=${TIER##*:}; D=09_M$PCT
  mkdir -p 08_MO/M$PCT $D/fasta $D/aln $D/phi $D/final

  # 5. Orthologs (monophyletic outgroup method)
  python $PDC/prune_paralogs_MO.py 07_cut subtree $MIN 08_MO/M$PCT $REPO/data/taxon_file.txt

  # 6. Ortholog alignments: MAFFT; pxclsq removes columns with > 70% missing data
  python $PDC/write_ortholog_fasta_files.py 06_blast/all.fa 08_MO/M$PCT $D/fasta $MIN
  for FA in $D/fasta/*.fa; do
    mafft --genafpair --maxiterate 1000 --thread $THREADS $FA \
          > $D/aln/$(basename $FA .fa).fa.mafft.aln
  done
  python $PDC/pxclsq_wrapper.py $D/aln 0.3 dna

  # 7. PHI test (PhiPack, default settings); loci with P <= 0.05 removed,
  #    loci with too few informative sites for the test kept
  for ALN in $D/aln/*.aln-cln; do
    N=$(basename $ALN)
    mkdir -p $D/phi/$N
    (cd $D/phi/$N && Phi -f ../../aln/$N > /dev/null 2>&1 || true)
    PV=$(grep 'PHI (Normal)' $D/phi/$N/Phi.log 2>/dev/null | awk '{print $NF}' || true)
    echo -e "$N\t${PV:-NA}" >> $D/phi_results.txt
    if [[ -z "$PV" ]] || awk -v p="$PV" 'BEGIN {exit !(p > 0.05)}'; then
      cp $ALN $D/final/
    fi
  done

  # 8. Supermatrix (loci with >= 100 sites) and partition file
  python $PDC/concatenate_matrices_phyx.py $D/final 100 $MIN $D/supermatrix_MO_mintaxa$MIN
  cp $D/supermatrix_MO_mintaxa$MIN.phy M${PCT}_supermatrix.phy
  cp $D/supermatrix_MO_mintaxa$MIN.model M${PCT}_supermatrix.model
done
