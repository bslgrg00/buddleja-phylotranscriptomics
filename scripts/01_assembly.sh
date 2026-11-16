#!/bin/bash
# 01_assembly.sh: read filtering, de novo assembly, transcript filtering, ORF prediction,
#                 and redundancy reduction for each accession
# Input:  reads/<code>_R1.fq.gz, reads/<code>_R2.fq.gz; databases/Atha.fa, databases/alternifolia.fa
# Output (Zenodo 01_assemblies): 04_transdecoder/<code>/<code>.cds.fa, <code>.pep.fa
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
THREADS=32
PDC=path/to/phylogenomic_dataset_construction/scripts

mkdir -p 01_qc 02_trinity 03_transrate 03_filtered 03_chimera 03_corset 03_reps \
         04_transdecoder 05_cdhit busco

# BLAST database for TransDecoder (Arabidopsis thaliana and Buddleja alternifolia proteomes)
cat databases/Atha.fa databases/alternifolia.fa > databases/db
makeblastdb -in databases/db -parse_seqids -dbtype prot -out databases/db

for ID in $(ls reads/*_R1.fq.gz | sed 's|.*/||; s|_R1.fq.gz||'); do
  Q1=01_qc/${ID}_R1.fq.gz; Q2=01_qc/${ID}_R2.fq.gz

  # 1. Read filtering
  fastp -i reads/${ID}_R1.fq.gz -I reads/${ID}_R2.fq.gz -o $Q1 -O $Q2 \
        -5 -3 -n 0 -l 60 -q 20 -u 30 -w $THREADS

  # 2. De novo assembly
  Trinity --seqType fq --left $Q1 --right $Q2 --max_memory 200G --CPU $THREADS \
          --output 02_trinity/${ID}_trinity
  cp 02_trinity/${ID}_trinity/Trinity.fasta 02_trinity/${ID}.Trinity.fasta

  # 3. TransRate filtering
  transrate --assembly 02_trinity/${ID}.Trinity.fasta --left $Q1 --right $Q2 \
            --threads $THREADS --output 03_transrate/${ID}
  python $PDC/filter_transcripts_transrate.py 02_trinity/${ID}.Trinity.fasta \
          $(find 03_transrate/${ID} -name contigs.csv) 03_filtered

  # 4. Chimera removal
  python $PDC/run_chimera_detection.py 03_filtered/${ID}.good_transcripts.short_name.fa \
          databases/alternifolia.fa $THREADS 03_chimera

  # 5. Corset clustering; longest transcript per cluster
  python $PDC/corset_wrapper.py 03_chimera/${ID}.filtered_transcripts.fa $Q1 $Q2 \
          $THREADS 03_corset salmon
  python $PDC/filter_corset_output.py 03_chimera/${ID}.filtered_transcripts.fa \
          03_corset/${ID}_salmon-clusters.txt 03_reps

  # 6. ORF prediction
  mkdir -p 04_transdecoder/${ID}
  python $PDC/transdecoder_wrapper.py 03_reps/${ID}.largest_cluster_transcripts.fa \
          $THREADS non-stranded 04_transdecoder/${ID}

  # 7. Sequence names to Genus_species_code@transcript
  NAME=$(awk -v id=$ID '$1 == id {print $2}' $REPO/data/taxon_names.txt)
  sed -i "s/^>${ID}@/>${NAME}@/" 04_transdecoder/${ID}/${ID}.cds.fa 04_transdecoder/${ID}/${ID}.pep.fa

  # 8. Redundancy reduction of coding sequences
  cd-hit-est -i 04_transdecoder/${ID}/${ID}.cds.fa -o 05_cdhit/${ID}.cds.cdhit \
             -c 0.99 -n 10 -r 0 -T $THREADS -M 0

  # 9. Completeness (BUSCO, eudicots_odb10)
  busco -i 02_trinity/${ID}.Trinity.fasta -m transcriptome -l eudicots_odb10 \
        -c $THREADS -o busco/${ID} --offline
done
