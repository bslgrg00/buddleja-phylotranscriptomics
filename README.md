# Phylotranscriptomics of Asian *Buddleja*

Scripts for the analyses in:

> Gurung, B., J. Liu, F. R. Rahaingoson, F. Yang, G. Chen, and J. Ge. Phylotranscriptomics reveals multiple sources of gene tree discordance in the rapid radiation of Asian *Buddleja*. *xxx* (submitted).

- Data (assemblies, alignments, trees, and results): Zenodo, https://doi.org/10.5281/zenodo.24182112
- Raw reads: NGDC Genome Sequence Archive, CRA052186 (BioProject PRJCA077110)

## Repository contents

| File | Analysis |
|---|---|
| `scripts/01_assembly.sh` | Read filtering, Trinity assembly, transcript filtering, TransDecoder, CD-HIT, BUSCO |
| `scripts/02_orthology.sh` | Homolog clustering and trees, tree cleaning, MO orthologs, PHI test, supermatrices (M10, M20, M30) |
| `scripts/03_phylogeny.sh` | Gene tree sets, IQ-TREE, RAxML-NG, ASTRAL-III, ASTRAL-Pro3, collapsed gene trees |
| `scripts/03_RF_distances.R` | Robinson-Foulds distances among the 12 species trees |
| `scripts/04_concordance.sh` | Phyparts; gene and site concordance factors |
| `scripts/05_MSCquartets.R` | Quartet tests (T1, T3) and NANUQ distances |
| `scripts/05_anomaly_zone.R` | Anomaly-zone internode pairs |
| `scripts/05_phytop.sh` | Phytop ILS and hybridization indices |
| `scripts/06_gene_flow.sh` | PhyloNet and HyDe |
| `scripts/07_relative_importance.sh` | ILS, gene flow, and estimation error predictors (Cai et al., 2021) |
| `scripts/07_relaimpo.R` | Relative importance (LMG, relaimpo) |
| `scripts/08_dating.sh` | treePL, SortaDate, and BEAST 2 |
| `helpers/compare_astral_trees.py` | Compares ASTRAL trees from collapsed and uncollapsed gene trees |
| `helpers/restore_bootstrap_support.py` | Restores bootstrap support to TreeShrink-trimmed trees |
| `helpers/mask_tips_by_taxonID_transcripts.py` | Yang and Smith masking script adapted for CDS |
| `data/taxon_names.txt` | Accession code and taxon name |
| `data/taxon_file.txt` | Ingroup (IN) and outgroup (OUT) taxa for `prune_paralogs_MO.py` |

## Usage

1. Install the software listed below and set the `path/to/` variables at the top of each script.
2. Place paired reads as `reads/<code>_R1.fq.gz` and `reads/<code>_R2.fq.gz`, and the *Arabidopsis thaliana* and *Buddleja alternifolia* proteomes as `databases/Atha.fa` and `databases/alternifolia.fa`.
3. Run the scripts in numerical order from one working directory, for example:

```bash
bash path/to/buddleja-phylotranscriptomics/scripts/01_assembly.sh
Rscript path/to/buddleja-phylotranscriptomics/scripts/05_MSCquartets.R
```

Each script header lists its output files and the Zenodo folder that holds them.

## Taxon codes

Tip names follow `Genus_species_code` (e.g., `Buddleja_colvilei_AB1`). AB = Asian *Buddleja*; OWAS, SNOD, UMUL, and WUZV = outgroups. The two non-Asian *Buddleja* accessions are labeled OUT1 and OUT2 in all files and correspond to NAB1 (*B. americana*) and NAB2 (*B. aromatica*) in the article.

## Software

| Step | Software |
|---|---|
| Assembly | fastp v0.21.0, Trinity v2.15.1, TransRate v1.0.3, Corset v1.09, Salmon, TransDecoder v5.7.1, CD-HIT v4.8.1, BUSCO v6.0.0 |
| Orthology | BLAST+, MCL v22-282, MAFFT v7.525, phyx v1.3.2, RAxML v8.2.12, TreeShrink v1.3.9, PhiPack |
| Phylogeny | IQ-TREE v2.1.4, RAxML-NG, ASTRAL-III v5.7.8, ASTRAL-Pro3 v1.24.3.7, Newick Utilities |
| Discordance and ILS | Phyparts, MSCquartets v3.3, Phytop |
| Gene flow | PhyloNet v3.8.2, HyDe v1.0.2, SplitsTree v4.19.2 |
| Relative importance | Phybase v1.1, Seq-Gen v1.3.5, relaimpo |
| Dating | treePL v1.0, SortaDate, BEAST v2.6.3 (BEAUti, LogCombiner, TreeAnnotator) |
| R packages | ape, phangorn, ggplot2, reshape2, MSCquartets, relaimpo |

Third-party scripts called by the pipeline:

- Yang and Smith pipeline: https://bitbucket.org/yanglab/phylogenomic_dataset_construction
- PhiPack: https://github.com/julianzaugg/phipack
- Phyparts: https://bitbucket.org/blackrim/phyparts
- phypartspiecharts.py: https://github.com/mossmatters/phyloscripts
- Phytop: https://github.com/zhangrengang/phytop
- HyDe: https://github.com/pblischak/HyDe
- Cai et al. (2021): https://github.com/lmcai/Coalescent_simulation_and_gene_flow_detection
- SortaDate: https://github.com/FePhyFoFum/SortaDate

## License

Code is released under the MIT License. Data on Zenodo are released under CC BY 4.0.

## Contact

Bishal Gurung (bsl.grg1996@gmail.com)
