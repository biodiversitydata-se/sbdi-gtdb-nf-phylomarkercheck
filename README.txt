# SBDI Sativa curated 16S GTDB database

## General information

Author: SBDI molecular data team
Contact e-mail: daniel.lundin@lnu.se, anders.andersson@scilifelab.se
DOI: 10.17044/scilifelab.14869077
License: CC BY 4.0
Version: R08-RS214-1
Categories: Bacteriology (310701), Microbial ecology (310703), Microbial genetics (310704), Medical bacteriology (320701), Medical microbiology not elsewhere classified (320799)
Item type: Dataset
Keywords: 16S rRNA, GTDB, DADA2, SBDI, Ampliseq
Funding: Curation of this data was funded by the Swedish Research Council (VR), grant number 2019-00242.

This readme file was last updated: 2024-04-24

Please cite as: Swedish Biodiversity Infrastructure (SBDI; 2021). SBDI Sativa curated 16S GTDB database. https://doi.org/10.17044/scilifelab.14869077

## Dataset description

The data in this [repository](https://doi.org/10.17044/scilifelab.14869077) is the result of vetting 16S sequences from the GTDB database release R08RS214 (r214) (https://gtdb.ecogenomic.org/; Parks et al. 2018) with the Sativa program (Kozlov et al. 2016)
using the [sbdi-phylomarkercheck](https://github.com/biodiversitydata-se/sbdi-phylomarkercheck) Nextflow pipeline.

Using Sativa [Kozlov et al. 2016], 16S sequences from GTDB were checked so that their phylogenetic signal is consistent with their taxonomy.

Before calling Sativa, sequences longer than 2000 nucleotides or containing Ns were removed, and the reverse complement of each is calculated.
Subsequently, sequences were aligned with HMMER [Eddy 2011] using the Barrnap [https://github.com/tseemann/barrnap] archaeal and bacterial 16S profiles respectively, 
and sequences containing more than 10% gaps were removed.
The remaining sequences were analyzed with Sativa, and sequences that were not phylogenetically consistent with their taxonomy were removed.

Files for the DADA2 (Callahan et al. 2016) methods `assignTaxonomy` and `addSpecies` are available, in three different versions each. 
The `assignTaxonomy` files contain taxonomy for domain, phylum, class, order, family, genus and species. 
(Note that it has been proposed that species assignment for short 16S sequences require 100% identity (Edgar 2018), so use species assignments from `assignTaxonomy` with caution.) 
The versions differ in the maximum number of genomes that we included per species: 1, 5 or 20, indicated by "1genome", "5genomes" and "20genomes" in the file names respectively. 
Using the version with 20 genomes per species should increase the chances to identify an exactly matching sequence by the `addSpecies` algorithm, while using a file with many genomes 
per species could potentially give biases in the taxonomic annotations at higher levels by `assignTaxonomy`.
Our recommendation is hence to use the "1genome" files for `assignTaxonomy` and "20genomes" for `addSpecies`.

All files are gzipped fasta files with 16S sequences, the assignTaxonomy associated with taxonomy hierarchies from domain to species whereas the `addSpecies` file have sequence identities and species names.
There is also a fasta files with the original GTDB sequence names: sbdi-gtdb-sativa.r08rs214.20genomes.fna.gz.

Taxonomical annotation of 16S amplicons using this data is available as an optional argument to the nf-core/ampliseq Nextflow workflow from version 2.1: `--dada_ref_taxonomy sbdi-gtdb` 
(https://nf-co.re/ampliseq; Straub et al. 2020).

In addition to the fasta files, the workflow outputs phylogenetic trees by optimizing branch-lengths of the original phylogenomic GTDB trees based on a 16S sequence alignment.
As not all species in GTDB will have correct 16S sequences, the GTDB trees are first subset to contain only species for which the species representative genome has a correct 16S sequence.
Subsequently, branch lengths for the tree are optimized based on the original alignment of 16S sequences using IQTREE [Nguyen et al. 2015] with a GTR+F+I+G4 model.
The alignment files end with `.alnfna`, the taxonomy files with `.taxonomy.tsv` and the tree files (newick-formatted) end with `.brlenopt.newick`.
They will be made available in nf-core/ampliseq for phylogenetic placement.

## Version history

* v8 (2024-06-25): Update to GTDB R09-RS220 from R08-RS214.

* v7 (2024-04-24): Replace manual procedure with Nextflow pipeline. Update to GTDB R08-RS214 from R07-RS207.

* v6 (2022-10-07): Add missing fasta file with original GTDB names.

* v5 (2022-09-02): Update README (this document)

* v4 (2022-08-31): Update to GTDB R07-RS207 from R06-RS202

## Citations

Cock, Peter J. A., Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg, et al. 2009. “Biopython: Freely Available Python Tools for Computational Molecular Biology and Bioinformatics.” Bioinformatics 25 (11): 1422–23. https://doi.org/10.1093/bioinformatics/btp163.

Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

Eddy, Sean R. 2011. “Accelerated Profile HMM Searches.” PLoS Comput Biol 7 (10): e1002195. https://doi.org/10.1371/journal.pcbi.1002195.

Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh. 2015. “IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies.” Molecular Biology and Evolution 32 (1): 268–74. https://doi.org/10.1093/molbev/msu300.

Shen, Wei, Shuai Le, Yan Li, and Fuquan Hu. 2016. “SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation.” PLOS ONE 11 (10): e0163962. https://doi.org/10.1371/journal.pone.0163962.

Rice, P., I. Longden, and A. Bleasby. 2000. “EMBOSS: The European Molecular Biology Open Software Suite.” Trends in Genetics: TIG 16 (6): 276–77. https://doi.org/10.1016/s0168-9525(00)02024-2.

Straub, Daniel, Nia Blackwell, Adrian Langarica-Fuentes, Alexander Peltzer, Sven Nahnsen, and Sara Kleindienst. 2020. “Interpretations of Environmental Microbial Community Studies Are Biased by the Selected 16S RRNA (Gene) Amplicon Sequencing Pipeline.” Frontiers in Microbiology 11. https://doi.org/10.3389/fmicb.2020.550420.

## Acknowledgements

The computations were enabled by resources in project [NAISS 2023/22-601, SNIC 2022/22-500 and SNIC 2021/22-263] provided by the National Academic Infrastructure for Supercomputing in Sweden (NAISS) at UPPMAX, funded by the Swedish Research Council through grant agreement no. 2022-06725.
Computations were also enabled by resources provided by Dr. Maria Vila-Costa, Institute of Environmental Assessment and Water Research (IDAEA-CSIC), Barcelona.
