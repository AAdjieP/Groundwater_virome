This directory contains workflows for (_i_) metagenomics processing and Metagenomic assembled genomes (MAGs) and (_ii_) virus bioinformatics

**Metagenomics processing and Metagenomics-assembled genomes (MAGs) workflow and scripts**
- Metagenome processing: quality control (BBDuk) and reporting (fastp)
- Metagenome assembly: MetaSpades, MEGAHIT
- Metagenome binning: Binsanity, abawaca, maxbin2, CONCOCT and metabat2
- MAG dereplication: drep
- MAG read-mapping: CoverM
- MAG Taxonomy analysis: GTDBTK
- MAG annotation: DRAM
- Long-reads processing: guppy_basecaller, assembly (hybrid)

**Virus bioinformatic workflow and scripts** (Virus_bioinformatics.md)
- Virus identification: DeepVirFinder, VIBRANT, geNomad, and VirSorter2
- Virus clustering (to get vOTU): MMseq2
- Virus genome quality assessment: CheckV
- Virus read-mapping: CheckM
- Virus taxonomy: vConTACT3
- Virus annotation and AMG identification: DRAM-v
- Virus microdiversity assessment: MetaPop
- Virus host prediction: iPhop
- CRISPR-Cas spacers identificationa and matches: minced and blastn
