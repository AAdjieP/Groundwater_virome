# Metagenomics processing and Metagenomics-assembled genomes (MAGs) workflow and scripts
- Metagenome processing: quality control (BBDuk) and reporting (fastp)
- Metagenome assembly: MetaSpades, MEGAHIT
- Metagenome binning: Binsanity, abawaca, maxbin2, CONCOCT and metabat2 
- MAG dereplication: drep
- MAG read-mapping: CoverM
- MAG Taxonomy analysis: GTDBTK
- MAG annotation: DRAM
- Long-reads processing: guppy_basecaller, assembly (hybrid)

## Metagenome processing
**_Quality control_**
```
# Read more about BBDuk: https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/
bbduk.sh -Xmx1g \
        in1=READS_FORWARD \
        in2=READS_REVERSE \
        out1=OUTPUT_1 \
        out2=OUTPUT_2 \
        statscolumns=5 \
        stats=SAMPLE_STATS.txt \
        ref=/PATH/adapters.fa \
        threads=INTEGER \
        ktrim=r \
```
**_Quality control reporting_**
```
# Read more about fastp: https://github.com/OpenGene/fastp
fastp -i READS-SAMPLE_R1.fastq.gz -I READS-SAMPLE_R2.fastq.gz -R SAMPLE_fastp_report -j  SAMPLE_fastp.json -h SAMPLE_fastp.html -w INTEGER
```

## Assemply process
**_MetaSpades_**
```
# Read more about Spades: [https://github.com/AnantharamanLab/VIBRANT](https://github.com/ablab/spades)
# Spade paper: [https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpbi.102)
module load SPAdes/3.15.2
spades.py --meta -o ${sample} -1 ${sample}_R1.fastq.gz -2 ${sample}_R2.fastq.gz -t INTEGER -k 21,33,55,77 -m 190
```

**_MEGAHIT_**
```
# Read more about MEGAHIT: https://github.com/voutcn/megahit
# MEGAHIT paper: https://academic.oup.com/bioinformatics/article/31/10/1674/177884
module load MEGAHIT/1.2.9-Python-2.7.18
megahit -1 SAMPLE_1.fastq.gz -2 SAMPLE}_2.fastq.gz -o OUTPUT -t INTEGER
```

## Binning process
**_MetaSpades_**
```
# Read more about MetaWrap: https://github.com/bxlab/metaWRAP
# MetaWrap paper: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0541-1
# 40 markers
metawrap binning -o /OUTPUT/ -t INT --universal -a /INPUT/assembly.fasta --maxbin2 --metabat2 --concoct /READS_DIR/reads_1.fastq /READS_DIR/reads_2.fastq
# 107 markers
metawrap binning -o /OUTPUT/ -t INT -a /INPUT/assembly.fasta --maxbin2 /READS_DIR/reads_1.fastq /READS_DIR/reads_2.fastq
```

**_BinSanity_**
```
# Read more about BinSanity: https://github.com/edgraham/BinSanity; https://github.com/edgraham/BinSanity/wiki
# BinSanity paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC5345454/
# Binsanity-profile uses featureCounts to produce the coverage profiles requires in Binsanity
Binsanity-profile -i /INPUT/*.fasta -s /metaWRAP_out/work_files/ -T 24 -o /cov_profiles/ --ids /PATH/SAMPLE.txt -c /cov_profiles/binsanity

# Binsanity-wf runs Binsanity and Binsanity-refine sequentially to optimize cluster results
Binsanity-wf -f /metaWRAP_out/metawrap_bins/ -l /PATH/*.fasta -c /cov_profiles/binsanity.cov.x100.lognorm -o /final_bins --threads INTEGER
```

**_Abawaca_**
```
# Read more about Abawaca: https://github.com/CK7/abawaca
# Running abawaca=1.0.7
# Generate ESOM (Emergent Self-Organizing Map) files that are data inputs used by ABAWACA to perform binning based on machine learning clustering
perl /global/apps/metagenomics/prepare_esom_files.pl \
  -m INTEGER \  # minimal contigs length We used 5000 and 10000
  --coverage \
  -sa $(ls /work/.../metawrap_binning/work_files/*.sorted.bam) \
  /work/.../abawaca_binning/ \
  /work/.../scaffolds_3000.fasta
# Run abawaca  
abawaca /PATH/esom.names /PATH/esom.lrn /PATH/*.fasta /OUTPUT/
```

**_Bins Refinements_**
```
# Refining ABAWACA + BINSANITY Bins
metawrap bin_refinement -o /OUTPUT -t  INTEGER -c 50 -x 10 -A /INPUT_abawaca_BINS_1 -B /INPUT_abawaca_BINS_2 -C /INPUT_binsanity_BINS_1
# FInal refinement
metawrap bin_refinement -o /OUTPUT -t  INTEGER -c 50 -x 10 -A /PATH/metabat2_bins -B /PATH/concoct_bins -C /PATH/metawrap_50_10_bins
```

## MAG DEREPLICATION 
```
# Read more about drep: https://github.com/MrOlm/drep
# drep paper: https://academic.oup.com/ismej/article/11/12/2864/7537826?login=false
dRep dereplicate --debug -p INTEGER -pa 0.90 -sa 0.99 -comp 50 -con 10 /PATH/OUTPUT -g /OUTPUT/*.fasta
```

## MAG Read-mapping
```
# Read more about CoverM: https://github.com/wwood/CoverM
# CoverM paper: https://doi.org/10.48550/arXiv.2501.11217
# CoverM version 0.4.0
coverm genome -1 READS-SAMPLE_R1.fastq.gz -2 READ-SAMPLE_R2.fastq.gz --genome-fasta-directory /PATH/DIR --genome-fasta-extension fasta --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 --min-covered-fraction 0.25 --methods trimmed_mean --bam-file-cache-directory BAMS_OUTPUT --discard-unmapped -t INTEGER > SAMPLE.txt
```

## Taxonomy assignment analysis
```
# Read more about GTDBTK: https://github.com/Ecogenomics/GTDBTk; https://gtdb.ecogenomic.org
# GTDBTK paper: https://doi.org/10.1093/bioinformatics/btac672
module load GTDB-Tk/1.4.1
gtdbtk classify_wf --genome_dir /PATH/DIR/ --extension fasta --cpus INTEGER --out_dir /PATH/GTDBTK_OUTPUT/
```

## MAG annotation
```
# Read more about DRAM-v: https://github.com/WrightonLabCSU/DRAM
# DRAM-v paper: https://academic.oup.com/nar/article/48/16/8883/5884738
module load DRAM/1.3.6
DRAM.py annotate -i MAG_INPUT.fasta -o OUTOUT_ANNOTATION --min_contig_size 1000 --threads INTEGER
DRAM.py distill -i /PATH/annotations.tsv -o DISTILL_OUTPUT
```

## Long-reads metagenome processing
**_Quality control_**
```
module load gcccuda/2020b
module load fosscuda/2020b
module load GCC/10.2.0
module load CUDA/11.1.1
guppy_basecaller -i fast5/ -s SAMPLE -c dna_r9.4.1_450bps_sup.cfg --compress_fastq  --recursive  -x "cuda:0" --trim_primers --trim_adapters  --do_read_splitting
guppy_basecaller -i /PATH/fast5 -s /PATH/DIR/
```

**_Hybrid assembly_**
```
module load SPAdes/3.15.2
spades.py --meta -o OUTPUT_Nanopore -1 READS-SAMPLE_R1.fastq.gz -2 READS-SAMPLE_R2.fastq.gz -t INTEGER -m 380 --nanopore /PATH/DIR/
```



