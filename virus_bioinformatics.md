# Virus bioinformatic workflow and scripts
- Virus identification: DeepVirFinder, VIBRANT, geNomad, and VirSorter2
- Virus clustering (to get vOTU): MMseq2
- Virus genome quality assessment: CheckV
- Virus read-mapping: CheckM
- Virus taxonomy: vConTACT3
- Virus annotation and AMG identification: DRAM-v
- Virus microdiversity assessment: MetaPop
- Virus host prediction: iPhop
- CRISPR-Cas spacers identificationa and matches: minced and blastn

## Virus identification:
**DeepVirFinder virus identification**
```
# Read more about DeepVirfinder: https://github.com/jessieren/DeepVirFinder
# DeepVirFinder paper: https://link.springer.com/article/10.1007/s40484-019-0187-4

module load DeepVirFinder/1.0
dvf.py -i /PATH/INPUT.fasta -o OUTPUT_DIR -l 1000 -c 28
```
**VIBRANT virus identification**
```
# Read more about VIBRANT: https://github.com/AnantharamanLab/VIBRANT
# VIBRANT paper: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0

module load "VIBRANT/1.2.1"
VIBRANT_run.py -i /PATH/INPUT.fasta -f nucl -t 28 -virome -d /PATH/VIBRANTDB/
```

**geNomad virus identification**
```
# Read more about geNomad: https://github.com/apcamargo/genomad
# geNomad paper: https://www.nature.com/articles/s41587-023-01953-y
genomad end-to-end --cleanup --splits 48 --threads 20 --min-virus-marker-enrichment 1 --min-virus-hallmarks 1 /PATH/INPUT.fasta genomad_output /PATH/genomad_db
```

**VirSrter2 virus identification**
```
# Read more about VirSorter2: https://github.com/jiarong/VirSorter2
# VirSrter2 SOP: https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3
# VirSrter2 paper: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00990-y
VirSorter2-2.2.3.sif run -w DIR -i /PATH/INPUT.fasta --include-groups dsDNAphage,ssDNA --keep-original-seq -j 48 --min-score 0.5 --min-length 5000 all
```

## Virus clustering:
Virus clustering at 95% ID and 80% coverage using MMseq2
```
# Read more about MMseq2: https://github.com/soedinglab/MMseqs2
# MMseq2 paper: https://www.nature.com/articles/nbt.3988
mmseqs easy-cluster /PATH/INPUT.fasta /PATH/output --min-seq-id 0.95 -c 0.8
```

## Virus genome quality assessment:
```
# Read more about CheckV: https://bitbucket.org/berkeleylab/checkv/src/master/
# CheckV paper: https://www.nature.com/articles/s41587-020-00774-7
module load CheckV/2021.02.03
checkv end_to_end -t 48 /PATH/INPUT.fasta OUTPUT
```

## Virus read-mapping:
vOTU read mapping performed with CoverM version 0.4.0
```
# Read more about CoverM: https://github.com/wwood/CoverM
# CoverM paper: https://doi.org/10.48550/arXiv.2501.11217
coverm contig --coupled /PATH/INPUT_R1.fastq.gz /PATH/INPUT_R2.fastq.gz --reference /PATH/INPUT.fasta --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 --min-covered-fraction 0.7 -m trimmed_mean --bam-file-cache-directory BAM_FILES --discard-unmapped -t 28 > OUTPUT_mapping_SAMPLE.txt
```

## Virus taxonomy analysis: 
```
# Read more about vConTACT3: https://bitbucket.org/MAVERICLab/vcontact3/src/master/
vcontact3 run --nucleotide /PATH/INPUT.fasta --db-domain prokaryotes --db-version 220 --output vc3_output_HT --exports cytoscape --db-path /PATH/vcontact3_dbs
```

## Virus annotation and AMG identification:
```
# Read more about DRAM-v: https://github.com/WrightonLabCSU/DRAM
# DRAM-v paper: https://academic.oup.com/nar/article/48/16/8883/5884738
#first, prepare the input data from VirSorter2 --prep-for-dramv
VirSorter2-2.2.3.sif run --seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv -i checkv/combined.fna -w DIR --include-groups dsDNAphage,ssDNA --min-length 5000 --min-score 0.5 -j 28 all

#DRAM-v annotation
module load DRAM
DRAM-v.py annotate -i /PATH/INPUT.fasta -v /PATH/viral-affi-contigs-for-dramv.tab -o HT_final --min_contig_size 1000 --threads 48

#DRAM-v distill
DRAM-v.py distill -i /PATH/annotations.tsv -o distill
```

## Virus microdiversity assessment:
```
# Read more about MetaPop: https://github.com/metaGmetapop/metapop 
# MetaPop paper: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-022-01231-0
module load metapop
#first, you need to split all the votus: FASTAS/
#second, prepared bam files, the mapping files of all the contigs to the reads: BAM_FILES/
#third, prepared a list of reads count (read_counts.txt) tab-seperated: sample reads
metapop --input_samples /PATH/BAM_FILES/ --reference /PATH/FASTAS/ --norm /PATH/read_counts.txt --threads 48 --min_cov 70
```
Using global_contig_microdiversity.tsv in "10.Microdiversity" directory, the final microdiversity value (average π) for each sample was determined by averaging π values from 100 randomly selected viral populations across 1,000 subsamplings

## Virus host prediction analysis:
```
# Read more about iPhop: https://bitbucket.org/srouxjgi/iphop/src/main/
# iPhop paper: https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002083
#first, prepare the gtdbtk analysis: MAG_gtdbtk/
module load GTDB-Tk/1.4.1
gtdbtk de_novo_wf --genome_dir /PATH/MAG_DIR/ --bacteria --outgroup_taxon p__Cyanobacteria --out_dir MAG_bac_gtdbtk/ --cpus 48 --force --extension fasta
gtdbtk de_novo_wf --genome_dir /PATH/MAG_DIR/ --archaea --outgroup_taxon p__Altarchaeota --out_dir MAG_arc_gtdbtk/ --cpus 48 --force --extension fasta

#second, add database (Metagenome-assembled genomes) as a symbolic links to the original iPhop database (/PATH/iPHoP/Sept_2021_pub_rw)
iphop add_to_db --fna_dir /PATH/MAG_DIR/ --gtdb_dir /PATH/MAG_gtdbtk/ --out_dir iphop_add_db --db_dir /fs/ess/PAS1117/apratama/iphop-db/Sept_2021_pub_rw/ -t 48

#third, run iPhop with default database
iphop predict --fa_file /PATH/INPUT.fasta --out_dir OUTPUT_DIR --db_dir /PATH/Sept_2021_pub_rw/ --num_threads 48

#forth, run iPhop with the added database
iphop predict --fa_file /PATH/INPUT.fasta --out_dir OUTPUT_DIR --db_dir /PATH/iphop_add_db/ --num_threads 48
```

## Identification of CRISPR-Cas spacers:
This analysis was done to identify pro-spacers from MAGs. We then assess whether the proto-spacer matches between viruses and both CPR/DPANN and non-CPR/non-DPANN bacterial and archaeal MAGs
```
# Read more about minced: https://github.com/ctSkennerton/minced
module load minced
minced -minNR 3 -spacers /PATH/INPUT.fasta OUTPUT_CRISPR-Cas_spacers

# Spacer blastn
module load blast

#first: makeblastdb -in FINAL_vOTU.fasta -dbtype nucl -out FINAL_vOTU_db
#see more on: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.options_common_to_all_blast/
blastn \
  -query MAG_spacers.fa \
  -db /PATH/FINAL_vOTU_db \
  -task blastn-short \
  -word_size 7 \
  -evalue 1e-5 \
  -reward 1 -penalty -1 \
  -ungapped \
  -outfmt "6 qseqid sseqid pident length mismatch qlen qstart qend sstart send evalue bitscore qcovs" \
  -dust no \
  -num_threads 4 \
  -out BLAST_OUTPUT.txt


#Filter for ≤1 Mismatch
awk '$3 >= 95 && $13 >= 90' BLAST_RESULTS.txt > filtered_BLAST_RESULTS.txt

#Filter for 0 Mismatch
awk '$5 == 0 && $13 >= 90' BLAST_RESULTS.txt > filtered_BLAST_RESULTS.txt

#Filter for full-length matches only (100% query coverage)
awk '$5 == 0 && $13 == 100' BLAST_RESULTS.txt > filtered_BLAST_RESULTS.txt
```


