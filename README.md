# Nextflow pipeline for APGAP Salmonella analysis - Cooper Lab

### Input
Paired-end Illumina reads:
*_R1_001.fastq.gz
*_R2_001.fastq.gz

### Run
nextflow run main.nf \
  --reads "*_{R1,R2}_001.fastq.gz" \
  --bakta_db /db/bakta \
  --checkm_db /db/checkm2 \
  -profile docker

