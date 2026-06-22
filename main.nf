#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads      = params.reads ?: "*_{R1,R2}_001.fastq.gz"
params.outdir     = params.outdir ?: "results"
params.bakta_db   = params.bakta_db ?: null
params.checkm_db  = params.checkm_db ?: null

/*
Expected input files:

Sample1_R1_001.fastq.gz
Sample1_R2_001.fastq.gz
Sample2_R1_001.fastq.gz
Sample2_R2_001.fastq.gz
*/

Channel
.fromFilePairs(params.reads, checkIfExists: true)
.set { reads_ch }

//
// FASTP
//

process FASTP {

```
tag "$sample"

publishDir "${params.outdir}/fastp", mode: 'copy'

cpus 8
memory '16 GB'

input:
tuple val(sample), path(reads)

output:
tuple val(sample),
      path("${sample}_R1.trim.fastq.gz"),
      path("${sample}_R2.trim.fastq.gz")

script:
"""
fastp \
    -i ${reads[0]} \
    -I ${reads[1]} \
    -o ${sample}_R1.trim.fastq.gz \
    -O ${sample}_R2.trim.fastq.gz \
    --detect_adapter_for_pe \
    --thread ${task.cpus} \
    --html ${sample}.html \
    --json ${sample}.json
"""
```

}

//
// FASTQC
//

process FASTQC {

```
tag "$sample"

publishDir "${params.outdir}/fastqc", mode: 'copy'

cpus 4

input:
tuple val(sample), path(r1), path(r2)

output:
path("*_fastqc.zip")
path("*_fastqc.html")

script:
"""
fastqc -t ${task.cpus} ${r1} ${r2}
"""
```

}

//
// SPADES
//

process SPADES {

```
tag "$sample"

publishDir "${params.outdir}/assemblies", mode: 'copy'

cpus 16
memory '64 GB'

input:
tuple val(sample), path(r1), path(r2)

output:
tuple val(sample), path("${sample}.fasta")

script:
"""
spades.py \
    --isolate \
    --careful \
    -1 ${r1} \
    -2 ${r2} \
    -t ${task.cpus} \
    -m 60 \
    -o spades

cp spades/scaffolds.fasta ${sample}.fasta
"""
```

}

//
// QUAST
//

process QUAST {

```
tag "$sample"

publishDir "${params.outdir}/quast", mode: 'copy'

cpus 8

input:
tuple val(sample), path(fasta)

output:
path("${sample}_quast")

script:
"""
quast.py ${fasta} \
    -o ${sample}_quast \
    -t ${task.cpus}
"""
```

}

//
// CHECKM2
//

process CHECKM2 {

```
tag "$sample"

publishDir "${params.outdir}/checkm2", mode: 'copy'

cpus 16
memory '32 GB'

input:
tuple val(sample), path(fasta)

output:
path("${sample}_checkm2")

script:
"""
mkdir genomes
cp ${fasta} genomes/

checkm2 predict \
    --threads ${task.cpus} \
    --input genomes \
    --output-directory ${sample}_checkm2 \
    --database_path ${params.checkm_db}
"""
```

}

//
// BAKTA
//

process BAKTA {

```
tag "$sample"

publishDir "${params.outdir}/bakta", mode: 'copy'

cpus 16
memory '32 GB'

input:
tuple val(sample), path(fasta)

output:
tuple val(sample),
      path("${sample}_bakta/${sample}.gff3")

script:
"""
bakta \
    --db ${params.bakta_db} \
    --threads ${task.cpus} \
    --output ${sample}_bakta \
    --prefix ${sample} \
    ${fasta}
"""
```

}

//
// SISTR
//

process SISTR {

```
tag "$sample"

publishDir "${params.outdir}/sistr", mode: 'copy'

cpus 8

input:
tuple val(sample), path(fasta)

output:
path("${sample}_sistr.csv")

script:
"""
sistr \
    --threads ${task.cpus} \
    --csv ${sample}_sistr.csv \
    ${fasta}
"""
```

}

//
// PANAROO
//

process PANAROO {

```
publishDir "${params.outdir}/panaroo", mode: 'copy'

cpus 32
memory '128 GB'

input:
path(gffs)

output:
path("core_gene_alignment.aln")
path("gene_presence_absence.csv")

script:
"""
mkdir gffs
cp ${gffs} gffs/

panaroo \
    -i gffs/*.gff3 \
    -o panaroo_out \
    --clean-mode strict \
    --core_threshold 0.99 \
    -t ${task.cpus}

cp panaroo_out/core_gene_alignment.aln .
cp panaroo_out/gene_presence_absence.csv .
"""
```

}

//
// RAXML
//

process RAXML {

```
publishDir "${params.outdir}/raxml", mode: 'copy'

cpus 32
memory '64 GB'

input:
path(aln)

output:
path("RAxML*")
path("*.raxml.*")

script:
"""
raxml-ng \
    --all \
    --msa ${aln} \
    --model GTR+G \
    --threads ${task.cpus} \
    --bs-trees 100 \
    --prefix salmonella_tree
"""
```

}

workflow {

```
trimmed = FASTP(reads_ch)

FASTQC(trimmed)

assemblies = SPADES(trimmed)

QUAST(assemblies)

CHECKM2(assemblies)

SISTR(assemblies)

annotations = BAKTA(assemblies)

gff_files = annotations.map { sample, gff -> gff }.collect()

panaroo = PANAROO(gff_files)

RAXML(panaroo.out[0])
```

}
