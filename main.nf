nextflow.enable.dsl=2

/*
 * -------------------------
 * PARAMETERS
 * -------------------------
 */

params.samples = null
params.outdir = "results"

/*
 * -------------------------
 * INPUT CHANNEL
 * -------------------------
 */

Channel.fromPath(params.samples)
    .splitCsv(header:true, sep:'\t')
    .map { row ->
        tuple(
            row.sample,
            file(row.fastq_1),
            file(row.fastq_2)
        )
    }
    .set { samples_ch }

/*
 * -------------------------
 * QC (fastp)
 * -------------------------
 */

process QC {

    tag "$sample"

    cpus 4
    memory '8 GB'

    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample),
          path("${sample}_R1.fq.gz"),
          path("${sample}_R2.fq.gz")

    script:
    """
    fastp \
        -i ${read1} \
        -I ${read2} \
        -o ${sample}_R1.fq.gz \
        -O ${sample}_R2.fq.gz \
        -q 20 -l 50 \
        --thread ${task.cpus}
    """
}

/*
 * -------------------------
 * ASSEMBLY (shovill)
 * -------------------------
 */

process ASSEMBLY {

    tag "$sample"

    cpus 4
    memory '8 GB'

    container 'staphb/shovill:1.1.0'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}.fa")

    script:
    """
    shovill \
        --R1 ${read1} \
        --R2 ${read2} \
        --outdir ${sample}_asm \
        --cpus ${task.cpus}

    cp ${sample}_asm/contigs.fa ${sample}.fa
    """
}

/*
 * -------------------------
 * CHECKM2
 * -------------------------
 */

process CHECKM2 {

    tag "$sample"

    cpus 4
    memory '16 GB'

    container 'staphb/checkm2:latest'

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("checkm2_out")

    script:
    """
    mkdir -p checkm2_db

    checkm2 database \
        --download \
        --path checkm2_db

    checkm2 predict \
        -i ${fasta} \
        -o checkm2_out \
        --database_path checkm2_db/CheckM2_database \
        -t ${task.cpus}
    """
}

/*
 * -------------------------
 * BAKTA
 * -------------------------
 */

process BAKTA {

    tag "$sample"

    cpus 4
    memory '16 GB'

    container 'staphb/bakta:latest'

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("bakta_out")

    script:
    """
    bakta \
        --input ${fasta} \
        --output bakta_out \
        --prefix ${sample} \
        --threads ${task.cpus}
    """
}

/*
 * -------------------------
 * SISTR
 * -------------------------
 */

process SISTR {

    tag "$sample"

    cpus 2
    memory '4 GB'

    container 'staphb/sistr_cmd:latest'

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("sistr_out")

    script:
    """
    mkdir sistr_out

    sistr \
        --input ${fasta} \
        --output sistr_out/${sample}.tsv \
        --alleles-output sistr_out/${sample}_alleles.json \
        --qc
    """
}

/*
 * -------------------------
 * WORKFLOW
 * -------------------------
 */

workflow {

    qc_ch = QC(samples_ch)

    asm_ch = ASSEMBLY(qc_ch)

    CHECKM2(asm_ch)
    BAKTA(asm_ch)
    SISTR(asm_ch)
}
