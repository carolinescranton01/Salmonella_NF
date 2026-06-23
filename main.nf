nextflow.enable.dsl = 2

/*
 * -------------------------
 * PARAMETERS
 * -------------------------
 */

params.samples = null
params.outdir  = "results"

/*
 * -------------------------
 * INPUT CHANNEL (SINGLE-END FIX)
 * -------------------------
 */

samples_ch =
    Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(row.sample, file(row.fastq))
        }

/*
 * -------------------------
 * QC STEP (fastp)
 * -------------------------
 */

process QC {

    tag "$sample"

    cpus 4
    memory '8 GB'

    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.fq.gz")

    script:
    """
    fastp \
        -i ${reads} \
        -o ${sample}.fq.gz \
        -q 20 -l 50 \
        --thread ${task.cpus}
    """
}

/*
 * -------------------------
 * ASSEMBLY STEP (SHOVILL FIXED SINGLE-END)
 * -------------------------
 */

process ASSEMBLY {

    tag "$sample"

    cpus 4
    memory '8 GB'

    container 'quay.io/biocontainers/shovill:1.1.0--hdfd78af_1'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.fa")

    script:
    """
    shovill \
        --R1 ${reads} \
        --outdir . \
        --cpus ${task.cpus}

    mv contigs.fa ${sample}.fa
    """
}

/*
 * -------------------------
 * CHECKM STEP
 * -------------------------
 */

process CHECKM {

    tag "$sample"

    cpus 4
    memory '16 GB'

    container 'quay.io/biocontainers/checkm2:1.0.1--pyhdfd78af_0'

    input:
    tuple val(sample), path(contigs)

    output:
    tuple val(sample), path("checkm.txt")

    script:
    """
    checkm lineage_wf \
        -x fa \
        . checkm_out \
        --threads ${task.cpus}

    touch checkm.txt
    """
}

/*
 * -------------------------
 * BAKTA ANNOTATION
 * -------------------------
 */

process BAKTA {

    tag "$sample"

    cpus 4
    memory '16 GB'

    container 'quay.io/biocontainers/bakta:1.9.3--pyhdfd78af_0'

    input:
    tuple val(sample), path(contigs)

    output:
    tuple val(sample), path("${sample}.gbk")

    script:
    """
    bakta \
        --db light \
        --output . \
        --prefix ${sample} \
        ${contigs}
    """
}

/*
 * -------------------------
 * WORKFLOW
 * -------------------------
 */

workflow {

    QC(samples_ch)
        .set { qc_ch }

    ASSEMBLY(qc_ch)
        .set { asm_ch }

    CHECKM(asm_ch)
    BAKTA(asm_ch)
}
