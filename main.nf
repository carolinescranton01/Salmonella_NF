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
            tuple(
                row.sample,
                file(row.fastq_1),
                file(row.fastq_2)
            )
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

    container 'quay.io/biocontainers/checkm-genome:1.2.2--pyhdfd78af_0'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}_R1.fq.gz"), path("${sample}_R2.fq.gz")

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
 * ASSEMBLY STEP (SHOVILL FIXED SINGLE-END)
 * -------------------------
 */

process ASSEMBLY {

    tag "$sample"

    cpus 4
    memory '8 GB'

    container 'quay.io/biocontainers/shovill:1.1.0--hdfd78af_1'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}.fa")

    script:
    """
    shovill \
        --R1 ${read1} \
        --R2 ${read2} \
        --outdir shovill_out \
        --cpus ${task.cpus}

    cp shovill_out/contigs.fa ${sample}.fa
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
    memory '8 GB'

    container 'staphb/checkm2:latest'

    input:
    tuple val(sample), path(assembly)

    output:
    path("checkm_out")

    script:
    """
    mkdir genomes
    cp ${assembly} genomes/

    checkm2 predict \
        -i genomes \
        -o checkm_out \
        --threads ${task.cpus}
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
 * SISTR ANNOTATION
 * -------------------------
 */

process SISTR {

    tag "$sample"

    cpus 2
    memory '4 GB'

    container 'staphb/sistr_cmd:latest'

    input:
    tuple val(sample), path(assembly)

    output:
    path("${sample}_sistr")

    script:
    """
    set -e

    echo "Checking assembly:"
    ls -lh ${assembly}
    head -n 5 ${assembly}

    rm -rf ${sample}_sistr
    mkdir ${sample}_sistr

    sistr \
        --input ${assembly} \
        --output Sample1_sistr \
        --qc \
        --alleles-output ${sample}_sistr/${sample}_alleles.json
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
    SISTR(asm_ch)
}
