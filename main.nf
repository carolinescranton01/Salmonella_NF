nextflow.enable.dsl = 2

params.samples = params.samples
params.outdir  = params.outdir ?: "results"

// ----------------------
// INPUT CHANNEL
// ----------------------
Channel
    .fromPath(params.samples)
    .splitCsv(header:true, sep:'\t')
    .map { row ->
        tuple(row.sample, file(row.fastq))
    }
    .set { READS }


// ----------------------
// QC (fastp)
// ----------------------
process QC {

    tag "$sample"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.fq.gz"), emit: reads_qc

    script:
    """
    fastp \
        -i ${reads} \
        -o ${sample}.fq.gz \
        -q 20 -l 50 \
        --thread 4
    """
}


// ----------------------
// ASSEMBLY (Shovill)
// ----------------------
process ASSEMBLY {

    tag "$sample"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.fa"), emit: contigs

    script:
    """
    shovill \
        --R1 ${reads} \
        --outdir . \
        --cpus 4

    mv contigs.fa ${sample}.fa
    """
}


// ----------------------
// CHECKM2
// ----------------------
process CHECKM {

    tag "$sample"

    input:
    tuple val(sample), path(contigs)

    output:
    path "checkm_out/**"

    script:
    """
    checkm lineage_wf \
        -x fa \
        ${contigs} \
        checkm_out
    """
}


// ----------------------
// BAKTA
// ----------------------
process BAKTA {

    tag "$sample"

    input:
    tuple val(sample), path(contigs)

    output:
    path "bakta_out/**"

    script:
    """
    bakta \
        --output bakta_out \
        ${contigs}
    """
}


// ----------------------
// WORKFLOW
// ----------------------
workflow {

    READS
        | QC
        | ASSEMBLY
        | CHECKM
        | BAKTA
}
