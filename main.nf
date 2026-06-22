#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { FASTP }     from './modules/fastp'
include { FASTQC }    from './modules/fastqc'
include { SPADES }    from './modules/spades'
include { QUAST }     from './modules/quast'
include { CHECKM2 }   from './modules/checkm2'
include { SISTR }     from './modules/sistr'
include { BAKTA }     from './modules/bakta'
include { PANAROO }   from './modules/panaroo'
include { RAXML }     from './modules/raxml'

params.reads  = null
params.outdir = "results"

workflow {

    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    trimmed_ch = FASTP(reads_ch)

    FASTQC(trimmed_ch)

    assemblies_ch = SPADES(trimmed_ch)

    QUAST(assemblies_ch)

    CHECKM2(assemblies_ch)

    SISTR(assemblies_ch)

    annotations_ch = BAKTA(assemblies_ch)

    gff_ch = annotations_ch.map { sample, gff -> gff }

    PANAROO(gff_ch.collect())

    RAXML(PANAROO.out.aln)
}
