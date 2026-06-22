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

    /*
     * INPUT CHANNEL (clean + explicit)
     */
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    /*
     * QC + TRIMMING
     */
    trimmed_ch = FASTP(read_pairs_ch)

    FASTQC(trimmed_ch)

    /*
     * ASSEMBLY
     */
    assemblies_ch = SPADES(trimmed_ch)

    /*
     * PER-SAMPLE ANALYSIS
     */
    QUAST(assemblies_ch)
    CHECKM2(assemblies_ch)
    SISTR(assemblies_ch)

    /*
     * ANNOTATION
     */
    bakta_ch = BAKTA(assemblies_ch)

    /*
     * PAN-GENOME STEP (SAFE AGGREGATION)
     */
    gff_ch = bakta_ch
        .map { sample, gff -> gff }
        .collect()

    panaroo_out = PANAROO(gff_ch)

    /*
     * PHYLOGENY
     */
    RAXML(panaroo_out.aln)
}
