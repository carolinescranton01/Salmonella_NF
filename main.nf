workflow {

    Channel
        .fromFilePairs(params.samples, flat: true)
        .set { reads }

    QC(reads)
    ASSEMBLE(QC.out.reads)
    CHECKM(ASSEMBLE.out.contigs)
    BAKTA(ASSEMBLE.out.contigs)
}
