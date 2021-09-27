#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) == 3) {

    stopifnot(file.exists(args[2]))
    stopifnot(file.exists(args[3]))

    library(DSS)
    require(bsseq)

    s1_fn <- args[2]
    s2_fn <- args[3]

    dmr_fn <- paste(args[1], "DMRs.txt", sep='.')
    dml_fn <- paste(args[1], "DMLs.txt", sep='.')

    s1 = read.table(s1_fn, header=TRUE)
    s2 = read.table(s2_fn, header=TRUE)
    BSobj = makeBSseqData(list(s1, s2), c('s1', 's2'))

    dmlTest = DMLtest(BSobj, group1=c('s1'), group2=c('s2'), smoothing=TRUE)
    dmrs = callDMR(dmlTest, p.threshold=0.05)

    write.table(dmrs, file=dmr_fn, quote=F, sep='\t')

    dmls = callDML(dmlTest, p.threshold=1)
    write.table(dmls, file=dml_fn, quote=F, sep='\t')

}

if(length(args) != 3) {
    cat("usage: dss.r <sample name> <DSS_input_1.txt> <DSS_input_2.txt>\n")
}
