
library(DiffBind)

run_diffbind <- function(id) {
  dbObj <- dba(sampleSheet=paste0("diffbind_design_",id))
  dbObj <- dba.count(dbObj)
  dbObj.contrast  <- dba.contrast(dbObj, categories=DBA_CONDITION, minMembers=2)
  dbObj.analyzed <- dba.analyze(dbObj.contrast, method=DBA_ALL_METHODS)
  res <- dba.report(dbObj.analyzed, th=0.1, bUsePval=TRUE, method=DBA_ALL_METHODS, fold=2)
  out <- as.data.frame(res)
  write.table(out, file=paste0(id,".csv"), sep="\t", quote=F, row.names=F)
}

ids = c("sicer_aso1_H3K27ac","sicer_aso2_H3K27ac","sicer_aso2_H3K36me3","sicer_aso1_H3K9me3","sicer_aso2_H3K9me3","macs_aso1_H3K27ac","macs_aso2_H3K27ac","macs_aso2_H3K36me3","macs_aso1_H3K9me3","macs_aso2_H3K9me3")

for(id in ids) {
  run_diffbind(id)
}
