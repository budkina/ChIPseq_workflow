library(ChIPpeakAnno)
library(biomaRt)
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl")#host="asia.ensembl.org"

ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
bm <- useDataset("hsapiens_gene_ensembl", mart=ensembl)

Annotation <- getAnnotation(mart)
library(GenomicFeatures)

txdb <- makeTxDbFromGFF('Homo_sapiens.GRCh38.105.gtf')
library(org.Hs.eg.db)
library(ggplot2)

no_replicates <- 1:16

no_replicates_names <- c('FGD5-AS1_04', 'FGD5-AS1_AD03', 'EMX2OS_AD02', 'EMX2OS_AD04', 'AC005592.2_02', 'AC005592.2_03', 'RP11-398K22.12_03', 'RP11-398K22.12_05', 'CTD-2587H24.5_01', 'CTD-2587H24.5_03', 'JPX_05', 'JPX_AD04', 'MAPKAPK5-AS1_06_H3K27me3', 'MAPKAPK5-AS1_06_H3K4me3', 'MAPKAPK5-AS1_AD04_H3K27me3', 'MAPKAPK5-AS1_AD04_H3K4me3')
macs2_peak_types <- c('narrowPeak', 'broadPeak')
for (peak_type in macs2_peak_types)
{
  for ( i in 1:length(no_replicates))
  {
    rep <- no_replicates[i]
    aso_name <- no_replicates_names[i]
    
    narrow_ranges <- toGRanges( paste0("peaks\\C",
                                       rep,
                                       "_peaks.",peak_type,".filtered"),
                                format=peak_type,
                                header=FALSE)
    
    seqlevelsStyle(narrow_ranges) <- seqlevelsStyle(Annotation)
    
    narrow_ranges.both.anno <- annotatePeakInBatch(narrow_ranges, 
                                                   AnnotationData=Annotation,
                                                   output="both",
                                                   maxgap=3000L)
    
    narrow_ranges.promoter.anno <- annotatePeakInBatch(narrow_ranges, 
                                                       AnnotationData=Annotation, 
                                                       output="nearestBiDirectionalPromoters",
                                                       bindingRegion=c(-3000, 3000))
    
    narrow_ranges.both.anno.out <- unname(narrow_ranges.both.anno)
    narrow_ranges.promoter.anno.out <- unname(narrow_ranges.promoter.anno)
    
    narrow_ranges.both.anno.out$peakNames <- NULL # remove the CharacterList to avoid error message.
    narrow_ranges.promoter.anno.out$peakNames <- NULL # remove the CharacterList to avoid error message.
    
    write.table(data.frame(narrow_ranges.both.anno.out), 
                paste0("ChIPpeakAnno//", aso_name, ".",peak_type,".both.csv"),
                sep="\t", 
                quote=F, 
                row.names=F)
    
    write.table(data.frame(narrow_ranges.promoter.anno.out), 
                paste0("ChIPpeakAnno//", aso_name, ".",peak_type,".promoter.csv"),
                sep="\t", 
                quote=F, 
                row.names=F)
    
  }
}

######### macs2

macs2_peak_types <- c('narrowPeak', 'broadPeak')

replicates_rep1 <- c(17,19,21,23,25,26)
replicates_rep2 <- c(18,20,22,24,27,28)

control_names <- c('NCA_H3K9me3', 'NCA_H3K36me3', 'NCA_H3K9ac', 'NCA_H3K27ac', 'NCA_H3K27me3', 'NCA_H3K4me3') 
for (peak_type in macs2_peak_types)
{
  for (i in 1:length(replicates_rep1))
  {
    rep1 <- replicates_rep1[i]
    rep2 <- replicates_rep2[i]
    control_name <- control_names[i]
    rep1_ranges <- toGRanges( paste0("peaks\\C",
                                       rep1,
                                       "_peaks.",peak_type,".filtered"),
                                format=peak_type,
                                header=FALSE)
    rep2_ranges <- toGRanges( paste0("peaks\\C",
                                     rep2,
                                     "_peaks.",peak_type,".filtered"),
                              format=peak_type,
                              header=FALSE)
    
    ol <- findOverlapsOfPeaks(rep1_ranges, rep2_ranges, connectedPeaks = 'keepAll')
    ol <- addMetadata(ol, colNames="score", FUN=mean) 
    
    seqlevelsStyle(rep1_ranges) <- seqlevelsStyle(Annotation)
    seqlevelsStyle(rep2_ranges) <- seqlevelsStyle(Annotation)
    overlaps <- ol$peaklist[["rep1_ranges///rep2_ranges"]]
    seqlevelsStyle(overlaps) <- seqlevelsStyle(Annotation)

    narrow_ranges.both.anno <- annotatePeakInBatch(overlaps, 
                                                   AnnotationData=Annotation,
                                                   output="both",
                                                   maxgap=3000L)
    
    narrow_ranges.promoter.anno <- annotatePeakInBatch(overlaps, 
                                                       AnnotationData=Annotation, 
                                                       output="nearestBiDirectionalPromoters",
                                                       bindingRegion=c(-3000, 3000))
    
    narrow_ranges.both.anno.out <- unname(narrow_ranges.both.anno)
    narrow_ranges.promoter.anno.out <- unname(narrow_ranges.promoter.anno)
    
    narrow_ranges.both.anno.out$peakNames <- NULL # remove the CharacterList to avoid error message.
    narrow_ranges.promoter.anno.out$peakNames <- NULL # remove the CharacterList to avoid error message.
    
    write.table(data.frame(narrow_ranges.both.anno.out), 
                paste0("ChIPpeakAnno//", control_name, ".",peak_type,".both.csv"),
                sep="\t", 
                quote=F, 
                row.names=F)
    
    write.table(data.frame(narrow_ranges.promoter.anno.out), 
                paste0("ChIPpeakAnno//", control_name, ".",peak_type,".promoter.csv"),
                sep="\t", 
                quote=F, 
                row.names=F)
    
  }
}

######### sicer
no_replicates <- c("1_","2_","3_","4_","5_","6_","7_","8_","9_","10","11","12","13","14","15","16")

for ( i in 1:length(no_replicates))
{
  rep <- no_replicates[i]
  aso_name <- no_replicates_names[i]

  df <- read.csv(paste0("peaks\\C-",rep,".sicer.filtered"), header = F,
                 sep = "\t", quote = "")
  df <-df[c('V1', 'V2', 'V3')]
  
  colnames(df) <- c('seqnames',  'start', 'end')
  
  narrow_ranges <- makeGRangesFromDataFrame(df,
                                        keep.extra.columns=FALSE,
                                        ignore.strand=T,
                                        seqinfo=NULL,
                                        seqnames.field=c("seqnames"),
                                        start.field="start",
                                        end.field="end")
  
  seqlevelsStyle(narrow_ranges) <- seqlevelsStyle(Annotation)
  
  narrow_ranges.both.anno <- annotatePeakInBatch(narrow_ranges, 
                                                 AnnotationData=Annotation,
                                                 output="both",
                                                 maxgap=3000L)
  
  narrow_ranges.promoter.anno <- annotatePeakInBatch(narrow_ranges, 
                                                     AnnotationData=Annotation, 
                                                     output="nearestBiDirectionalPromoters",
                                                     bindingRegion=c(-3000, 3000))
  
  narrow_ranges.both.anno.out <- unname(narrow_ranges.both.anno)
  narrow_ranges.promoter.anno.out <- unname(narrow_ranges.promoter.anno)
  
  narrow_ranges.both.anno.out$peakNames <- NULL # remove the CharacterList to avoid error message.
  narrow_ranges.promoter.anno.out$peakNames <- NULL # remove the CharacterList to avoid error message.
  
  write.table(data.frame(narrow_ranges.both.anno.out), 
              paste0("ChIPpeakAnno//", aso_name, ".sicer.both.csv"),
              sep="\t", 
              quote=F, 
              row.names=F)
  
  write.table(data.frame(narrow_ranges.promoter.anno.out), 
              paste0("ChIPpeakAnno//", aso_name, ".sicer.promoter.csv"),
              sep="\t", 
              quote=F, 
              row.names=F)
  
}


replicates_rep1 <- c('17','19','21','23','25','26')
replicates_rep2 <- c('18','20','22','24','27','28')

for ( i in 1:length(replicates_rep1))
{
  rep1 <- replicates_rep1[i]
  rep2 <- replicates_rep2[i]
  control_name <- control_names[i]
  df <- read.csv(paste0("peaks\\C-",rep1,".sicer.filtered"), header = F,
                 sep = "\t", quote = "")
  df <-df[c('V1', 'V2', 'V3')]
  
  colnames(df) <- c('seqnames',  'start', 'end')
  
  rep1_ranges <- makeGRangesFromDataFrame(df,
                                            keep.extra.columns=FALSE,
                                            ignore.strand=T,
                                            seqinfo=NULL,
                                            seqnames.field=c("seqnames"),
                                            start.field="start",
                                            end.field="end")
  
  df <- read.csv(paste0("peaks\\C-",rep2,".sicer.filtered"), header = F,
                 sep = "\t", quote = "")
  df <-df[c('V1', 'V2', 'V3')]
  
  colnames(df) <- c('seqnames',  'start', 'end')
  
  rep2_ranges <- makeGRangesFromDataFrame(df,
                                          keep.extra.columns=FALSE,
                                          ignore.strand=T,
                                          seqinfo=NULL,
                                          seqnames.field=c("seqnames"),
                                          start.field="start",
                                          end.field="end")
  
  ol <- findOverlapsOfPeaks(rep1_ranges, rep2_ranges, connectedPeaks = 'keepAll')

  seqlevelsStyle(rep1_ranges) <- seqlevelsStyle(Annotation)
  seqlevelsStyle(rep2_ranges) <- seqlevelsStyle(Annotation)
  overlaps <- ol$peaklist[["rep1_ranges///rep2_ranges"]]
  seqlevelsStyle(overlaps) <- seqlevelsStyle(Annotation)
  
  narrow_ranges.both.anno <- annotatePeakInBatch(overlaps, 
                                                 AnnotationData=Annotation,
                                                 output="both",
                                                 maxgap=3000L)
  
  narrow_ranges.promoter.anno <- annotatePeakInBatch(overlaps, 
                                                     AnnotationData=Annotation, 
                                                     output="nearestBiDirectionalPromoters",
                                                     bindingRegion=c(-3000, 3000))
  
  narrow_ranges.both.anno.out <- unname(narrow_ranges.both.anno)
  narrow_ranges.promoter.anno.out <- unname(narrow_ranges.promoter.anno)
  
  narrow_ranges.both.anno.out$peakNames <- NULL # remove the CharacterList to avoid error message.
  narrow_ranges.promoter.anno.out$peakNames <- NULL # remove the CharacterList to avoid error message.
  
  write.table(data.frame(narrow_ranges.both.anno.out), 
              paste0("ChIPpeakAnno//", control_name, ".sicer.both.csv"),
              sep="\t", 
              quote=F, 
              row.names=F)
  
  write.table(data.frame(narrow_ranges.promoter.anno.out), 
              paste0("ChIPpeakAnno//", control_name, ".sicer.promoter.csv"),
              sep="\t", 
              quote=F, 
              row.names=F)
}

de<-read.csv('fantom6_signatures//RP11-398K22.12_ASO_G0229852_05.csv', header = T, sep = '\t')
de_genes <- de$geneID
anno_files <- c('RP11-398K22.12_05.broadPeak.both', 'RP11-398K22.12_05.broadPeak.promoter', 'RP11-398K22.12_05.narrowPeak.both', 'RP11-398K22.12_05.narrowPeak.promoter', 'RP11-398K22.12_05.sicer.both', 'RP11-398K22.12_05.sicer.promoter')
for (file in anno_files)
{
  anno <- read.csv(paste0("ChIPpeakAnno//",file,".csv"), sep = '\t', header = T, quote = "")
  de_anno <- anno[anno$feature %in% de_genes,]
  write.table(de_anno, 
              paste0("ChIPpeakAnno_de//",file,".de_only.csv"),
              sep="\t", 
              quote=F, 
              row.names=F)
  
}


asos <- c('aso1', 'aso2')
lncrnas <- c('JPX', 'AC005592', 'FGD5', 'EMX2OS', 'CTD', 'RP11', 'MAPKAPK5_8580', 'MAPKAPK5_6002')
types <- c('narrow','broad','sicer')
for (rna in lncrnas)
{
  for (aso in asos)
  {
    for (t in types)
    {
      df <- read.csv(paste0("diffbind\\",rna ,"_", aso,"_no_adj.diffbind_", t),
                     sep = "\t", quote = "")

      diff_bind <- makeGRangesFromDataFrame(df,
                                            keep.extra.columns=FALSE,
                                            ignore.strand=T,
                                            seqinfo=NULL,
                                            seqnames.field=c("seqnames"),
                                            start.field="start",
                                            end.field="end")
      
      seqlevelsStyle(diff_bind) <- seqlevelsStyle(Annotation)
      
      narrow_ranges.both.anno <- annotatePeakInBatch(diff_bind, 
                                                     AnnotationData=Annotation,
                                                     output="both",
                                                     maxgap=3000L)
      
      narrow_ranges.promoter.anno <- annotatePeakInBatch(diff_bind, 
                                                         AnnotationData=Annotation, 
                                                         output="nearestBiDirectionalPromoters",
                                                         bindingRegion=c(-3000, 3000))
      
      narrow_ranges.both.anno.out <- unname(narrow_ranges.both.anno)
      narrow_ranges.promoter.anno.out <- unname(narrow_ranges.promoter.anno)
      
      narrow_ranges.both.anno.out$peakNames <- NULL # remove the CharacterList to avoid error message.
      narrow_ranges.promoter.anno.out$peakNames <- NULL # remove the CharacterList to avoid error message.
      
      write.table(data.frame(narrow_ranges.both.anno.out), 
                  paste0("diffbind_ChIPpeakAnno//", rna, "_", aso,"_",t, ".both.csv"),
                  sep="\t", 
                  quote=F, 
                  row.names=F)
      
      write.table(data.frame(narrow_ranges.promoter.anno.out), 
                  paste0("diffbind_ChIPpeakAnno//", rna, "_", aso,"_",t, ".promoter.csv"),
                  sep="\t", 
                  quote=F, 
                  row.names=F)
    }
  }
}



de<-read.csv('fantom6_signatures//MAPKAPK5-AS1_ASO_G0234608_AD_04.csv', header = T, sep = '\t')
de_genes <- de$geneID
anno_files <- c("MAPKAPK5-AS1_AD04_H3K4me3_broad.both", "MAPKAPK5-AS1_AD04_H3K4me3_broad.promoter", "MAPKAPK5-AS1_AD04_H3K4me3_narrow.both", "MAPKAPK5-AS1_AD04_H3K4me3_narrow.promoter", "MAPKAPK5-AS1_AD04_H3K4me3_sicer.both", "MAPKAPK5-AS1_AD04_H3K4me3_sicer.promoter", "MAPKAPK5-AS1_AD04_H3K27me3_broad.both", "MAPKAPK5-AS1_AD04_H3K27me3_broad.promoter", "MAPKAPK5-AS1_AD04_H3K27me3_narrow.both", "MAPKAPK5-AS1_AD04_H3K27me3_narrow.promoter", "MAPKAPK5-AS1_AD04_H3K27me3_sicer.both", "MAPKAPK5-AS1_AD04_H3K27me3_sicer.promoter")
for (file in anno_files)
{
  anno <- read.csv(paste0("diffbind_ChIPpeakAnno//",file,".csv"), sep = '\t', header = T, quote = "")
  de_anno <- anno[anno$feature %in% de_genes,]
  write.table(de_anno, 
              paste0("diffbind_ChIPpeakAnno_de//",file,".de_only.csv"),
              sep="\t", 
              quote=F, 
              row.names=F)
}
