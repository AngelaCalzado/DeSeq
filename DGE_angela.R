#DGE
#Run a Differentially Expressed Genes Analysis usng as input the count files 
#produced in previous step

suppressPackageStartupMessages({
  library(rstudioapi)
  library(gplots, quietly = T)
  library(ggplot2, quietly = TRUE)
  library(pheatmap, quietly = TRUE)
  library(DESeq2, quietly = TRUE)
  library('org.Hs.eg.db', quietly = TRUE, character.only = TRUE)
  library(EnhancedVolcano, quietly = TRUE)
})
require("ggrepel", quietly = TRUE)

#Paths
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

configFile <- 'Archivo_configuracion_mPFC.txt'
#Outputs
resD <- 'DEG_results/'
rawCountsF <- paste0(resD,"counts_raw.tsv")
normCountsF <- paste0(resD,"counts_normalized.tsv")
PCAF <- paste0(resD,"PCA.jpeg")
distancesF <- paste0(resD,"distances.jpeg")
dispersionF <- paste0(resD,"dispersion.tiff")
MAplotF <- paste0(resD,"maplot.jpeg")
genesTSV <- paste0(resD,"all_genes.tsv")
sigTSV <- paste0(resD,"sig_pval.tsv")
sigPCAF <- paste0(resD,"PCA_sig.jpeg")
alphasigTSV <- paste0(resD,"0.05_sig_padj.tsv")
volcanoF <- paste0(resD,"volcanoPlot.jpeg")
heatmapF <- paste0(resD,"heatmap.jpeg")
DESEqResultsF <- paste0(resD, 'deseq_objects.RData')

#Parameters
cutoff <- 0.05 #significancy p value adjusted
FCthres <- 4 #fold change threshold to consider in graphs
covs <- T

#Functions
make_italics <- function(x){
  as.expression(lapply(x, function(y) bquote(italic(.(y)))))
}


sampleTable <- read.table(configFile, header=TRUE
                          ,colClasses= c('factor','character','factor',
                                         'factor','factor')
)


data <- DESeqDataSetFromHTSeqCount(sampleTable, directory=".", 
                                   design = ~ Familia + Sexo + Grupo)

keep <- rowSums(counts(data)) >= 10
data <- data[keep,]

dds <- estimateSizeFactors(data)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit = 10000)

dds_raw <- counts(dds, normalized=FALSE)
dds_normalized <- counts(dds, normalized=TRUE)

#save file with counts and normalized counts
write.table(dds_raw, file=rawCountsF, quote=FALSE, 
            sep = "\t", col.names=NA)
write.table(dds_normalized, file=normCountsF, quote=FALSE, 
            sep = "\t", col.names=NA)
save(dds, file = DESEqResultsF)

#PCA
vst <- varianceStabilizingTransformation(dds, blind = FALSE)
mat <- assay(vst)
mm <- model.matrix(~ Sexo + Familia + Grupo, colData(vst))
mat <- limma::removeBatchEffect(mat, batch=vst$Sexo, 
                                batch2=vst$Familia, 
                                group=vst$Grupo)
assay(vst) <- mat

jpeg(filename = PCAF, width=900, height=900, quality=300)
pca <- plotPCA(vst, intgroup = "Sexo")
title <- "Principal Components Plot"
pca + ggtitle(title) + 
  geom_point(size = 6) +
  theme(plot.title = element_text(size=40, hjust = 0.5, face = "bold"), axis.title=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=15)) +
  geom_text_repel(aes(label=colnames(vst)), size=5, point.padding = 0.6)
dev.off()

distRL <- dist(t(mat))
distMat <- as.matrix(distRL)
hc <- hclust(distRL)
hmcol <- colorRampPalette(c("white", "blue"))(299)
title <- "Distances matrix"
jpeg(filename = distancesF, width=900, height=900, quality=300)
heatmap.2(distMat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col=rev(hmcol),
          margin=c(10, 6), main=title, key.title=NA)
invisible(dev.off())

