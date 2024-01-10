workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

library(DESeq2)
library(GOplot)

input <- 'DEG_results_males_sinFamilia/all_genes.csv'
res <- read.csv(input, header = TRUE, sep = ";")

colores <- c('red', 'green', 'blue', 'orange', 'yellow')

allgenes <- as.data.frame(res$log2FoldChange)
rownames(allgenes) <- res$x
colnames(allgenes) <- 'logFC'
allgenes$gene <- res$x

#Bubble plot of 68 categories validated from GSEA
circ0 <- read.csv('tabla_males_bubble.txt', header = T, sep="\t")

#To compute the Z-score
cats <- as.data.frame(circ0$term)
cats$gene <- circ0$gene
cats$set_size <- circ0$count
colnames(cats) <- c('term', 'gene', 'set_size')

# Assuming you have loaded your data frames 'cats' and 'allgenes'

# Create an empty vector to store Z-scores
z_scores <- numeric(nrow(cats))

# Iterate through each row of 'cats' and calculate the Z-score
for (i in 1:nrow(cats)) {
  # Extract the list of genes in the current term and split it into a vector
  genes_in_term <- unlist(strsplit(cats$gene[i], "/"))
  
  # Filter 'allgenes' to get log2 fold changes for genes in the current term
  log2foldchanges_in_term <- allgenes$logFC[allgenes$gene %in% genes_in_term]

  # Calculate the Z-score for the term
  upregulated_count <- sum(log2foldchanges_in_term > 0)
  downregulated_count <- sum(log2foldchanges_in_term < 0)
  set_size <- length(log2foldchanges_in_term)
  
  z_score <- (upregulated_count - downregulated_count) / sqrt(set_size)
  
  # Store the calculated Z-score in the vector
  z_scores[i] <- z_score
}

# Add the Z-scores to the 'cats' data frame
cats$ZScore <- z_scores

#Relate to the circ data frame
circ0$zscore <- cats$ZScore[cats$term == circ0$term]
circ <- circ0

#Z-score P value GSEA
trace(GOBubble, edit = T)

jpeg(file = 'prueba2.jpeg', units = 'in', width = 15, height = 20, res = 300)
par(mar = c(2, 2, 2, 5)) 
GOBubble(circ, labels = 2.5, colour = colores, ID =T, table.col=F, table.legend = F)

invisible(dev.off())
