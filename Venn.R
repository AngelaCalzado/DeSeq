install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)

category <- "C2_KEGG" 

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

input_males <- paste("results_GSEA_males_sinFamilia/stat/", category, "/GSEA_results_males_sinFamilia.txt", sep = "")
input_females <- paste("results_GSEA_females_sinFamilia/stat/", category, "/GSEA_results_females_sinFamilia.txt", sep = "")
#input_all <- paste("results_GSEA_sinFamilia/stat/", category, "/GSEA_results_sinFamilia.txt", sep = "")

males <- read.table(input_males, header = TRUE, sep="\t", row.names = 1)
females <- read.table(input_females, header = TRUE, sep="\t", row.names = 1)
#all <- read.table(input_all, header = TRUE, sep="\t", row.names = 1)

cat_males <- row.names(males)
cat_females <- row.names(females)
#cat_all <-row.names(all)

# Diagrama de Venn
#data <- list(Males = cat_males, Females = cat_females, All = cat_all)
data <- list(Males = cat_males, Females = cat_females)
Venn <- ggVennDiagram(data)
ggsave(paste("Venn_GSEA/", category, "_Venn_males_females.jpeg", sep = ""), Venn)

# Ver que categoria tienen en comun
males_females <- intersect(rownames(males), rownames(females))
#males_all <- intersect(rownames(males), rownames(all))
#females_all <- intersect(rownames(females), rownames(all))
#males_females_all <- Reduce (intersect, list(rownames(males), rownames(females), rownames(all)))

common_males <- males[row.names(males) %in% males_females, ]
write.table(common_males, file = paste("Venn_GSEA/", category,"_common_males.tsv"), sep = "\t", row.names = TRUE)

common_females <- females[row.names(females) %in% males_females, ]
write.table(common_females, file = paste("Venn_GSEA/", category,"_common_females.tsv"), sep = "\t", row.names = TRUE)


# ignorar
gene_list_males <- unlist(strsplit(males["KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS","core_enrichment"], "/"))
gene_list_females <- unlist(strsplit(females["KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS","core_enrichment"], "/"))
num_common <- length(intersect(gene_list_males, gene_list_females))
