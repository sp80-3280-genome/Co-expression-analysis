
#=====================================================================================
#
#  Co-expression analysis: Souza and Van Sluys at al. 2019 (GigaScience)
#
#=====================================================================================


# Loading sample infromation
sample_annot <- read.csv("sample_annot.csv", header = T)

# Loading transcript annotation
gmt_in = read.csv("annot_sucestfun.csv", header=T, sep=";")

# Loading expression matrix
caneData = read.delim("G&M_C2_new_T.csv", header = T, sep=";", dec = ",")

# Editing the expression matrix for co-expression analysis
rownames(caneData) = caneData[,3]
caneData = caneData[,-c(1:4)]
names(caneData) = make.names(sample_annot$SampleName)
save(caneData, file = "caneData.RData")

# Run cemitool
library("CEMiTool")
library(ggplot2)

cem <- cemitool(caneData, sample_annot, gmt_in, apply_vst = T, cor_method = "spearman",
                filter=TRUE, plot=TRUE, verbose=TRUE, diss_thresh = 0.9)

save(cem, file = "cem.RData")

# Create report as pdf and html documents
generate_report(cem, directory="./Report",
                output_format=c("pdf_document", "html_document"))

# Create diagnostic report
diagnostic_report(cem, directory="./Diagnostics")

# Write analysis results into files
write_files(cem, directory="./Tables")

# Save all plots
save_plots(cem, "all", directory="./Plots")

# Module GO enrichment analysis using TopGO

module = read.delim("./Tables/module.tsv", header = T)
MIG = as.vector(module$genes)

SAS_annot = read.delim("SAS_annotation_complete_cane_reg_net_annex.txt", header = T)

module_dscrptn = subset(SAS_annot, SeqName %in% MIG)
module_dscrptn = module_dscrptn[,2:3]

write.csv(module_dscrptn, file = "module_annot.csv")

library('topGO')
geneID2GO = readMappings('SAS_GO_IDs_complete_cane_reg_net_edit')
geneNames<- names(geneID2GO)

# Create a folder named 'topGo' before running this!

for (i in 1:(length(unique(module$modules)) - 1)) {
  M = as.vector(subset(module, modules == paste("M", i, sep = "")))
  myInterestingGenes = M$genes
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  #str(geneList)
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  first <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GOdata, classic = first)
  write.csv(allRes, file = paste("./topGo/M", i, "_topGO.csv", sep=""))
}

# Module heatmap

library(gplots)
load("cem.RData")
filt_expr = t(cem@expression)

info_add = as.matrix(read.csv("info_add.csv", header = T))
info_add[info_add == "L1"] = "#018571"
info_add[info_add == "I1"] = "#80cdc1"
info_add[info_add == "I5"] = "#dfc27d"
info_add[info_add == "I9"] = "#756bb1"
info_add = as.data.frame(info_add)

# Select modules
module = read.delim("./Tables/module.tsv", header = T)

# Create a folder named 'Heatmap' before running this!

for(i in 1:(length(unique(module$modules)) - 1)){
  M = as.vector(subset(module, modules == paste("M", i, sep = "")))
  myInterestingGenes = M$genes
  # Scale expresson values
  DOG = filt_expr[,colnames(filt_expr) %in% myInterestingGenes]
  write.csv(DOG, file = paste("./Heatmap/Expr_table_", paste("M", i, sep=""), ".csv", sep=""))
  DOG = as.matrix(DOG)
  #DOG[is.na(DOG)] <- 0
  #DOG = t(DOG)
  # Export heatmap pdf file
  pdf(file = paste("./Heatmap/heatmap_", paste("M", i, sep=""), ".pdf", sep=""))
  heatmap.2(DOG, trace = 'none', key = T, density.info = "none", #scale = "none",
            #dendrogram='row',
            col=colorRampPalette(c("#ffffcc", "#fed976", "#fd8d3c", "#e31a1c", "#800026")),
            main = paste("M", i, sep = ""), margins=c(10,10),
            RowSideColors = as.character(info_add$Class),
            distfun = function(x) dist(x,method = 'euclidean'))
  dev.off()
}

# Co-expression network files for visualization using Cytoscape

annot = read.delim("SAS_annotation_complete_cane_reg_net_annex.txt", header = T)
annot = annot[,2:3]

module = read.delim("./Tables/module.tsv", header = T)

module_annot = subset(annot, SeqName %in% module$genes)
names(module_annot) = c("genes", "Description")

library("igraph")
library("psych")

load("cem.RData")
filt_expr = cem@expression

module = read.delim("./Tables/module.tsv", header = T)

M = subset(module, modules != "Not.Correlated")
myInterestingGenes = as.vector(M$genes)

M_expr = subset(filt_expr, rownames(filt_expr) %in% myInterestingGenes)

M_corr <- corr.test( t(M_expr), method="spearman", ci=F, adjust="fdr")

M_corr$p[lower.tri( M_corr$p,diag=TRUE)]=NA
Pval.adj <- as.data.frame(as.table(M_corr$p))

M_corr$r [lower.tri( M_corr$r,diag=TRUE)]=NA
Correlation <- as.data.frame(as.table(M_corr$r))

Cor.table <- na.exclude(cbind( Correlation, Pval.adj))[,c(1,2,3,6)]

colnames(Cor.table) <- c("gene1","gene2","cor","p.adj")

Cor.table.filt <- Cor.table [(abs(Cor.table[,3])>0.8 & Cor.table[,4] <0.01 ),]

p.adj <- Cor.table.filt[,4]
p.adj[p.adj==0] <- as.numeric(unlist(format(.Machine)))[1]
Cor.table.filt <- cbind(Cor.table.filt, log10(p.adj))
write.table(Cor.table.filt, "General_Cor.table.filter.txt",
            sep="\t", row.names=F, quote=F)

g <- graph.data.frame( Cor.table.filt[,1:2], directed=FALSE)
degree <- degree(g)
betweenness <- betweenness(g)
Node_nw_st <- data.frame( degree, betweenness)

Rank_stat <- rowMeans(cbind(rank(Node_nw_st[,1]), rank(Node_nw_st[,2])))
Node_nw_st <- cbind(Node_nw_st, Rank_stat)
write.table(Node_nw_st, file = "General_Node_nw_st.txt", sep="\t", col.names = NA, quote=F)

rm(list = ls())
