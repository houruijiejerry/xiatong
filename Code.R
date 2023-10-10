# DESCRIPTION ####
# The following code will install all dependencies in R allowing for analysis to be performed in the paper Xin Zhang et al.
# This code was written and developed by Dr.Xin Zhang
# Install and run R v4.2.2 (x64 bit) and RStudio 2023.03.0+386 "Cherry Blossom" (x64 bit) for running this code

# Seurat is a multimodal single Cell RNA seq analysis algorithm created by
# The Satija Lab. For more information please see: https://satijalab.org/seurat/

# INSTALLATION ####
install.packages("CellChat") 
install.packages("ggalluvial")
install.packages("nichenetr")
install.packages("Seurat")
install.packages("SeuratObject")
install.packages("tidyverse")
install.packages("AUCell")
install.packages("clusterProfiler")
install.packages("ggplot2")
install.packages("enrichplot")
install.packages("fgsea")
install.packages("msigdbr")
install.packages("Scissor")
install.packages("clustermole")
install.packages("singscore")
install.packages("foreach")
install.packages("reshape2")
install.packages("psych")
install.packages("harmony")
install.packages("umap9")
install.packages("RColorBrewer")
install.packages("org.Hs.eg.db")

# load librarys ####
library(CellChat)
library(ggalluvial)
library(nichenetr)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(AUCell)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(fgsea)
library(msigdbr)
library(Scissor)
library(clustermole)
library(singscore)
library(foreach)
library(reshape2)
library(psych)
library(harmony)
library(umap9)
library(RColorBrewer)
library(org.Hs.eg.db)

# OBJECT SETUP AND NORMALIZATION ####
sample_cols=unique(c(brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3")))
defined_cols = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', sample_cols)


# STEP 1: Load data ####
# load 10x data
Data1 <- Read10X(data.dir ="Data1")
# load h5 data
Data2 <- Read10X_h5("Data2.h5")

# STEP 2: Create Seurat objects ####
Data1<- CreateSeuratObject(counts =Data1, min.cells = 3, min.features = 200)
# Group specific Metadata addition
Data1@meta.data$group='Control' ### Healthy control
Data2@meta.data$group='Covid'   ### COVID-19 patients
Data3@meta.data$group='DM'      ### COVID-19 patients with diabetes
Data4@meta.data$group='HTN'     ### COVID-19 patients with hypertension
Data5@meta.data$group='HD'      ### COVID-19 patients with diabetes and hypertension

# Step 3: Merge Datasets ####
ALL=merge(Data1,c(Data2,Data3,Data4,Data5,Data6,Data7,...))

# Step 4: QC and selecting cells for further analysis ####
# The operator can add columns to object metadata. This is a great place to stash QC stats
ALL[["percent.mt"]] <- PercentageFeatureSet(ALL, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(ALL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# RNA based cell thresholding
ALL <- subset(ALL, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA> 400 & nCount_RNA < 40000)

# Step 5: Normalizing the data ####
ALL = NormalizeData(ALL, verbose = FALSE)

# Step 6:Identification of highly variable features (feature selection) ####
ALL = FindVariableFeatures(ALL, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ALL), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ALL)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Step 7: Scaling the data ####
ALL = ScaleData(ALL, verbose = FALSE)

# Step 8: Perform linear dimensional reduction ####
ALL = RunPCA(ALL,npcs = 100, verbose = FALSE)

# Step 9: Determine the ‘dimensionality’ of the dataset ####
Elbowplot(ALL)

# Step 10: Batch correction using Harmony ####
# Batch correction based on samples
ALL = RunHarmony(ALL, group.by.vars="sample")

# Step 11: Run non-linear dimensional reduction (UMAP) ####
# Run UMAP calculations
ALL = RunUMAP(ALL, reduction="harmony", dim=1:20)

# Step 12: Cluster analysis ####
ALL = FindNeighbors(ALL, reduction="harmony", dim=1:20)
ALL = FindClusters(ALL, resolution=0.6)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(ALL, reduction = "umap", label = TRUE, cols=defined_cols)
DimPlot(ALL, reduction = "umap", label = FALSE, group.by="group", cols=sample_cols)
DimPlot(ALL, reduction = "umap", label = FALSE, group.by="orig.ident") + guides(color=FALSE)

# Cell type identification
# find markers for every cluster compared to all remaining cells, report only the positive ones
ALL.markers <- FindAllMarkers(ALL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ALL.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# prediction using R packages (clustermole, singscore)
avg_exp_mat <- AverageExpression(ALL)
avg_exp_mat <- as.matrix(avg_exp_mat$RNA)
avg_exp_mat <- log1p(avg_exp_mat)
avg_exp_mat[1:5, 1:5]
enrich_tbl <- clustermole_enrichment(expr_mat = avg_exp_mat, species = "hs")
enrich_tbl %>% filter(cluster == "1") %>% head(15)
write.csv(enrich_tbl,file="cellmarker.csv")

# classical cell marker genes ####
## Endothelials（PECAM1,VWF,CLDN5）
## Macrophage（MARCO,MSR1,MRC1）
## Tcells（CD3E,CD3D,GZMH）
## Mast cells（MS4A2,CPA3,TPSAB1）
## Bcells（MS4A1,BANK1,CD79A）
## Plasma (CD27, SLAMF7)
## Monocytes（CD14,S100A8,FCN1）
## Neutrophil (S100A9, IFITM2, FCGR3B)
## DCs（LILRB4, IRF8,LILRA4)
## Fibroblasts (COL1A1,PDGFRA)
## SMC(CNN1, ACTA2, TAGLN)
## AT1 (AGER,PDPN,CLIC5)
## AT2（SFTPB,SFTPD,ETV5）
## Ciliated（FOXJ1,TP73,CCDC78）
## Basal(KRT5, KRT14, TP63)
## Goblet(MUC5B, MUC5AC, SPDEF)

# Rename Idents
ALL <- RenameIdents(ALL, `0` = "Fibroblasts", `1` = "Tcells", `2` = "Macrophage",  `3` = "Endothelials", 
                    `4` = "Macrophage", `5` = "AT2", `6` = "Macrophage", `7` = "AT1", `8` = "Fibroblasts", 
                    `9` = "Neutrophils", `10` = "Plasma", `11` = "Monocytes", `12` = "Ciliated", `13` = "SMC", 
                    `14` = "Goblet", `15` = "Basal", `16` = "Bcells", `17` = "Endothelials", `18` = "Tcells", `19` = "Mast", 
                    `20` = "mix1", `21` = "mix2", `22` = "Tcells",  `23` = "EndMT", `24` = "mix3", `25` = "AT2",
                    `26` = "Macrophage",`27` = "DCs", `28` = "Endothelials", `29` = "AT2")
DimPlot(ALL, reduction = "umap", label = TRUE, cols=defined_cols)

# visualize gene expression 
features=c(PECAM1,VWF,CLDN5, MARCO,MSR1,MRC1, CD3E,CD3D,GZMH, MS4A2,CPA3,TPSAB1,
           MS4A1,BANK1,CD79A, CD27, SLAMF7, CD14,S100A8,FCN1, S100A9, IFITM2, FCGR3B,
           LILRB4, IRF8,LILRA4, COL1A1,PDGFRA, CNN1, ACTA2, TAGLN, AGER,PDPN,CLIC5,
           SFTPB,SFTPD,ETV5,FOXJ1,TP73,CCDC78, KRT5, KRT14, TP63, MUC5B, MUC5AC, SPDEF)

# Four visualizations of marker feature expression
# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(ALL,features=features,dot.scale =12)

# Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(ALL, features = features, pt.size = F, ncol = 4,raster=FALSE)

# Feature plot - visualize feature expression in low-dimensional space
FeaturePlot(ALL, features = features)

# Single cell heatmap of feature expression
DoHeatmap(subset(ALL, downsample = 100), features = features, size = 3)

# Proportion of cells vary by group ###
r=prop.table(table(Idents(ALL), ALL$orig.ident), margin = 2)
write.csv(r,file="ALL.csv")
saveRDS(Data,file="Data.harmony.rds")

### subset ###
# Endothelials #
setwd("C:\\Users\\lexb4\\Desktop\\scRNA")
Data=readRDS("Data.harmony.rds")
endo=subset(Data, idents = 'endo')
endo = NormalizeData(endo, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(npcs = 100, verbose = FALSE)
dev.off()
ElbowPlot(endo)
endo <- RunUMAP(endo, reduction = "pca",dims = 1:10,min.dist = 0.8)
endo <- FindNeighbors(endo, reduction = "pca", dims = 1:10)
endo <- FindClusters(endo, resolution = 0.8)
DimPlot(endo,label = T)
DotPlot(endo,features = c("GJA5","DKK2","HEY1","ACKR1","SELE","SELP","PROX1","PDPN","IL7R","CA4"),dot.scale = 8)
endo <- RenameIdents(endo, `0` = "rest ECs1", `1` = "rest ECs1", `2` = "rest ECs1", 
                     `3` = "rest ECs2", `4` = "rest ECs3", `5` = "activated ECs", 
                     `6` = "rest ECs4", `7` = "rest ECs1", `8` = "rest ECs2",
                     `9` = "rest ECs4", `10` = "rest ECs1", `11` = "rest ECs1", 
                     `12` = "rest ECs4", `13` = "rest ECs1", `14` = "rest ECs2")
# proportions of Endothelials subtypes
r=prop.table(table(Idents(endo), endo$orig.ident),margin = 2)
write.csv(r,file="endocellratio.csv")

# Tcells #
Tcells=subset(ALL,idents = 'Tcells')
Tcells = NormalizeData(Tcells, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(npcs = 50, verbose = FALSE)
Idents(Tcells)="seurat_clusters"
Tcells<- RunUMAP(Tcells, reduction = "harmony", dims = 1:15)
Tcells <- FindNeighbors(Tcells,reduction="harmony", dim=1:15)
Tcells <- FindClusters(Tcells, resolution = 0.8)
DimPlot(Tcells, reduction = "umap", label = T, cols=defined_cols)
DotPlot(Tcells, features = c("CD8A", "CD4",### CD8
                             "CXCR3", "TNF", "CCR5","CCR3",##Th1
                             "KLRF1", "GNLY", "NCAM1",##NK
                             "FOXP3","CTLA4", "IL2RA", "IL7R",##Treg
                             "CCR8", "CCR4"))### Th2
Tcells <- RenameIdents(Tcells, `0` = "Trg", `1` = "Th2", `2` = "mix",  `3` = "Th2", `4` = "Th2", `5` = "NK/NKT", `6` = "Th1", 
                       `7` = "NK/NKT", `8` = "Th2", `9` = "CD8+")
Tcells<- subset(Tcells,idents=c("mix"),invert=TRUE)
# proportions of T cells subtypes
r=prop.table(table(Idents(Tcells), Tcells$orig.ident),margin = 2)
write.csv(r,file="Tcellsratio.csv")

# FIbroblasts #
Fibro=subset(ALL,idents = 'Fibro')
Fibro = NormalizeData(Fibro, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(npcs = 100, verbose = FALSE)
ElbowPlot(Fibro)
Fibro<- RunUMAP(Fibro, reduction = "harmony", dims = 1:5)
Fibro <- FindNeighbors(Fibro,reduction="harmony", dim=1:5)
Fibro <- FindClusters(Fibro, resolution = 0.8)

Fibro <- RenameIdents(Fibro, `0` = "fibro", `1` = "fibro", `2` = "Mesothelial",  `3` = "Lipofibroblast", 
                      `4` = "Myofibroblast", `5` = "Mesothelial", `6` = "Lipofibroblast", `7` = "fibro", 
                      `8` = "Myofibroblast", `9` = "Proliferative", `10` = "fibro", `11` = "Lipofibroblast", 
                      `12` = "Myofibroblast", `13` = "Mesothelial", `14` = "Myofibroblast", `15` = "mix")
# proportions of Fibroblasts subtypes
r=prop.table(table(Idents(Fibro), Fibro$orig.ident),margin = 2)
write.csv(r,file="Fibroratio.csv")

# Macrophage #
Macro=subset(ALL,idents = 'Macrophage')
Macro = NormalizeData(Macro, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(npcs = 50, verbose = FALSE)
dev.off()
ElbowPlot(Macro)
Macro<- RunUMAP(Macro, reduction = "harmony", dims = 1:7)
Macro <- FindNeighbors(Macro,reduction="harmony", dim=1:7)
Macro <- FindClusters(Macro, resolution = 0.4)
DimPlot(Macro, reduction = "umap", label = T, cols=defined_cols)
DotPlot(Macro, features = c("CD64C","IDO1","SOCSI","CXCL10","IL1B","TNF"),dot.scale = 8)#####M1
DotPlot(Macro, features = c("MRC1","TGM2","CD23",",CCL22","FOLR2","SLC46A1"),dot.scale = 8)######M2
Macro <- RenameIdents(Macro, `0` = "M2", `1` = "M2", `2` = "M1M2",  `3` = "M2", 
                      `4` = "M2", `5` = "M1", `6` = "mix1", `7` = "M2", `8` = "M1", 
                      `9` = "mix", `10` = "M2", `11` = "M2", `12` = "mix", `13` = "M1")
# proportions of Macrophage subtypes
r=prop.table(table(Idents(Macro), Macro$orig.ident),margin = 2)
write.csv(r,file="Macrophage_ratio.csv")
 
# Differentially expressed genes (DEGs)
# Differential gene testing across conditions
# Subtype
Subtype=subset(data,idents = 'Subtype')
DM.subtypeDEG = FindMarkers(Subtype,ident.1 = 'DM',ident.2 = 'Covid',group.by = "group")
write.csv(DM.subtypeDEG,file="DM.subtypeDEG.csv")

HTN.subtypeDEG = FindMarkers(Subtype,ident.1 = 'HTN',ident.2 = 'Covid',group.by = "group")
write.csv(HTN.subtypeDEG,file="HTN.subtypeDEG.csv")

HD.subtypeDEG = FindMarkers(Subtype,ident.1 = 'HD',ident.2 = 'Covid',group.by = "group")
write.csv(HD.subtypeDEG,file="HD.subtypeDEG.csv")

# scissor
phenotype <- read.table(("phenotype.txt"), sep="\t", header=T, row.names=1)
bulk_dataset <- read.table(("bulk_dataset.txt"), sep="\t", header=T, row.names=1)

head(bulk_dataset[,1:10])
head(phenotype)

bulk_dataset=t(bulk_dataset)

colnames(bulk_dataset) == row.names(phenotype)
Idents(sc_dataset)="RNA"
sc_dataset = FindNeighbors(sc_dataset, dims = 1:10, verbose = F)
names(sc_dataset)

class(sc_dataset)

names(phenotype)
names(sc_dataset)
sc_dataset@graphs$RNA_snn

tag <- c('0', '1')

infos4 <- Scissor(as.matrix(bulk_dataset), sc_dataset, phenotype$DLCO, tag=tag, alpha = 0.5, 
                  family = "binomial") 

names(infos4)
length(infos4$Scissor_pos)
length(infos4$Scissor_neg)

Scissor_select <- rep(0, ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos4$Scissor_pos] <- 1
Scissor_select[infos4$Scissor_neg] <- 2

sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 0.1, order = c(2,1))
DimPlot(sc_dataset, reduction = 'umap', group.by = 'group', pt.size = 0.2, order = c(2,1))

# cell proption
Idents(sc_dataset)="scissor"
table(Idents(sc_dataset), sc_dataset$group)
prop.table(table(Idents(sc_dataset), sc_dataset$group), margin = 2)
cell.prop<-as.data.frame(prop.table(table(Idents(sc_dataset), sc_dataset$group)))
colnames(cell.prop) <- c("Var1","Var2","Freq")
ggplot(cell.prop,aes(Var2,Freq,fill=Var1))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,"cm"))+
  guides(fill=guide_legend(title=NULL))  

r=prop.table(table(Idents(sc_dataset), sc_dataset$orig.ident),margin = 2)
write.csv(r,file="scissorfibroratio_sample.csv")

# DEGs of Scissor+ vs Scissor-
Idents(sc_dataset)="scissor"
scissor.marker= FindMarkers(sc_dataset,ident.1 = '1',ident.2 = '2',group.by = 'scissor')
write.csv(scissor.marker,file="scissor.marker.csv")
head(scissor.marker)

# GO and KEGG analysis
rt=read.table("symbol.txt",sep="\t",check.names=F,header=T)    
genes=as.vector(rt[,1])

entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)    

rt=read.table("id.txt",sep="\t",header=T,check.names=F)      
rt=rt[is.na(rt[,"entrezID"])==F,]                           
gene=rt$entrezID
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)   
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                          
barplot(kk, drop = TRUE, showCategory = 20)

rt=read.table("id.txt",sep="\t",header=T,check.names=F)         
rt=rt[is.na(rt[,"entrezID"])==F,]                                
gene=rt$entrezID
go <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(go,file="Data.txt",sep="\t",quote=F,row.names = F)                
barplot(go, drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()


# GSEA analysis ####
# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

# reading in data from deseq2
df = read.csv("DEGs.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)
head(gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", ####keytype of gene(SYMBOL or ENSEMBL)##
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
write.table(gse,file="GSEA_GO.txt",sep="\t",quote=F,row.names = F)            

# GSEA Plot
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[188], geneSetID = "GO:0048511")

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)

# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$X %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$avg_logFC

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "kegg")
write.table(kk2,file="GSEA_KEGG.txt",sep="\t",quote=F,row.names = F)            
dotplot(kk2, showCategory = 10) 

# Extraction of gene-positive cell subtypes ###
gene_count = GetAssayData(Subtype,slot="counts")["gene",]
Subtype@meta.data$gene = "gene-"
Subtype@meta.data$gene[gene_count>0] = "gene+"
Idents(Subtype) = Subtype@meta.data$gene
r=prop.table(table(Idents(Subtype), Subtype$orig.ident),margin = 2)
write.csv(r,file="genepositiveSubtype.csv")

# AddModuleScore ####
genes = c("")
DefaultAssay(Subtype) <- 'RNA'
Subtype$CellType <- Idents(Subtype)
genesets = list(diffgenes = genes)
signature_names = c("genesets")
names(signature_names) = paste("Addmodulegenes", 1:length(genesets), sep="")
st = AddModuleScore(object = Subtype, features = genesets, ctrl = 100, name = "Addmodulegenes")
Idents(st)="seurat_clusters"
DimPlot(st, reduction = "umap", label = F, cols=defined_cols)
dev.off()


# AUCell ####
# Seurat RDS object
Data$CellType <- Idents(Data)
cells_rankings <- AUCell_buildRankings(Data@assays$RNA@data)  #
# set gene set of interest
genes <- c("")
geneSets <- list(geneSet1=genes)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
par(mfrow=c(3,3))
set.seed(123)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
aucs <- as.numeric(getAUC(cells_AUC))
Data$AUC <- aucs
dev.off()

# cell–cell communication ####
# NicheNet
# Perform NicheNet analysis starting from a Seurat object

ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5]
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
head(weighted_networks$lr_sig)
head(weighted_networks$gr)
ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

# receiver
receiver = "Myofibroblasts"
expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

# sender
sender_celltypes = c("Endothelials")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

seurat_obj_receiver= subset(seuratObj, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]])

condition_oi = "Covid"
condition_reference = "Control" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
igand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
DotPlot(immune.combined, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

# Active target gene inference
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

# Receptors of top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

# CellChat ####
# Load data
load("data.rdata")
# Prepare input data for CelChat analysis
data.input <- GetAssayData(data, assay = "RNA", slot = "data")
identity <- subset(data@meta.data, select = "celltype")
# Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "celltype")

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human

# Show the structure of the database
unique(CellChatDB$interaction$annotation)
# "Secreted Signaling" ，"ECM-Receptor"， "Cell-Cell Contact" 

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")# use Secreted Signaling
# set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis ####
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

# Inference of cell-cell communication network ####
# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 3)

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)
write.csv(df.net,file="df.net.csv")
# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network ####
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

mat <- cellchat@net$count
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways ####
# Bubble plot
levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use = 21, targets.use = c(1,7,11,14,17,19,23,24,25,26), remove.isolate = FALSE)
saveRDS(cellchat,file='cellchat.rds')


################################################### #
# END
################################################### #
################################################### #

#### Session Info ####
sessionInfo()
## R version 4.2.2 (2022-10-31)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.5 LTS

## Matrix products: default
## BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

## locale:
## [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
## [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
## [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

## other attached packages:
## [1] harmony_1.0.3                     Rcpp_1.0.11                       RColorBrewer_1.1-3               
## [4] dplyr_1.1.3                       future_1.33.0                     patchwork_1.1.3                  
## [7] Seurat_4.3.0.1                    Signac_1.10.0                     SeuratObject_4.1.4               
## [10] sp_2.1-0                          R.utils_2.12.2                    R.oo_1.25.0                      
## [13] R.methodsS3_1.8.2                 Rbowtie2_2.4.2                    Rsubread_2.12.3                  
## [16] BSgenome.Hsapiens.UCSC.hg38_1.4.5 BSgenome_1.66.3                   rtracklayer_1.58.0               
## [19] Rsamtools_2.14.0                  Biostrings_2.66.0                 XVector_0.38.0                   
## [22] ChIPQC_1.34.1                     BiocParallel_1.32.6               DiffBind_3.8.4                   
## [25] SummarizedExperiment_1.28.0       Biobase_2.58.0                    MatrixGenerics_1.10.0            
## [28] matrixStats_1.0.0                 GenomicRanges_1.50.2              GenomeInfoDb_1.34.9              
## [31] IRanges_2.32.0                    S4Vectors_0.36.2                  BiocGenerics_0.44.0              
## [34] ggplot2_3.4.3                   

##loaded via a namespace (and not attached):
##[1] utf8_1.2.3                                spatstat.explore_3.2-1                   
##[3] reticulate_1.31                           tidyselect_1.2.0                         
##[5] RSQLite_2.3.1                             AnnotationDbi_1.60.2                     
##[7] htmlwidgets_1.6.2                         grid_4.2.2                               
##[9] Rtsne_0.16                                munsell_0.5.0                            
##[11] codetools_0.2-19                          ica_1.0-3                                
##[13] interp_1.1-4                              systemPipeR_2.4.0                        
##[15] miniUI_0.1.1.1                            withr_2.5.1                              
##[17] spatstat.random_3.1-5                     colorspace_2.1-0                         
##[19] progressr_0.14.0                          filelock_1.0.2                           
##[21] rstudioapi_0.15.0                         ROCR_1.0-11                              
##[23] tensor_1.5                                listenv_0.9.0                            
##[25] bbmle_1.0.25                              GenomeInfoDbData_1.2.9                   
##[27] polyclip_1.10-4                           mixsqp_0.3-48                            
##[29] hwriter_1.3.2.1                           bit64_4.0.5                              
##[31] coda_0.19-4                               parallelly_1.36.0                        
##[33] vctrs_0.6.3                               TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2  
##[35] generics_0.1.3                            BiocFileCache_2.6.1                      
##[37] R6_2.5.1                                  apeglm_1.20.0                            
##[39] invgamma_1.1                              locfit_1.5-9.8                           
##[41] spatstat.utils_3.0-3                      bitops_1.0-7                             
##[43] cachem_1.0.8                              DelayedArray_0.24.0                      
##[45] promises_1.2.1                            BiocIO_1.8.0                             
##[47] scales_1.2.1                              gtable_0.3.4                             
##[49] globals_0.16.2                            goftest_1.2-3                            
##[51] rlang_1.1.1                               RcppRoll_0.3.0                           
##[53] splines_4.2.2                             lazyeval_0.2.2                           
##[55] spatstat.geom_3.2-4                       abind_1.4-5                              
##[57] yaml_2.3.7                                reshape2_1.4.4                           
##[59] TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2 GenomicFeatures_1.50.4                   
##[61] httpuv_1.6.11                             tools_4.2.2                              
##[63] ellipsis_0.3.2                            gplots_3.1.3                             
##[65] ggridges_0.5.4                            plyr_1.8.9                               
##[67] progress_1.2.2                            zlibbioc_1.44.0                          
##[69] purrr_1.0.2                               RCurl_1.98-1.12                          
##[71] prettyunits_1.1.1                         deldir_1.0-9                             
##[73] pbapply_1.7-2                             ashr_2.2-63                              
##[75] cowplot_1.1.1                             zoo_1.8-12                               
##[77] chipseq_1.48.0                            ggrepel_0.9.3                            
##[79] cluster_2.1.4                             magrittr_2.0.3                           
##[81] scattermore_1.2                           data.table_1.14.8                        
##[83] TxDb.Hsapiens.UCSC.hg18.knownGene_3.2.2   lmtest_0.9-40                            
##[85] RANN_2.6.1                                truncnorm_1.0-9                          
##[87] mvtnorm_1.2-3                             SQUAREM_2021.1                           
##[89] amap_0.8-19                               fitdistrplus_1.1-11                      
##[91] TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2   hms_1.1.3                                
##[93] mime_0.12                                 xtable_1.8-4                             
##[95] XML_3.99-0.14                             emdbook_1.3.13                           
##[97] jpeg_0.1-10                               gridExtra_2.3                            
##[99] compiler_4.2.2                            biomaRt_2.54.1                           
##[101] bdsmatrix_1.3-6                           tibble_3.2.1                             
##[103] KernSmooth_2.23-21                        crayon_1.5.2                             
##[105] htmltools_0.5.6                           later_1.3.1                              
##[107] tidyr_1.3.0                               DBI_1.1.3                                
##[109] dbplyr_2.3.3                              MASS_7.3-60                              
##[111] rappdirs_0.3.3                            ShortRead_1.56.1                         
##[113] Matrix_1.6-1.1                            cli_3.6.1                                
##[115] parallel_4.2.2                            igraph_1.5.1                             
##[117] pkgconfig_2.0.3                           TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2  
##[119] GenomicAlignments_1.34.1                  numDeriv_2016.8-1.1                      
##[121] spatstat.sparse_3.0-2                     plotly_4.10.2                            
##[123] TxDb.Celegans.UCSC.ce6.ensGene_3.2.2      xml2_1.3.5                               
##[125] stringr_1.5.0                             digest_0.6.33                            
##[127] sctransform_0.3.5                         RcppAnnoy_0.0.21                         
##[129] spatstat.data_3.0-1                       leiden_0.4.3                             
##[131] fastmatch_1.1-3                           uwot_0.1.16                              
##[133] restfulr_0.0.15                           GreyListChIP_1.30.0                      
##[135] curl_5.0.2                                shiny_1.7.5                              
##[137] gtools_3.9.4                              rjson_0.2.21                             
##[139] nlme_3.1-162                              lifecycle_1.0.3                          
##[141] jsonlite_1.8.7                            viridisLite_0.4.2                        
##[143] limma_3.54.2                              fansi_1.0.4                              
##[145] pillar_1.9.0                              lattice_0.21-8                           
##[147] Nozzle.R1_1.1-1.1                         KEGGREST_1.38.0                          
##[149] fastmap_1.1.1                             httr_1.4.7                               
##[151] survival_3.5-5                            glue_1.6.2                               
##[153] png_0.1-8                                 bit_4.0.5                                
##[155] stringi_1.7.12                            blob_1.2.4                               
##[157] TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0 latticeExtra_0.6-30                      
##[159] caTools_1.18.2                            memoise_2.0.1                            
##[161] irlba_2.3.5.1                             future.apply_1.11.0  
