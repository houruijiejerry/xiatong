library(Seurat)
library(foreach)
library(reshape2)
library(ggplot2)
library(psych)
library(RColorBrewer)
library(harmony)
library(umap)
sample_cols = unique(c(brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3")))
defined_cols = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', sample_cols)
Data = readRDS("../Data.rds")
ALL=merge(Data1,c(Data2,Data3,Data4,Data5,Data6,Data7,))
ALL[["percent.mt"]] <- PercentageFeatureSet(ALL, pattern = "^MT-")
ALL <- subset(ALL, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA> 400 & nCount_RNA < 40000)
ALL = NormalizeData(ALL, verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(npcs = 100, verbose = FALSE)
###########
ALL = RunHarmony(ALL, group.by.vars="stim") %>% 
	RunUMAP(reduction="harmony", dim=1:20) %>% 
	FindNeighbors(reduction="harmony", dim=1:20) %>%
	FindClusters(resolution=0.6)
DimPlot(ALL, reduction = "umap", label = TRUE, cols=defined_cols)
DimPlot(ALL, reduction = "umap", label = FALSE, group.by="group", cols=sample_cols)
DimPlot(ALL, reduction = "umap", label = FALSE, group.by="stim") + guides(color=FALSE)
###########
DotPlot(ALL, features = c("MS4A2","CPA3","TPSAB1","MS4A1", "BANK1", "CD79A","CD3E","CD3D","GZMH","MARCO",
                                  "MSR1", "MRC1", "CD14","S100A8","FCN1","AGER","CLIC5","SFTPB","SFTPD","ETV5",
                                  "FOXJ1","TP73","CCDC78","CLDN5","VWF","LYVE1","PROX1","COL1A1", "PDGFRA", "DCN","ACTA2"),dot.scale =12)
###########
features = c("FGFR1", "FGFR2", "FGFR3","FLT1","FLT4","KDR",
             "PDGFRA","PDGFRB")
VlnPlot(ALL, features = features, pt.size = F, ncol = 4,raster=FALSE)
##########
FeaturePlot(object = ALL, features = c("COL1A1", "FGFR1"), 
            reduction = 'umap', pt.size = 2, blend = T, cols = c("grey", "red", "green", "blue"), min.cutoff = 0, max.cutoff = 5)
############cell ratio
table(Idents(st), st$group)
prop.table(table(Idents(st), st$group), margin = 2)
cell.prop<-as.data.frame(prop.table(table(Idents(st), st$group)))
colnames(cell.prop) <- c("Var1","Var2","Freq")
ggplot(cell.prop,aes(Var2,Freq,fill=Var1))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,"cm"))+
  guides(fill=guide_legend(title=NULL))  
#######################
saveRDS(ALL,file="ALL.rds")
#############################################AddModuleScore##############################
library(Seurat)
library(ggplot2)
########################
sensitive.genes = c("ARF4","PHLDA1","TMED2","RPN1","TNFAIP6",
                    "PRDX6","C5orf15","EGR1","TPBG","CALR",
                    "KCTD12","SSR1","FKBP14","DUSP6","THBS1")
##
sensitive.genes = c("YIPF5","CALR","PDIA6","RPN1","TRAM1","KCTD12","ARF4","TMED2",
                    "PHLDA1","IKBIP","TMEM167A","HOXB5","MCL1","PRDX6","SRPRB","AQP1",
                    "DNAJB11","PDIA4","DYNLT3","C5orf15","TGM2","TPBG","SSR1","FKBP14",
                    "CD82","SLC38A5","ENC1","TNFAIP6","EGR1","DUSP6")
########################
DefaultAssay(Fibro) <- 'RNA'
Fibro$CellType <- Idents(Fibro)
sensitive_genesets = list(diffgenes = sensitive.genes)
signature_names = c("sensitive_genesets")
names(signature_names) = paste("sensitive", 1:length(sensitive_genesets), sep="")
st = AddModuleScore(object = Fibro, features = sensitive_genesets, ctrl = 100, name = "sensitive")
Idents(st)="seurat_clusters"
DimPlot(st, reduction = "umap", label = F, cols=defined_cols)
########################
umaps = FeaturePlot(object = st, names(signature_names), reduction = "umap", combin=T, min.cutoff="q1")
##############################################GO analysis###############################
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
rt=read.table("symbol.txt",sep="\t",check.names=F,header=T)    
genes=as.vector(rt[,1])
#################
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  ##########human
##################
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)    
#######################KEGG############################
rt=read.table("id.txt",sep="\t",header=T,check.names=F)      
rt=rt[is.na(rt[,"entrezID"])==F,]                           
gene=rt$entrezID
#####
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)   
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                          
barplot(kk, drop = TRUE, showCategory = 20)
#######################GO###########################
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
