library(Seurat)
library(nichenetr)
library(tidyverse) 

hra<-read.csv('celseq_matrix_ru10_molecules.tsv.725585', sep='\t', row.names=1)
meta<-read.csv('celseq_meta.tsv.725591', sep = '\t', row.names = 1)
hra[is.na(hra)]<-0
ra<-CreateSeuratObject(hra, project = 'RA', min.cells = 3, min.features = 500, meta.data = meta)

ra[['percent_mt_molecules']]<-as.numeric(ra$percent_mt_molecules)

ra<-subset(ra, subset= percent_mt_molecules < 0.2 & nCount_RNA < 20000)

ra <- NormalizeData(ra, normalization.method = "LogNormalize", scale.factor = 10000)

ra <- FindVariableFeatures(ra, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(ra)
ra <- ScaleData(ra, features = all.genes)

ra <- RunPCA(ra, features = VariableFeatures(object = ra))

VizDimLoadings(ra, dims = 1:2, reduction = "pca")

ra <- JackStraw(ra, num.replicate = 100)
ra <- ScoreJackStraw(ra, dims = 1:20)

ra <- FindNeighbors(ra, dims = 1:20)
ra <- FindClusters(ra, resolution = 0.5)
ra <- RunUMAP(ra, dims = 1:20)

Idents(ra)<- ra$type

ligand_target_matrix <- readRDS("ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5]
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05

lr_network = readRDS("lr_network.rds")
head(lr_network)
##   from  to    source         database
##   <chr> <chr> <chr>          <chr>   
## 1 CXCL1 CXCR2 kegg_cytokines kegg    
## 2 CXCL2 CXCR2 kegg_cytokines kegg    
## 3 CXCL3 CXCR2 kegg_cytokines kegg    
## 4 CXCL5 CXCR2 kegg_cytokines kegg    
## 5 PPBP  CXCR2 kegg_cytokines kegg    
## 6 CXCL6 CXCR2 kegg_cytokines kegg

weighted_networks = readRDS("weighted_networks.rds")
head(weighted_networks$lr_sig)
##   from  to     weight
##   <chr> <chr>   <dbl>
## 1 A1BG  ABCC6  0.422 
## 2 A1BG  ACE2   0.101 
## 3 A1BG  ADAM10 0.0970
## 4 A1BG  AGO1   0.0525
## 5 A1BG  AKT1   0.0855
## 6 A1BG  ANXA7  0.457
head(weighted_networks$gr)
##   from  to     weight
##   <chr> <chr>   <dbl>
## 1 A1BG  A2M    0.0294
## 2 AAAS  GFAP   0.0290
## 3 AADAC CYP3A4 0.0422
## 4 AADAC IRF8   0.0275
## 5 AATF  ATM    0.0330
## 6 AATF  ATR    0.0355

seuratObj = ra

DimPlot(seuratObj, group.by='disease', pt.size=1.5)
DimPlot(SeuratObj, group.by='cellType', pt.size=1.5)


###FIBROBLAST

nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = "Fibroblast", 
  condition_colname = "disease", condition_oi = "RA", condition_reference = "OA", 
  sender = c("T cell", "Monocyte", "B cell"), 
  geneset="up",
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "human")

DotPlot(seuratObj, features=c('TGFB1','IL1B','SPP1','ANXA1','PIK3CB','IFNG','CSF1','IL1RN','ICAM1','CTGF','GAS6','SEMA4D','ITGB7','ADAM17','LGALS3','YARS','ITGAM','HGF'), assay='RNA', cols='RdYlBu', group.by='type')+RotatedAxis()

nichenet_output$ligand_activity_target_heatmap + scale_fill_gradient2(low = "whitesmoke",  high = "royalblue", breaks = c(0,0.003,0.006)) + xlab("RA Response in Monocytes") + ylab("Prioritized immmune cell ligands")


###T CELL

nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = "T cell", 
  condition_colname = "disease", condition_oi = "RA", condition_reference = "OA", 
  sender = c("Monocyte", "Fibroblast", "B cell"), 
  geneset = "up",
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "human")

DotPlot(seuratObj, features=c('SPP1','HLA-DRA','PLAU','DUSP18','TGFB1','CD86','ICAM1','CCL2','APOE','IL16','HMGB1','HGF','PTHLH','LRPAP1','ADAM17','HAS2','FGF2','PTDSS1'), assay='RNA', cols='RdYlBu', group.by='type')+RotatedAxis()

nichenet_output$ligand_activity_target_heatmap + scale_fill_gradient2(low = "whitesmoke",  high = "royalblue", breaks = c(0,0.003,0.006)) + xlab("RA Response in Monocytes") + ylab("Prioritized immmune cell ligands")


###MONOCYTE

nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = "Monocyte", 
  condition_colname = "disease", condition_oi = "OA", condition_reference = "RA", 
  sender = c("B cell", "Fibroblast", "T cell"), 
  geneset = "up",
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "human")

DotPlot(seuratObj, features=c('IFNG','HAS2','TNF','ICAM1','APOE','ANXA1','CCL2','PIK3CB','ADAM17','DUSP18','HLA-DRA','SELP','VEGFA','NOV','TNFSF12','TFPI','KITLG','CD28', 'IL6', 'SPN'), assay='RNA', cols='RdYlBu', group.by='type')+RotatedAxis()

nichenet_output$ligand_activity_target_heatmap + scale_fill_gradient2(low = "whitesmoke",  high = "royalblue", breaks = c(0,0.003,0.006)) + xlab("RA Response in Monocytes") + ylab("Prioritized immmune cell ligands")
