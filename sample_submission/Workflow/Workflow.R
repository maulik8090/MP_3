
library(dplyr)
library(Seurat)
library(hdf5r)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
library(SeuratDisk)
library(harmony)
library(zeallot)

# clear the environment to avoid confusing other objects with this project's objects, and release unused memory
rm(list=ls())
gc()

# Set the working directory of the files with "setwd", or navigate to it in the Files tab and set it under 
# the "More" drop-down button.

# ------------------------------------------------ All-sample merged AHFD -----------------------------------------------------

# To process multiple samples' files at the same time, uniformly, they can be read into R as a list with
# the "list.files" function, and make them into a vector list using "vector" function. 

# Read-in loop(s) for each file:
setwd("Data")

#  # For .h5 extension files out of CellRanger (Note: if filenames are sorted by a number in the name (e.g. s1, s2),
#  # double-digits will preceed single digits -> 1, 10, 11, 2, 3... Rename single digits as "01, 02..." to avoid)
project="AHFD"
files <- list.files()
files
data <- vector(mode = "list", length = length(files))
for(file in files){
  indivs  <- Read10X_h5(file)
  indivs <- CreateSeuratObject(counts = indivs, project=project, min.cells=3, min.feature=200)
  data[files %in% file] <- indivs
}
data


#  # TO read in multiple .h5seurat file types (saved with SeuratDisk package)
setwd("Data")
files <- list.files()
files
data <- vector(mode = "list", length = length(files))
for(file in files){
  indivs  <- LoadH5Seurat(file)
  data[files %in% file] <- indivs
}
rm(indivs)
data




# -------------- Merge all of the samples files together into one database in which sample names are annotated as vector info ------------------
samplist <- c("Smpl1", "Smpl2", "Smpl3", "Smpl4", "Smpl5", "Smpl6", "Smpl7", "Smpl8", "Smpl9", "Smpl10", "Smpl11", "Smpl12")
ahfbase <- merge(x = samps6_formerging[[1]], y = samps6_formerging[2:length(samps6_formerging)], merge.data=TRUE, add.cell.ids=sampnames)
ahfbase <- SetIdent(ahfbase, value = "RNA")
ahfbase <- subset(ahfbase, subset = DF.class == "Singlet")  # subsetting to remove doubletfinder-predicted doublets. "DF.class" is meta.data category I made

# Basic sample/data metrics
ahfbase[["percent.mt"]] <-PercentageFeatureSet(ahfbase, pattern ="^mt-")
ahfbase[["percent.hb"]] <-PercentageFeatureSet(ahfbase, pattern ="^Hb")
ahfbase[["percent.gm"]]<- PercentageFeatureSet(ahfbase, pattern="^Gm") # many lncRNA names are preceded by "Gm"
ahfbase[["percent.rp"]]<- PercentageFeatureSet(ahfbase, pattern="^Rp")


nF_Vln<-VlnPlot(ahfbase, features=c("nFeature_RNA"),pt.size =0)
nC_Vln<-VlnPlot(ahfbase, features=c("nCount_RNA"),pt.size =0)
mthb_Vln<-VlnPlot(ahfbase, features=c("percent.mt", "percent.hb"),pt.size =0, ncol=2,y.max = 1)
rpgm_Vln<-VlnPlot(ahfbase, features=c("percent.rp","percent.gm"),pt.size =0, ncol=2)

plot1 <- FeatureScatter(ahfbase, feature1 = "percent.hb", feature2 = "percent.mt",pt.size = 1)
plot2 <- FeatureScatter(ahfbase, feature1 = "nCount_RNA", feature2 = "percent.gm")
plot3 <- FeatureScatter(ahfbase, feature1 = "nCount_RNA", feature2 = "percent.hb")
plot4 <- FeatureScatter(ahfbase, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

nF_Vln
nC_Vln
mthb_Vln
rpgm_Vln
plot1
plot2
plot3
plot4
rm(nCnF_Vln,
   mthb_Vln,
   nF_Vln,
   plot1,
   plot2,
   plot3,
   plot4,
   plots)

# save plots shown in the plot tab (to the right of here -->), from the temporary folder that R makes for plots
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="Data")
dev.off(dev.list()["RStudioGD"])
 
# ---------------------- Subset to exclude unwanted characteristics, and process for clustering ------------------------------ # 

ahfbase <- subset(ahfbase, subset = nFeature_RNA > 200 )
ahfbase <- NormalizeData(ahfbase, normalization.method = "LogNormalize", scale.factor = 10000)
ahfbase <- FindVariableFeatures(ahfbase, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ahfbase)
ahfbase <- ScaleData(ahfbase, features = all.genes)  # scale data for dimension reduction pre-processing...
ahfbase <- RunPCA(ahfbase, features = VariableFeatures(object = ahfbase))

# skip harmony section if not utilizing...

# ------------------------------------------------------ Harmony integration --------------------------------------------------------- #
# Note: installation into R 2023.06.01 with updated packages throws error:
# Error: element-wise multiplication: incompatible matrix dimensions: 100x11 and 100x1. This problem was fixed by converting the
# target group for harmony function with the as.factor function (e.g. object$sample <- as.factor(object$sample))

ahfbase$Sample <-as.factor(ahfbase$Sample)
data.class(ahfbase$Sample)

# pre-harmony visual
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = ahfbase, reduction = "pca", pt.size = .1, group.by = "Age",raster=FALSE)
p2 <- VlnPlot(object = ahfbase, features = "PC_1", group.by = "Age", pt.size = .1)
p1 + p2

options(repr.plot.height = 2.5, repr.plot.width = 6)
gc(full = T) # circumvents a warning message regarding quick-transfer stage steps exceeding memory
ahfbase <- ahfbase %>% 
  RunHarmony("Sample", plot_convergence = TRUE)

# post-harmony visual
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = ahfbase, reduction = "harmony", pt.size = .1, group.by = "Age", raster=FALSE)
p2 <- VlnPlot(object = ahfbase, features = "harmony_1", group.by = "Age", pt.size = .1)
p1 + p2

# -------------------------------------------- Neighbors, clustering, and Visuals ------------------------------------------------ #

# Include reduction = "harmony" for findneighbors and runumap if harmony has been used.
ahfbase <- FindNeighbors(ahfbase,dims = 1:25)
ahfbase <- FindClusters(ahfbase, resolution = .3)
ahfbase <- RunUMAP(ahfbase, dims = 1:25) # If you haven't installed UMAP, you can do so via reticulate::py_install(packages ='umap-learn')

ElbowPlot(ahfbase, ndims=50) # ElbowPlot can be used to determine dims# for cluster neighboring, instead of the JackStraw command below
top10 <- head(VariableFeatures(ahfbase), 10)  # View variable features (features = expressed genes) for PCA
top10
varplot1 <- VariableFeaturePlot(ahfbase)
varplot2 <- LabelPoints(plot = varplot1, points = top10, repel = TRUE)
varplot2
# VizDimLoadings(orig, dims = 1:2, reduction = "pca")           # Examine and visualize PCA results a few different ways
# DimHeatmap(ahfdall2, dims = 1, cells = 500, balanced = TRUE)   # Examine and visualize PCA results a few different ways

ahfbase <- SetIdent(ahfbase, value = "RNA")

DotPlot(ahfbase,features=c("Nphs1","Fhl2","Emcn","Lrp2","Pth1r","Slc34a1","Agt","Slc5a2","Slc5a12","Slc22a6","Slc7a13","Epha7","Corin","Slc14a2","Aqp1","Slc9a3","Slc12a1","Clcnka","Umod","Slc12a3","Pvalb", "Slc8a1","Aqp2","Kit", "Slc26a4", "Ptprc"))+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=8))+theme(axis.text.y=element_text(size=8))+theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(colour = "gray",linetype ="dashed",size=0.35))
DimPlot(ahfbase, reduction = "umap", group.by="Age", label = TRUE, repel = TRUE,raster=TRUE)
FeaturePlot(ahfbase, reduction = "umap", features = "Igkc", cols=c("pink","black"), raster=TRUE)
VlnPlot(ahfbase, features = "Fabp1", group.by="Age", split.by="Genotype", pt.size = .1)

ggplot(ahfbase@meta.data$Sample, aes(x=RNA_snn_res.0.5, fill=sample)) + geom_bar()
# to barplot clusters by sample content, sample values cannot be integer/numeric.
# need to change the numbers to characters (or factors?), somehow...

# Save plots and clear all plots
Sys.sleep(5)
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics",full.names=TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png",full.names=TRUE)
file.copy(from=plots.png.paths, to="Data")
dev.off(dev.list()["RStudioGD"])

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
VlnPlot(ahfbase, features=c("Nphs1","Nphs2", "Wt1", "Synpo","Cdh5", "Pecam1", "Emcn","Flt1"),pt.size=0, ncol=4)
FeaturePlot(orig, features = "Nphs1", reduction = "umap")
FeaturePlot(ahfbase, features = "Fabp1", split.by= "Genotype", reduction = "umap",cols=c("red","grey"),raster=TRUE)
DotPlot(orig,features=c("wt1","Fhl2","Emcn","Slc5a12","Slc7a13","Epha7","Slc12a1", "Slc12a3","Slc8a1", "Aqp2","Kit", "Slc26a4", "Ptprc"))+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=8))+theme(axis.text.y=element_text(size=8))+theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(colour = "gray",linetype ="dashed",size=0.35), main="Dims1:30, Reso 0.2")

# find all of the markers in each cluster with minimum quantity of .25% of total. WHat is
# the logfc.threshold? -> check the ahfdall data
all.markers <- FindAllMarkers(ahfbase, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers0 <- as.data.frame(all.markers)
write.csv(all.markers0, "2023.3.18_AHFD_base_cluster markers!=0.csv")

# stacked barplot for percentage of samples in clusters
sampcts_inclstrs <- table %>%
  melt(table, id.vars="sample", variable.name="cluster")
  table$sample <- rownames(table) %>%
  as.data.table(table) %>%
  cbind(rownames(table),table) %>%
  apply(table, 2, function(x){x/sum(x)}) %>%
  table(ahfbase$Sample, ahfbase$seurat_clusters)
  
  
#For enrich and DGE analysis 
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichR)
  
  # Set Enrichr site
  setEnrichrSite("Enrichr")
  
  # Define databases for enrichment
  enrichr_dbs <- c("KEGG_2019_Mouse", 
                   "GO_Molecular_Function_2023",
                   "GO_Cellular_Component_2023",
                   "GO_Biological_Process_2023",
                   "WikiPathways_2019_Mouse")
  
  # Set identity to cell type
  Idents(ahfbase) <- ahfbase$celltype  # Replace with correct metadata column if needed
  
  # Define condition column and groups
  condition_col <- "Age"
  group1 <- "24mo"
  group2 <- "3mo"
  
  # Get unique cell types
  celltypes <- unique(ahfbase$celltype)
  
  # Create output lists
  marker_list <- list()
  enrichr_results <- list()
  
  for (ct in celltypes) {
    message("Processing cell type: ", ct)
    
    # Subset to current cell type
    subset_obj <- subset(ahfbase, idents = ct)
    
    # Set condition as identity
    Idents(subset_obj) <- subset_obj[[condition_col]]
    
    # Run differential expression
    markers <- FindMarkers(subset_obj, ident.1 = group1, ident.2 = group2, logfc.threshold = 0.25)
    markers <- markers %>%
      tibble::rownames_to_column("gene") %>%
      arrange(desc(avg_log2FC)) %>%
      mutate(fold = ifelse(avg_log2FC > 0, 2^avg_log2FC, -1 * 2^abs(avg_log2FC))) %>%
      filter(p_val_adj < 0.05) %>%
      filter(!grepl("^Gm", gene)) %>%
      filter(!grepl("Rik", gene))
    
    # Save markers
    marker_list[[ct]] <- markers
    fwrite(markers, paste0("DGE_", ct, "_", group1, "_vs_", group2, ".csv"))
    
    # Top up/down genes
    top_up <- head(markers$gene[markers$avg_log2FC > 0], 50)
    top_down <- head(markers$gene[markers$avg_log2FC < 0], 50)
    
    # Enrichr analysis
    enrich_up <- enrichr(top_up, enrichr_dbs)
    enrich_down <- enrichr(top_down, enrichr_dbs)
    
    # Save results
    enrichr_results[[paste0(ct, "_up")]] <- enrich_up
    enrichr_results[[paste0(ct, "_down")]] <- enrich_down
    
    # Write each database result to CSV
    for (db in enrichr_dbs) {
      write.csv(enrich_up[[db]], paste0("Enrichr_", ct, "_UP_", db, ".csv"))
      write.csv(enrich_down[[db]], paste0("Enrichr_", ct, "_DOWN_", db, ".csv"))
    }
  }
  
