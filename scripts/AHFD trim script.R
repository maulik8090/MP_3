## merge samples 
MPC2 <- merge(SR006706_21, y = c(SR006706_22, SR006706_24, SR006706_25,SR006706_27,SR006706_30 ), 
                        add.cell.ids = c("SR006706_21", "SR006706_22", "SR006706_24","SR006706_25","SR006706_27","SR006706_30"),
                        project = "MPC2")
MPC2$samples <- MPC2$orig.ident                        
## to check the metadata for each sample
table(MPC2$samples)  


# Add group information (example: assign "Control" or "Treatment" based on orig.ident)
MPC2$group <- ifelse(MPC2$orig.ident %in% c("SR006706_21", "SR006706_24"), "Other_group")
MPC2$group[MPC2$orig.ident %in% c("SR006706_22", "SR006706_25")] <- "Group_B"
MPC2$group[MPC2$orig.ident %in% c("SR006706_27", "SR006706_30")] <- "Group_C"
 
# assign experimental groups  as appropriate

SR006706_21[["group"]]<-"Group A"

# Verify the new grouping
table(MPC2$group)

# Step 4: Normalize and preprocess the merged data
MPC2 <- NormalizeData(MPC2, normalization.method = "LogNormalize", scale.factor = 10000)
MPC2 <- FindVariableFeatures(MPC2,selection.method = "vst", nfeatures = 2000)
MPC2 <- ScaleData(MPC2)
MPC2 <- RunPCA(object = MPC2) 
MPC2 <- RunUMAP(object= MPC2, dims = 1:21) 
MPC2 <- FindNeighbors(object = MPC2, dims =1:21) 
MPC2 <- FindClusters(object = MPC2, resolution = .4) 
DimPlot(object = MPC2, reduction = "umap", label = TRUE, raster = TRUE, alpha = 0.7)
ggsave(file= "C:/Users/maulik/Box/Gewin Lab Folder/Maulik/Ageing_MPC2/MPC2/groups_col3_Dimplot_final.pdf")
DotPlot(MPC2,features=Inflammatory_PT,col.min = -0.1,col.max = 3,dot.min = 0,dot.scale = 10,scale.max = 100, cols =
          c("lightyellow","red"))+ scale_x_discrete(limit=rev)+ labs(x=NULL,
                                                                y=NULL,title=print(paste("samples.merged")))+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=14,colour
                                 = "black"), axis.text.y=element_text(size=14,colour = "black"),
        legend.text = element_text(size=14,colour = "black"),
        panel.grid.minor=element_blank(), panel.grid.major=element_line(colour =
                                                                        "gray",linetype ="dashed",size=0.35))
ggsave(file= "C:/Users/maulik/Box/Gewin Lab Folder/Maulik/Ageing_MPC2/MPC2/cellype_Dotplot_refine1.pdf")

saveRDS(MPC2,  "C:/Users/maulik/Box/Gewin Lab Folder/Maulik/Ageing_MPC2/MPC2/MPC2_renamefewcluster.rds")

##Pseudobulk analysis


function (X) 
{
  conditions <- unlist(strsplit(MPC2_finalisedcluster, split = "-", fixed = T))
  for (n in unique(pseudo_bulk$seurat_clusters)) {
    cat(paste("Cluster", n, conditions[[1]], "-", conditions[[2]], 
      "Results"), "\n")
    prev_dir <- getwd()
    cluster_dir <- paste("Cluster_", n, "_", conditions[[1]], 
      "-", conditions[[2]], "_Pseudo_Bulk_Results", sep = "")
    dir.create(cluster_dir)
    setwd(cluster_dir)
    if ((paste(n, conditions[[1]], sep = "_") %in% passing_clusters_by_bulk.condition) & 
      (paste(n, conditions[[2]], sep = "_") %in% passing_clusters_by_bulk.condition)) {
      diff <- FindMarkers(pseudo_bulk, ident.1 = paste(n, 
        conditions[[1]], sep = "_"), ident.2 = paste(n, 
        conditions[[2]], sep = "_"), test.use = "DESeq2", 
        verbose = FALSE)
      diff <- data.frame(gene_name = rownames(diff), diff)
      OUT <- merge(GENE_ANNOTATIONS, diff, by = "gene_name", 
        all.y = TRUE)
      OUT <- OUT[order(OUT$p_val_adj), ]
      outname <- paste("Cluster_", n, "_", conditions[[1]], 
        "-", conditions[[2]], "_diffMarkers", sep = "")
      write.table(OUT, file = paste(outname, ".tsv", sep = ""), 
        row.names = FALSE, col.names = TRUE, sep = "\t", 
        quote = FALSE)
      try(enrich.MSigDb(diff[diff$p_val_adj <= 0.05, ], 
        species = opt$species, out.prefix = outname))
    }
    setwd(prev_dir)
  }
}


## design pairs
function (levels = colnames(design)) 
{
  n <- length(levels)
  design <- matrix(0, n, choose(n, 2))
  rownames(design) <- levels
  colnames(design) <- 1:choose(n, 2)
  k <- 0
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    k <- k + 1
    design[i, k] <- 1
    design[j, k] <- -1
    colnames(design)[k] <- paste(levels[i], "-", levels[j], 
      sep = "")
  }
  design
}


##DEX_by_contrast

function (x, MPC2) 
{
  conditions <- unlist(strsplit(x, split = "-", fixed = T))
  for (n in levels(sample_cluster.table$group)) {
    cat(paste("Cluster", n, conditions[[1]], "-", conditions[[2]], 
      "Results"), "\n")
    prev_dir <- getwd()
    cluster_dir <- paste("Cluster_", n, "_", conditions[[1]], 
      "-", conditions[[2]], "_Results", sep = "")
    dir.create(cluster_dir)
    setwd(cluster_dir)
    if ((paste(n, conditions[[1]], sep = "_") %in% passing_clusters_by_condition) & 
      (paste(n, conditions[[2]], sep = "_") %in% passing_clusters_by_condition)) {
      diff <- FindMarkers(samples.merged.condition, assay = "SCT", 
        ident.1 = paste(n, conditions[[1]], sep = "_"), 
        ident.2 = paste(n, conditions[[2]], sep = "_"), 
        verbose = FALSE)
      diff <- data.frame(gene_name = rownames(diff), diff)
      OUT <- merge(GENE_ANNOTATIONS, diff, by = "gene_name", 
        all.y = TRUE)
      OUT <- OUT[order(OUT$p_val_adj), ]
      outname <- paste("Cluster_", n, "_", conditions[[1]], 
        "-", conditions[[2]], "_diffMarkers", sep = "")
      write.table(OUT, file = paste(outname, ".tsv", sep = ""), 
        row.names = FALSE, col.names = TRUE, sep = "\t", 
        quote = FALSE)
      if (nrow(OUT[OUT$p_val_adj < 0.05, ]) >= 5) {
        tiff(paste(outname, "Top5_genes_tSNE_plots.tiff", 
          sep = "_"), res = 300, height = 5000, width = 5000)
        print(FeaturePlot(object = samples.merged.condition, 
          features = OUT[1:5, 1], split.by = "condition"))
        dev.off()
        tiff(paste(outname, "downSampled_FDR_Signficant_gene_markers_heatmap.tiff", 
          sep = "_"), res = 300, height = 3500, width = 5000)
        try(print(DoHeatmap(subset(samples.merged.condition, 
          downsample = 100), features = OUT[OUT$p_val_adj < 
          0.05, 1], group.by = "cluster.condition", 
          angle = 90, size = 2) + NoLegend()))
        dev.off()
      }
      try(enrich.MSigDb(diff[diff$p_val_adj <= 0.05, ], 
        species = opt$species, out.prefix = outname))
    }
    setwd(prev_dir)
  }
}

## Enrichment analysis

function (x, species = opt$species, out.prefix, clusters = sample.table) 
{
  require(clusterProfiler)
  require(msigdbr)
  if (is.null(species)) {
    orgdb <- NULL
    cat("Species not specified, skipping enrichment analysis....\n")
  }
  else if (species == "human") {
    orgdb <- "org.Hs.eg.db"
    msigdb_species <- "Homo sapiens"
  }
  else if (species == "mouse") {
    orgdb <- "org.Mm.eg.db"
    msigdb_species <- "Mus musculus"
  }
  else if (species == "fly") {
    orgdb <- "org.Dm.eg.db"
    msigdb_species <- "Drosophila melanogaster"
  }
  else if (species == "rat") {
    orgdb <- "org.Rn.eg.db"
    msigdb_species <- "Rattus norvegicus"
  }
  else if (species == "zebrafish") {
    orgdb <- "org.Dr.eg.db"
    msigdb_species <- "Danio rerio"
  }
  else if (species == "worm") {
    orgdb <- "org.Ce.eg.db"
    msigdb_species <- "Caenorhabditis elegans"
  }
  else if (species == "yeast") {
    orgdb <- "org.Sc.eg.db"
    msigdb_species <- "Saccharomyces cerevisiae"
  }
  else {
    orgdb <- NULL
    cat("Species not supported at this time, skipping enrichment analysis....\n")
  }
  MSigDb.categories <- data.frame(collection = c("H", "C1", 
    "C2", "C2", "C2", "C2", "C2", "C2", "C2", "C3", "C3", 
    "C3", "C3", "C4", "C4", "C5", "C5", "C5", "C5", "C6", 
    "C7", "C7", "C8"), category = c("Hallmark_Gene_Sets", 
    "Positional_Gene_Sets", "Chemical_and_Genetic_Perturbations", 
    "Canonical_Pathways", "Canonical_Pathways.Biocarta", 
    "Canonical_Pathways.KEGG", "Canonical_Pathways.PID", 
    "Canonical_Pathways.Reactome", "Canonical_Pathways.Wikipathways", 
    "miRDB_microRNA_Targets", "Legacy_microRNA_Targets", 
    "Transcription_Factor_GTRD_Targets", "Transcription_Factor_Motif_Targets", 
    "Cancer_Gene_Neighborhoods", "Cancer_Modules", "Gene_Ontology.BiologicalProcess", 
    "Gene_Ontology.CellularComponent", "Gene_Ontology.MolecularFunction", 
    "Human_Phenotype_Ontology", "Oncogenic_Signatures", 
    "Immunologic_Signatures", "Vaccine_Response", "Cell_Type_Signatures"), 
    sub.category = c("NULL", "NULL", "CGP", "CP", "CP:BIOCARTA", 
      "CP:KEGG", "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS", 
      "MIR:MIRDB", "MIR:MIR_Legacy", "TFT:GTRD", "TFT:TFT_Legacy", 
      "CGN", "CM", "GO:BP", "GO:CC", "GO:MF", "HPO", "NULL", 
      "IMMUNESIGDB", "VAX", "NULL"), stringsAsFactors = F)
  gList <- x$gene_name
  if (!is.null(orgdb)) {
    for (gs in 1:nrow(MSigDb.categories)) {
      out.prefix <- gsub(" ", "", out.prefix)
      out.suffix <- paste("FDR_genes_MSigDb", MSigDb.categories[gs, 
        "category"], "Enrichment_Results", sep = "_")
      if (MSigDb.categories[gs, "sub.category"] == "NULL") {
        gDb <- msigdbr(species = msigdb_species, category = MSigDb.categories[gs, 
          "collection"]) %>% dplyr::select(gs_name, 
          gene_symbol) %>% as.data.frame()
      }
      else {
        gDb <- msigdbr(species = msigdb_species, category = MSigDb.categories[gs, 
          "collection"], subcategory = MSigDb.categories[gs, 
          "sub.category"]) %>% dplyr::select(gs_name, 
          gene_symbol) %>% as.data.frame()
      }
      enrich <- enricher(gene = gList, maxGSSize = 3000, 
        TERM2GENE = gDb)
      enrich_OUT <- as.data.frame(enrich)
      write.table(enrich_OUT, file = paste(out.prefix, 
        out.suffix, "tsv", sep = "."), quote = F, sep = "\t", 
        row.names = F)
    }
  }
  detach(package:msigdbr, unload = TRUE)
}


##
# Step 1: Check existing metadata
# View sample identities in the "orig.ident" column
unique(MPC2$orig.ident)

# Step 2: Add group information
# Assign multiple samples to a single group (e.g., "Group_A" includes Sample1, Sample2, Sample3)
MPC2$group <- ifelse(MPC2$orig.ident %in% c("Sample1", "Sample2", "Sample3"), 
                               "Group_A", 
                               "Other_Group")

# Step 3: Verify the grouping
table(MPC2$group)

# Optional: Add more groups if needed
MPC2$group[MPC2$orig.ident %in% c("Sample4", "Sample5")] <- "Group_B"

# Step 4: Use the group information in downstream analysis
# Example: Visualize UMAP by groups
DimPlot(MPC2, reduction = "umap", group.by = "group", label = TRUE)



# cells selection and deletion 
plot<-DimPlot(object = MPC2, reduction = "umap", label = TRUE)
plot$layyers
MPC2<- CellSelector(plot = plot, object =MPC2 , ident = "selected")
MPC2<-subset(MPC2, idents = 'selected', invert = TRUE)
levels(MPC2)
Idents(MPC2) <- "seurat_clusters"
DimPlot(object = MPC2, reduction = "umap", label = TRUE)

#Find all markers 
MPC2.markers <- FindAllMarkers(MPC2_finalisedcluster, only.pos = TRUE)
MPC2.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

#heatmap
MPC2.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(MPC2, features = top10$gene) + NoLegend()
setwd("C:/Users/maulik/Box/Gewin Lab Folder/Maulik/Ageing_MPC2/MPC2")
allmarkers <- FindAllMarkers(MPC2)
fwrite(MPC2.markers, "MPC2_allmarkers.csv")
write.csv(data, "C:/Users/maulik/Box/Gewin Lab Folder/Maulik/Ageing_MPC2/MPC2/'output_file.csv', raw.names = FALSE")

# pseudobulk cells only by cell type

bulk <- AggregateExpression(MPC2, group.by = "seurat_annotations", return.seurat = TRUE)
Cells(bulk)

# pseudobulk cells by stimulation condition AND cell type
bulk <- AggregateExpression(MPC2, group.by = c("group", "seurat_annotations"), return.seurat = TRUE)
Cells(bulk)

# pseudobulk cells by stimulation condition AND cell type AND donor
bulk <- AggregateExpression(MPC2, group.by = c("group", "seurat_annotations", "donor_id"), return.seurat = TRUE)
Cells(bulk)

pseudobulk_counts <- aggregate.Matrix(
  x = MPC2@assays$RNA@counts,   # Counts matrix
  groupings = MPC2$group,       # Grouping variable
  fun = "sum"                   # Summation function
)


#cell type name 
Mpc2_ls <- Mpc2_ls %>% RenameIdents(
  "0" = "CD_PC",
  "1" = "MTAL",
  "2" = "CD_ICA",
  "3" = "DTL",
  "4" = "PCT-S3",
  "5" = "PCT-S1 & S2",
  "6" = "Inj.EC",
  "7" = "MCPH",
  "8" = "PCT-S2",
  "9" = "PCT-S1",
  "10" = "B & T Cells",
  "11" = "EC",
  "12" = "MC",
  "13" = "CTAL",
  "14" = "DCT",
  "15" = "Inj.PT",
  "16" = "PODO"
  
)
levels <- c("PODO","MC","PCT-S1","PCT-S2","PCT-S3","PCT-S1&S2","DTL", "MTAL","CTAL", "DCT","CCD","CD-ICA","EC","Inj.EC","B&T Cell","MCPH","Inj.PT" ) # for cluster names
Idents(Mpc2_ls) <- factor(Idents(Mpc2_ls), levels=levels)
MPC2 <- MPC2 %>% RenameIdents("8" = "NK-myeloid")
Idents(Mpc2_ls)

new_idents <-c("0" = "Enc", "1"="PCT-S1","2" = "PCT-S2","3"= "CTAL","4"= "DCT","5"= "PCT-S3","6" ="MTAL","7"= "MC","8"= "NK-myeloid","9"= "InFl","10"= "CNT","11"= "CD-PC","12"= "DTL","13"= "CD-ICA","14"= "MCPHG","15"= "CD-ICB", "16"="T-cells","17"="B-cells","18"="Podo")
levels<-c("Enc", "PCT-S1","PCT-S2","CTAL","DCT","PCT-S3","MTAL","MC","NK-myeloid","InFl","CNT","CD-PC","DTL", "CD-ICA","MCPHG","CD-ICB", "T-cells","B-cells","Podo")
levels(Idents(MPC2)) <- new_idents
updated_idents <- Idents(MPC2)
print(Idents)
Idents(MPC2)<-new_idents


## 
seurat_annotations <- c("0" = "Enc", "1"="PCT-S1","2" = "PCT-S2","3"= "CTAL","4"= "DCT","5"= "PCT-S3","6" ="MTAL","7"= "MC","8"= "NK-myeloid","9"= "InFl","10"= "CNT","11"= "CD-PC","12"= "DTL","13"= "CD-ICA","14"= "MCPHG","15"= "CD-ICB", "16"="T-cells","17"="B-cells","18"="Podo")
MPC2_finalisedcluster$seurat_annotations <- seurat_annotations[as.character(Idents(MPC2_finalisedcluster))]

head(MPC2_finalisedcluster@meta.data)


####
Mpc2_ls <- NormalizeData(Mpc2_ls, normalization.method = "LogNormalize", scale.factor = 10000)
Mpc2_ls <- FindVariableFeatures(Mpc2_ls,selection.method = "RNA", nfeatures = 2000)
Mpc2_ls <- ScaleData(Mpc2_ls)
Mpc2_ls <- RunPCA(object = Mpc2_ls) 
Mpc2_ls <- RunUMAP(object= Mpc2_ls,dims = 1:13) 
Mpc2_ls <- FindNeighbors(object = Mpc2_ls, dims =1:13) 
Mpc2_ls <- FindClusters(object = Mpc2_ls, resolution = 0.3) 
DimPlot(object = Mpc2_ls, reduction = "umap", label = TRUE, raster = TRUE, alpha = 0.7)
#ggsave(file= "C:/Users/maulik/Box/Gewin Lab Folder/Maulik/Ageing_Mpc2_ls/Mpc2_ls/groups_col3_Dimplot_final.pdf")
DotPlot(Mpc2_ls,features=dotgenes,col.min = -0.1,col.max = 3,dot.min = 0,dot.scale = 10,scale.max = 100, cols =
          c("lightyellow","red"))+ scale_x_discrete(limit=rev)+ labs(x=NULL,
                                                                     y=NULL,title=print(paste("samples.merged")))+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=14,colour
                                 = "black"), axis.text.y=element_text(size=14,colour = "black"),
        legend.text = element_text(size=14,colour = "black"),
        panel.grid.minor=element_blank(), panel.grid.major=element_line(colour =
                                                                          "gray",linetype ="dashed",size=0.35))


## Sell selector
DimPlot(object = Mpc2_ls, reduction = "umap", label = TRUE, raster = TRUE, alpha = 0.7)
select.cells <- CellSelector(plot = DimPlot(Mpc2_ls))
Mpc2_ls <- subset(Mpc2_ls, cells = select.cells, invert = TRUE)
DimPlot(object = Mpc2_ls, reduction = "umap", label = TRUE, raster = TRUE, alpha = 0.7, split.by = "condition")
DotPlot(Mpc2_ls,features=dotgenes2,col.min = -0.1,col.max = 3,dot.min = 0,dot.scale = 10,scale.max = 100, cols =
          c("lightyellow","red"))+ scale_x_discrete(limit=rev)+ labs(x=NULL,
                                                                     y=NULL,title=print(paste("samples.merged")))+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=14,colour
                                 = "black"), axis.text.y=element_text(size=14,colour = "black"),
        legend.text = element_text(size=14,colour = "black"),
        panel.grid.minor=element_blank(), panel.grid.major=element_line(colour =
                                                                          "gray",linetype ="dashed",size=0.35))
#### Merge cluster
Mpc2_ls$seurat_clusters <- as.character(Idents(Mpc2_ls))  # Convert to character
Mpc2_ls$seurat_clusters[Mpc2_ls$seurat_clusters %in% c("2","8")] <- "2"
Idents(Mpc2_ls) <- Mpc2_ls$celltype


#### New cluster 


DimPlot(object = Mpc2_ls, reduction = "umap", label = TRUE, raster = TRUE, alpha = 0.7)
select.cells <- CellSelector(plot = DimPlot(Mpc2_ls))
# Step 2: Assign a new identity to the selected cells
Idents(Mpc2_ls)[select.cells] <- "NewCluster"
# Step 3: Visualize the updated identities
DimPlot(Mpc2_ls, reduction = "umap", label = TRUE, raster = TRUE, group.by = "ident")


## Rename the ident
Mpc2_ls <- RenameIdents(Mpc2_ls, "0" = "0", NewCluster" = "14")


##New cluster numbering 

# Step 1: Get current cluster identities
current_clusters <- levels(Idents(Mpc2_ls))

# Step 2: Create a new numeric mapping
# This will assign new numbers starting from 0
new_cluster_ids <- setNames(as.character(seq_along(current_clusters) - 1), current_clusters)

# Step 3: Apply the new cluster numbering
Mpc2_ls <- RenameIdents(Mpc2_ls, new_cluster_ids)

# Step 4: Visualize the updated clusters
DimPlot(Mpc2_ls, reduction = "umap", label = TRUE, raster = TRUE)



featu####

# Step 1: Extract current identities
idents_vec <- Idents(Mpc2_ls)

# Step 2: Find cells with NA identity
na_cells <- names(idents_vec[is.na(idents_vec)])

# Step 3: Convert identities to character vector
idents_char <- as.character(idents_vec)

# Step 4: Assign a new cluster label to NA cells
idents_char[na_cells] <- "14" # or any label you prefer

# Step 5: Set the updated identities back to the Seurat object
Idents(Mpc2_ls) <- factor(idents_char)

# Step 6: Check the result
table(Idents(Mpc2_ls))


#### Find Markers of the clusters
cluster17.markers <- FindMarkers(Mpc2_ls, ident.1 = 17)
head(cluster17.markers, n = 5)

###Find all markerss
Mpc2_ls.markers <- FindAllMarkers(Mpc2_ls, only.pos = TRUE)
(Mpc2_ls.markers %>%
group_by(cluster) %>%
)dplyr::filter(avg_log2FC > 1) %>%
slice_head(n = 5) %>%
ungroup() -> top5
DoHeatmap(Mpc2_ls, features = top5$gene) + NoLegend()

## to add NA to new clusters 
# Convert Idents to a data frame
df <- data.frame(Idents = Idents(Mpc2_ls))

# Add "NewCluster" to the levels using forcats
df$Idents <- fct_expand(df$Idents, "NewCluster")

# Replace NA values with "NewCluster"
df <- df %>% mutate(Idents = ifelse(is.na(Idents), "NewCluster", as.character(Idents)))

# Convert back to factor and assign to Idents(Mpc2_ls)
Idents(Mpc2_ls) <- as.factor(df$Idents)
Idents(Mpc2_ls)

