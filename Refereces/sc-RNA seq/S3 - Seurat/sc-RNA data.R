#run libraries
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(patchwork)
library(ggplot2)
#data is from: https://www.ncbi.nlm.nih.gov/geo/ , GSE132771


#reading files : each file represent one sample
nml1 <- Read10X(data.dir = "/Users/mozhganoroujlu/Desktop/MOZHGUN/Laleh/Data/real data/GSE132771/GSE132771_RAW/NLM1")
nml2 <- Read10X(data.dir = "/Users/mozhganoroujlu/Desktop/MOZHGUN/Laleh/Data/real data/GSE132771/GSE132771_RAW/NLM2")
nml3 <- Read10X(data.dir = "/Users/mozhganoroujlu/Desktop/MOZHGUN/Laleh/Data/real data/GSE132771/GSE132771_RAW/NLM3")

#view raw data
nml1
nml2
nml3

#creating Seurat object  
nml1_obj <- CreateSeuratObject(counts = nml1, project = "nml1",
                               min.cells = 3, min.features = 200) #min.cells: filters out genes that are detected in fewer than 3 cells.
                                                                  #min.features: filters out cells that have fewer than 200 detected genes
nml2_obj <- CreateSeuratObject(counts = nml2, project = "nml2",
                               min.cells = 3, min.features = 200) 


nml3_obj <- CreateSeuratObject(counts = nml3, project = "nml3",
                               min.cells = 3, min.features = 200) 

#view Seurat objects
view(nml1_obj@meta.data)
view(nml2_obj@meta.data)
view(nml3_obj@meta.data)

#saving Seurat objects in computer
saveRDS(nml1_obj, file = "/Users/mozhganoroujlu/Desktop/MOZHGUN/Laleh/Data/real data/nml1_obj.rds")
saveRDS(nml2_obj, file = "/Users/mozhganoroujlu/Desktop/MOZHGUN/Laleh/Data/real data/nml2_obj.rds")
saveRDS(nml3_obj, file = "/Users/mozhganoroujlu/Desktop/MOZHGUN/Laleh/Data/real data/nml3_obj.rds")
                               
#reading RDS files in R
sample_1 <- readRDS(file = "/Users/mozhganoroujlu/Desktop/MOZHGUN/Laleh/Data/real data/GSE132771/Seurat Objects/nml1_obj.rds")
sample_2 <- readRDS(file = "/Users/mozhganoroujlu/Desktop/MOZHGUN/Laleh/Data/real data/GSE132771/Seurat Objects/nml2_obj.rds")
sample_3 <- readRDS(file = "/Users/mozhganoroujlu/Desktop/MOZHGUN/Laleh/Data/real data/GSE132771/Seurat Objects/nml3_obj.rds")

#merge 3 samples in one RSD file
mergedSamples <- merge(x = sample_1, y = c(sample_2, sample_3), 
                       add.cell.ids = c("1__","2__","3__"), project = "MergedSamples")

#view merged sample
view(mergedSamples) # shows first 1000 rows

#view all of the merged samples
view(mergedSamples@meta.data) #meta.data shows all of the merged samples

#quality control

range(mergedSamples$nFeature_RNA) # shows number of unique genes detected per cell.
range(mergedSamples$nCount_RNA) # shows the total number of RNA molecules detected.

#mitochondrial genes
mergedSamples_MTgenes <- PercentageFeatureSet(mergedSamples, pattern = "MT-", 
                                      col.name = "percentage.mt.genes" ) #add percentage of mitochondrial genes to merged samples

#show the range of % mitochondrial genes per cell
range(mergedSamples_MTgenes$percentage.mt.genes)

#plot distribution of gene expression
VlnPlot(mergedSamples_MTgenes, features = c("nFeature_RNA", "nCount_RNA", "percentage.mt.genes"), ncol = 3)

#filter abnormal cells (quality control)
CorrectedMergedSamples <- subset(mergedSamples_MTgenes, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 &
                                   nCount_RNA<23000 & percentage.mt.genes <10 )

#show new range of features, counts and MT genes
range(CorrectedMergedSamples$nFeature_RNA) 
range(CorrectedMergedSamples$nCount_RNA) 
range(CorrectedMergedSamples$percentage.mt.genes)

#data normalization
NormalizedData <- NormalizeData(CorrectedMergedSamples, 
                                        
                                        normalization.method ="LogNormalize", scale.factor = 10000) # Log-normalized value(ij)=ln(/X(ij)∑iX(ij)/Xij × scale factor+1)

#Find variable features
 NormalizedData <- FindVariableFeatures(NormalizedData, selection.method = "vst") #identify genes whose expression varies significantly across cells. 
 
# Scale data
NormalizedData <- ScaleData(NormalizedData)

# Extract the mean, variance, and fitted variance from the Seurat object
hvf_info <- HVFInfo(NormalizedData, selection.method = "vst")

# Create the plot for finding variables
ggplot(hvf_info, aes(x = mean, y = variance)) +
  geom_point(aes(color = "Non-variable"), size = 1.5, alpha = 0.6) +  # Plot all genes
  geom_line(aes(y = variance.expected, color = "Fitted curve"), size = 1) +  # Add fitted curve
  geom_point(
    data = hvf_info[VariableFeatures(NormalizedData), ],  # Highlight variable genes
    aes(color = "Variable"), size = 1.5, alpha = 0.8
  ) +
  scale_x_log10() +  # Log-transform x-axis
  scale_y_log10() +  # Log-transform y-axis
  scale_color_manual(
    values = c("Non-variable" = "black", "Variable" = "red", "Fitted curve" = "blue")
  ) +
  labs(
    x = "Mean expression (log10)",
    y = "Variance (log10)",
    title = "Mean-Variance Plot with Fitted Curve",
    color = "Legend"
  ) +
  theme_minimal()

# Perform PCA
FinalData <- RunPCA(NormalizedData) #reduce dimension to 10 dim

# Cluster cells
FinalData <- FindNeighbors(FinalData, dims = 1:10) #KNN
FinalData <- FindClusters(FinalData, resolution = 0.5) # default algorithm : Louvain

# Run UMAP for visualization
FinalData <- RunUMAP(FinalData, dims = 1:10) #specifies that the first 10 principal components should be used as input for UMAP.

# Visualize clusters
DimPlot(FinalData, reduction = "umap")

# save merged sample
saveRDS(FinalData, file = "/Users/mozhganoroujlu/Desktop/MOZHGUN/Laleh/Data/real data/GSE132771/Seurat Objects/Final Data.rds")


                               