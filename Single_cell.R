suppressMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(Matrix)
  library(valiDrops)
  library(future)
  library(NMF)
  library(ggalluvial)
  library(loupeR)
  library(openxlsx)
})

options(future.globals.maxSize = 60000 * 1024^2) 
plan(multisession,workers=4)
options(future.seed = TRUE)
setwd("/home/genomics/single_cell/Raw from SB/")
# 定義樣本列表和分組資訊
sample.table <- read.table("./sample.list",header = T)

## -------- 建立樣本input list -------------------------------------------------
sample <- list()
sample <- purrr::pmap(sample.table,
                      ~list(file = ..1, id = ..2, group = ..3))

outdir  <- file.path(getwd(),"single_cell/restult_20251009/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
summary.table <- NULL
process_sample <- function(sample) {
  
  st_file <- file.path(paste0(sample$file,"/",sample$file,"_Expression_Data.st"))
  message("[Loading] ", st_file)
 
    ## ---------- 讀取 H5 ------------------------------------------------------
    st_table <- read.table(st_file, header = TRUE)
    st_matrix <- sparseMatrix(i=as.numeric(factor(st_table$Bioproduct)),
                              j=as.numeric(factor(st_table$Cell_Index)),
                              x=st_table$RSEC_Adjusted_Molecules)
    colnames(st_matrix) <- levels(factor(st_table$Cell_Index))
    rownames(st_matrix) <- levels(factor(st_table$Bioproduct))
    origin_count <- length(st_matrix@Dimnames[[2]])
    # ---------- 建立 raw 的 seurat 物件 -------------------------------------
     seurat.obj_raw <- CreateSeuratObject(
       st_matrix,
       project      =  sample$id,
       min.cells    = 1,
       min.features = 1
     )
     seurat.obj_raw[["percent.mt"]] <- PercentageFeatureSet(seurat.obj_raw, pattern = "^mt-")
     seurat.obj_raw <- AddMetaData(seurat.obj_raw, metadata = sample$group, col.name = "group")
     png(paste0(outdir,sample$id,"_sample_vlnplot_valid_raw.png"),width=4000,height = 2000,res=300)
     p <- VlnPlot(seurat.obj_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                  fill.by = orig.ident)
     print(p)
     dev.off()
    ## ---------- mt < 15 % 篩選 ----------------------------------------------
     seurat_obj <- subset(seurat.obj_raw, subset = nFeature_RNA > 200 & percent.mt < 15)
     
  ## ---------- 標準流程 ------------------------------------------------------
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
    seurat_obj <- FindClusters(seurat_obj)
    origin_count <- length(st_matrix@Dimnames[[2]])
    ## ---------- Doublet 檢測 --------------------------------------------------
    
    sce <- as.SingleCellExperiment(seurat_obj)
    sce <- scDblFinder(sce, cluster = TRUE)
    
    seurat_obj$doublet <- as.factor(sce$scDblFinder.class)
    Idents(seurat_obj)  <- seurat_obj$doublet
    singlet_cells       <- WhichCells(seurat_obj, idents = "singlet")
    seurat_obj          <- subset(seurat_obj, cells = singlet_cells)
    
    ## ---------- 追蹤 singlet / doublet 統計 ----------------------------------
    sample.qc <- data.frame(sample  = sample$id,
                            origin = origin_count,
                            singlet = sum(sce$scDblFinder.class == "singlet"),
                            doublet = sum(sce$scDblFinder.class == "doublet"))
    
    
    message("[Done] ", sample$id,
            " | singlet: ", sample.qc$singlet,
            " | doublet: ", sample.qc$doublet)
    
    sample.obj <- list( seurat_obj_raw = seurat.obj_raw, seurat_obj = seurat_obj, summary.table = sample.qc)
  
    return(sample.obj)
    
  }


## --------- 執行 --------------------------------------------------------------
library(future); library(future.apply)
seurat_objects <- lapply(sample, process_sample)

sample.qc <- data.frame(sample  = as.character(),
                        origin = as.numeric(),
                        singlet = as.numeric(),
                        doublet = as.numeric())
seurat_obj <- list()
seurat_raw <- list()
for (s in 1:length(seurat_objects)){
  seurat_obj[[s]] <- seurat_objects[[s]]$seurat_obj
  seurat_raw[[s]] <- seurat_objects[[s]]$seurat_obj_raw
  sample.qc <- rbind(sample.qc, seurat_objects[[s]]$summary.table)
}



write.xlsx(sample.qc ,"./sample_qc.xlsx")
saveRDS(seurat_obj,"./seurat_object.rds")
saveRDS(seurat_raw,"./seurat_raw.rds")


