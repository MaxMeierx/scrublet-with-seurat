#ensure seurat_list is in the environment

files <- list.files(path = "result fold", pattern = "\\.txt$", full.names = TRUE)



for (file in files) {
  # 提取样本名
  file_name <- basename(file)
  sample_name <- gsub("expression_matrix_(.*)_doublet.txt", "\\1", file_name )
  
  # 读取文件内容
  doublet_data <- read_csv(file)
  doublet_data <- column_to_rownames(doublet_data ,colnames(doublet_data)[1])
  rownames(doublet_data) <- gsub("^X", "", rownames(doublet_data))
  rownames(doublet_data) <- gsub("\\.", "-", rownames(doublet_data))
  # 找到对应的Seurat对象
  if (sample_name %in% names(seurat_list)) {
    seurat_object <- seurat_list[[sample_name]]
    meta=seurat_object@meta.data
    meta=merge(x=meta,y=doublet_data,by.x=0,by.y=0)
    meta=column_to_rownames(meta,colnames(meta)[1])
    # 将更新后的Seurat对象保存回列表
    seurat_object@meta.data=meta
    seurat_list[[sample_name]] <- seurat_object
  } else {
    cat("Sample name", sample_name, "not found in seurat_list.\n")
  }
}

for (i in names(seurat_list)) {
  seurat_object <- seurat_list[[i]]
  seurat_object <- subset(seurat_object ,subset=predicted_doublets=='FALSE')
  seurat_list[[i]]=seurat_object
}
