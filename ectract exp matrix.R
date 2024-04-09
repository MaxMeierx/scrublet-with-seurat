#save the expression matrix 
#seurat_list is a list which contains the seurat object
#finding  double should before merge

#expression matrix filefold
fold='./data/表达矩阵/'
for (i in names(seurat_list)) {
  # 提取每个 Seurat 对象的表达矩阵
  # 默认情况下，GetAssayData函数返回的是counts数据
  expr_matrix <- GetAssayData(seurat_list[[i]], slot = "counts")
  
  # 构造一个文件名，使用 Seurat 对象的名字
  filename <- paste0(fold,"expression_matrix_", i, ".csv")
  
  # 保存为CSV文件
  write.csv(expr_matrix, filename, row.names = TRUE)
}
