# GO富集分析脚本

# 加载必要的包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("clusterProfiler", quietly = TRUE))
    BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))  # 人类基因注释包，根据物种修改
    BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("enrichplot", quietly = TRUE))
    BiocManager::install("enrichplot")
if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr")

library(clusterProfiler)
library(org.Hs.eg.db)  # 人类基因注释包，根据研究物种替换
library(dplyr)

# 设置输入输出路径
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("使用方法: Rscript go_enrichment.R <基因列表文件> <输出结果文件>\n")
  cat("示例: Rscript go_enrichment.R genes.txt go_results.csv\n")
  quit(save = "no", status = 1)
}

gene_list_file <- args[1]
output_file <- args[2]

# 读取基因列表
genes <- readLines(gene_list_file)
genes <- unique(genes)  # 确保基因ID唯一

# 过滤掉无效基因（可选步骤）
gene_map <- bitr(genes, fromType = "SYMBOL",  # 根据输入基因ID类型修改fromType
                 toType = c("ENTREZID", "ENSEMBL"), 
                 OrgDb = org.Hs.eg.db)  # 根据物种修改OrgDb

# 检查是否有映射的基因
if (nrow(gene_map) == 0) {
  stop("没有找到有效的基因ID映射，请检查输入基因ID类型和物种注释包是否匹配。")
}

# 使用ENTREZ ID进行富集分析
mapped_genes <- gene_map$ENTREZID

# 进行GO富集分析 (BP: 生物过程, MF: 分子功能, CC: 细胞组分)
ego_BP <- enrichGO(gene          = mapped_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "BP",  # 可选 "BP", "MF", "CC", 或 "ALL"
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

ego_MF <- enrichGO(gene          = mapped_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

ego_CC <- enrichGO(gene          = mapped_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

# 合并所有GO分类的结果
if (!is.null(ego_BP)) {
  ego_BP$ONTOLOGY <- "BP"
}
if (!is.null(ego_MF)) {
  ego_MF$ONTOLOGY <- "MF"
}
if (!is.null(ego_CC)) {
  ego_CC$ONTOLOGY <- "CC"
}

all_results <- rbind(as.data.frame(ego_BP), 
                    as.data.frame(ego_MF),
                    as.data.frame(ego_CC))

# 如果没有富集结果，创建空表格
if (nrow(all_results) == 0) {
  all_results <- data.frame(
    ID = character(),
    Description = character(),
    ONTOLOGY = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    pAdjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer()
  )
}

# 保存结果到CSV文件
write.csv(all_results, file = output_file, row.names = FALSE)

cat(paste("GO富集分析完成，结果已保存至:", output_file, "\n"))