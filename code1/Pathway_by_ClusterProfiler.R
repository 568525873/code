# KEGG富集分析脚本

# 加载必要的包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("clusterProfiler", quietly = TRUE))
    BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))  # 人类基因注释包，根据物种修改
    BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr")
if (!requireNamespace("pathview", quietly = TRUE))
    BiocManager::install("pathview")

library(clusterProfiler)
library(org.Hs.eg.db)  # 人类基因注释包，根据研究物种替换
library(dplyr)

# 设置输入输出路径
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("使用方法: Rscript kegg_enrichment.R <基因列表文件> <输出结果文件>\n")
  cat("示例: Rscript kegg_enrichment.R genes.txt kegg_results.csv\n")
  quit(save = "no", status = 1)
}

gene_list_file <- args[1]
output_file <- args[2]

# 读取基因列表
genes <- readLines(gene_list_file)
genes <- unique(genes)  # 确保基因ID唯一

# 基因ID映射为ENTREZID
gene_map <- bitr(genes, fromType = "SYMBOL",  # 根据输入基因ID类型修改fromType
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)  # 根据物种修改OrgDb

# 检查是否有映射的基因
if (nrow(gene_map) == 0) {
  stop("没有找到有效的基因ID映射，请检查输入基因ID类型和物种注释包是否匹配。")
}

# 使用ENTREZ ID进行富集分析
mapped_genes <- gene_map$ENTREZID

# 进行KEGG富集分析
kk <- enrichKEGG(gene         = mapped_genes,
                 organism     = 'hsa',  # 根据物种修改，如小鼠为'mmu'
                 keyType      = 'kegg',
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)  # 宽松的q值阈值，可根据需要调整

# 转换KEGG ID为可读格式
if (!is.null(kk)) {
  kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
} else {
  # 如果没有富集结果，创建空表格
  kk <- data.frame(
    ID = character(),
    Description = character(),
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
write.csv(as.data.frame(kk), file = output_file, row.names = FALSE)

cat(paste("KEGG富集分析完成，结果已保存至:", output_file, "\n"))