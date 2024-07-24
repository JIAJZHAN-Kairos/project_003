setwd("C:/Users/76064/Desktop/work/projrct1/DEGs")
library(TCGAbiolinks)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(pheatmap)

query <- GDCquery(project = "TCGA-CHOL", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification")

GDCdownload(query)
data <- GDCprepare(query)

# 提取表达计数矩阵
countData <- assay(data)

# 提取样本信息
colData <- colData(data)
# 根据样本类型添加分组信息
colData$group <- ifelse(colData$sample_type == "Primary Tumor", "high_expression", "low_expression")
# 创建DESeq2对象
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ group)
# 过滤低表达基因
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# 数据标准化
vsd <- vst(dds, blind = FALSE)

# PCA分析并绘图
pcaData <- plotPCA(vsd, intgroup = c("group", "shortLetterCode"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaPlot <- ggplot(pcaData, aes(PC1, PC2, color = group, shape = shortLetterCode)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  geom_text_repel(aes(label = rownames(pcaData)))

print(pcaPlot)

# 提取标准化后的表达数据
normCounts <- assay(vsd)

# 将数据转为长格式
normCounts_long <- melt(normCounts)
colnames(normCounts_long) <- c("Gene", "Sample", "Expression")

# 添加批次信息
normCounts_long$Batch <- colData$shortLetterCode[match(normCounts_long$Sample, rownames(colData))]

# 绘制箱线图
boxplot <- ggplot(normCounts_long, aes(x = Batch, y = Expression, fill = Batch)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplot of Normalized Gene Expression by Batch",
       x = "Batch",
       y = "Normalized Expression")

print(boxplot)

# 识别离群样本（根据PCA图观察手动标记离群样本）
outliers <- c("TCGA-3X-AAV9-01A-72R-A41I-07", "TCGA-W5-AA2Q-11A-11R-A41I-07", "TCGA-W5-AA2H-01A-31R-A41I-07")

# 从countData和colData中剔除离群样本
countData_clean <- countData[, !(colnames(countData) %in% outliers)]
colData_clean <- colData[!(rownames(colData) %in% outliers),]

# 创建新的DESeq2对象
dds_clean <- DESeqDataSetFromMatrix(countData = countData_clean, colData = colData_clean, design = ~ group)

# 数据过滤
keep <- rowSums(counts(dds_clean) >= 10) >= 3
dds_clean <- dds_clean[keep,]

# 数据标准化
vsd_clean <- vst(dds_clean, blind = FALSE)

# PCA分析并绘图（剔除离群样本后）
pcaData_clean <- plotPCA(vsd_clean, intgroup = c("group", "shortLetterCode"), returnData = TRUE)
percentVar_clean <- round(100 * attr(pcaData_clean, "percentVar"))

pcaPlot_clean <- ggplot(pcaData_clean, aes(PC1, PC2, color = group, shape = shortLetterCode)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_clean[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_clean[2], "% variance")) +
  theme_minimal() +
  geom_text_repel(aes(label = rownames(pcaData_clean)))

print(pcaPlot_clean)
# 差异表达分析
dds_clean <- DESeq(dds_clean)
res <- results(dds_clean)

# 筛选差异表达基因 |log2(FC)| > 1 and adj.P < 0.05
DEGs <- subset(res)

# 检查DEGs
head(DEGs)


library(biomaRt)
# 使用biomaRt连接Ensembl数据库
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# 提取Ensembl基因ID
ensembl_ids <- sub("\\..*", "", rownames(DEGs))

# 获取基因名称
genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
               filters = "ensembl_gene_id", 
               values = ensembl_ids, 
               mart = mart)

# 将结果转换为数据框
genes_df <- as.data.frame(genes)
head(genes_df)

# 确保DEGs为数据框
DEGs_df <- as.data.frame(DEGs)
DEGs_df$ensembl_gene_id <- sub("\\..*", "", rownames(DEGs_df))

# 合并数据框
DEGs_named <- merge(DEGs_df, genes_df, by = "ensembl_gene_id")
# 去除没有基因名称的行
DEGs_named <- DEGs_named[DEGs_named$external_gene_name != "", ]
# 添加唯一标识到重复基因名称
unique_gene_names <- make.unique(DEGs_named$external_gene_name)

# 使用基因名称作为行名
rownames(DEGs_named) <- unique_gene_names
DEGs_named <- DEGs_named[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

# 检查结果
head(DEGs_named)
# 添加分组信息
DEGs_named$group <- "NOT"
DEGs_named$group[DEGs_named$log2FoldChange > 1 & DEGs_named$padj < 0.05] <- "UP"
DEGs_named$group[DEGs_named$log2FoldChange < -1 & DEGs_named$padj < 0.05] <- "DOWN"

# 选择需要标注的基因
top_genes <- rownames(DEGs_named)[order(DEGs_named$padj)][1:15]  # 根据需要调整选择条件
DEGs_named$label <- ifelse(rownames(DEGs_named) %in% top_genes, rownames(DEGs_named), "")

# 火山图
volcano <- ggplot(DEGs_named, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = group), alpha = 0.6, size = 1.5) +  # 调整alpha和size
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_text_repel(aes(label = label), size = 3) +
  theme_minimal() +
  labs(title = "Volcano Plot of DEGs", x = "log2 Fold Change", y = "-log10 Adjusted P-value") +
  theme(legend.position = "top", legend.title = element_blank())

print(volcano)
# 检查极端值
extreme_values <- DEGs_named[DEGs_named$log2FoldChange < -10 | DEGs_named$log2FoldChange > 10, ]
print(extreme_values)
# 移除极端值
DEGs_named <- DEGs_named[DEGs_named$log2FoldChange >= -10 & DEGs_named$log2FoldChange <= 10, ]

# 使用清理后的数据绘制火山图
volcano <- ggplot(DEGs_named, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = group), alpha = 0.6, size = 1.5) +  # 调整alpha和size
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_text_repel(aes(label = label), size = 3) +
  theme_minimal() +
  labs(title = "Volcano Plot of DEGs", x = "log2 Fold Change", y = "-log10 Adjusted P-value") +
  theme(legend.position = "top", legend.title = element_blank())

print(volcano)


# 筛选和排序上调和下调的TOP 10-15基因
top_up_genes <- rownames(DEGs_named[DEGs_named$log2FoldChange > 0, ][order(DEGs_named$padj[DEGs_named$log2FoldChange > 0]), ])[1:15]
top_down_genes <- rownames(DEGs_named[DEGs_named$log2FoldChange < 0, ][order(DEGs_named$padj[DEGs_named$log2FoldChange < 0]), ])[1:15]
top_genes <- c(top_up_genes, top_down_genes)

# 标注TOP基因
DEGs_named$label <- ifelse(rownames(DEGs_named) %in% top_genes, rownames(DEGs_named), "")

# 提取vsd_clean中的基因ID并去除版本号
vsd_ensembl_ids <- rownames(assay(vsd_clean))
vsd_ensembl_ids_no_version <- sub("\\..*", "", vsd_ensembl_ids)

# 获取基因名称
vsd_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                   filters = "ensembl_gene_id", 
                   values = vsd_ensembl_ids_no_version, 
                   mart = mart)

# 将结果转换为数据框
vsd_genes_df <- as.data.frame(vsd_genes)

# 确保vsd_clean为数据框
vsd_counts <- as.data.frame(assay(vsd_clean))
vsd_counts$ensembl_gene_id <- sub("\\..*", "", rownames(vsd_counts))

# 合并数据框
vsd_counts_named <- merge(vsd_counts, vsd_genes_df, by = "ensembl_gene_id", all.x = TRUE)

# 检查缺失值并移除缺失值行
vsd_counts_named <- vsd_counts_named[!is.na(vsd_counts_named$external_gene_name), ]
# 去除没有基因名称的行
vsd_counts_named <- vsd_counts_named[vsd_counts_named$external_gene_name != "", ]

# 添加唯一标识到重复基因名称
unique_vsd_gene_names <- make.unique(vsd_counts_named$external_gene_name)

# 使用基因名称作为行名
rownames(vsd_counts_named) <- unique_vsd_gene_names


# 移除ensembl_gene_id列并转换为矩阵
vsd_counts_named_matrix <- as.matrix(vsd_counts_named[, -ncol(vsd_counts_named)])


# 提取差异表达基因表达值
DEG_counts <- vsd_counts_named_matrix[top_genes, ]

# 确保 DEG_counts 是一个矩阵，并且行名和列名正确
if (!is.matrix(DEG_counts)) {
  DEG_counts <- as.matrix(DEG_counts)
}

# 检查 colData_clean 的列名和 DEG_counts 的列名是否匹配
colData_clean <- as.data.frame(colData_clean)  # 确保 colData_clean 是一个数据框
if (!all(colnames(DEG_counts) %in% rownames(colData_clean))) {
  stop("列名不匹配，请检查 colData_clean 和 DEG_counts 的列名。")
}

# 提取 colData_clean 中匹配的列
annotation_col <- colData_clean[colnames(DEG_counts), , drop = FALSE]

# 检查 annotation_col 是否有多维数组问题
if (dim(annotation_col)[2] > 1) {
  annotation_col <- annotation_col[, 1, drop = FALSE]
}

# 绘制热图
pheatmap(DEG_counts, cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = annotation_col, show_rownames = TRUE, show_colnames = FALSE,
         main = "Heatmap of Top DEGs")

write.csv(DEGs_named, "DEGs_results.csv", row.names = TRUE)
