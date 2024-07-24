setwd("C:/Users/76064/Desktop/work/projrct1/SIX1_ICC")
library(TCGAbiolinks)
library(dplyr)
library(survival)
library(survminer)
library(SummarizedExperiment)
library(biomaRt)
library(ggplot2)
query <- GDCquery(project = "TCGA-CHOL", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")
GDCdownload(query)
data <- GDCprepare(query)
# 提取SIX1基因表达数据
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 获取SIX1基因的Ensembl ID
gene_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                   filters = 'hgnc_symbol', 
                   values = 'SIX1', 
                   mart = ensembl)

ensembl_id <- gene_info$ensembl_gene_id
print(ensembl_id)
# 提取SIX1基因表达数据
SIX1_expr <- assay(data, "tpm_unstrand")[sub("\\..*", "", rownames(assay(data))) %in% ensembl_id, ]

group <- ifelse(colData(data)$sample_type == "Primary Tumor", "Tumor", "Normal")

# 创建数据框用于可视化
expr_data <- data.frame(
  Expression = SIX1_expr,
  Group = group
)

# Wilcoxon检验
wilcox.test(Expression ~ Group, data = expr_data)

# 绘制箱线图，使用log缩放
ggplot(expr_data, aes(x = Group, y = Expression, fill = Group)) + 
  geom_boxplot() +
  scale_y_log10() +
  labs(title = "SIX1 Expression in Normal vs Tumor Samples (Log Scale)", 
       x = "Subtype", 
       y = "Log(Expression)") +
  scale_fill_manual(values = c("Normal" = "blue", "Tumor" = "pink"))

# 生存分析

# 计算SIX1表达的中值
median_SIX1 <- median(SIX1_expr[group == "Tumor"])

# 分组
SIX1_group <- ifelse(SIX1_expr > median_SIX1, "High", "Low")

# 提取生存时间和生存状态信息
surv_time <- colData(data)$days_to_death
surv_status <- ifelse(colData(data)$vital_status == "Dead", 1, 0)

# 创建生存对象
surv_object <- Surv(surv_time, surv_status)
# 进行生存分析
fit <- survfit(surv_object ~ SIX1_group)
# 绘制Kaplan-Meier曲线
ggsurvplot(fit, data = data.frame(SIX1_group, surv_time, surv_status), 
           pval = TRUE, 
           risk.table = TRUE, 
           conf.int = TRUE,
           palette = c("#E7B800", "#2E9FDF"),
           legend.title = "SIX1 Expression",
           legend.labs = c("High-expression", "Low-expression"),
           title = "Kaplan-Meier Curve for SIX1 Expression",
           xlab = "Time (months)",
           ylab = "Survival Probability")
# 使用log-rank检验
surv_diff <- survdiff(surv_object ~ SIX1_group)
print(surv_diff)
