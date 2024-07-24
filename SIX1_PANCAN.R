setwd("C:/Users/76064/Desktop/work/projrct1")
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)

annotation <- fread("probeMap_gencode.v23.annotation.gene.probemap")
# 读取表达矩阵
expression <- fread("TcgaTargetGtex_rsem_gene_tpm")
# 读取表型数据
phenotype <- fread("TcgaTargetGTEX_phenotype.txt")
# 提取SIX1基因ID
six1_gene_id <- annotation %>%
  filter(gene == "SIX1") %>%
  select(id) %>%
  unique() %>%
  pull()
# 提取SIX1的表达数据
six1_expression <- expression %>%
  filter(sample == six1_gene_id)
# 转置表达数据
six1_expression_transposed <- t(six1_expression[, -1, with = FALSE])
colnames(six1_expression_transposed) <- six1_expression$sample
six1_expression_transposed <- as.data.frame(six1_expression_transposed)
six1_expression_transposed$sample_id <- rownames(six1_expression_transposed)
# 合并表型数据
combined_data <- merge(six1_expression_transposed, phenotype,
                       by.x = "sample_id", by.y = "sample")
# 去除ICC样本
combined_data <- combined_data %>%
  filter(detailed_category != "Cholangiocarcinoma")
# 重新分类样本类型
combined_data <- combined_data %>%
  mutate(sample_type = case_when(
    `_sample_type` %in% c("Primary Solid Tumor", "Recurrent Solid Tumor", "Metastatic", "Primary Tumor", "Recurrent Tumor", "Additional - New Primary", "Additional Metastatic") ~ "Tumor",
    `_sample_type` %in% c("Normal Tissue", "Solid Tissue Normal") ~ "Normal",
    TRUE ~ NA_character_
  ))
# 移除 NA 值的样本类型
combined_data <- combined_data %>%
  filter(!is.na(sample_type))

# 对每种癌症类型分别进行Wilcoxon检验
cancer_types <- unique(combined_data$detailed_category)
results <- data.frame(Cancer_Type = character(), P_Value = numeric(), Tumor_Count = integer(), Normal_Count = integer(), stringsAsFactors = FALSE)

for (cancer in cancer_types) {
  cancer_data <- combined_data %>%
    filter(detailed_category == cancer)
  
  tumor_count <- sum(cancer_data$sample_type == "Tumor")
  normal_count <- sum(cancer_data$sample_type == "Normal")
  
  if (tumor_count > 0 & normal_count > 0) { # Ensure both Normal and Tumor samples are present
    wilcox_test <- wilcox.test(as.numeric(cancer_data[[six1_gene_id]]) ~ cancer_data$sample_type)
    p_value <- wilcox_test$p.value
    results <- rbind(results, data.frame(Cancer_Type = cancer, P_Value = p_value, Tumor_Count = tumor_count, Normal_Count = normal_count))
  }
}

# 标记显著性
results$Significance <- cut(results$P_Value, breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("****", "***", "**", "*", "ns"))

# 准备绘图数据
combined_data <- combined_data %>%
  mutate(Cancer_Type = factor(detailed_category, levels = results$Cancer_Type))
# 检查并移除含有 NA 的行
combined_data <- combined_data %>%
  filter(!is.na(Cancer_Type))
# 确定每种癌症类型的最大表达值，以便放置显著性标记
max_expression <- combined_data %>%
  group_by(Cancer_Type) %>%
  summarise(Max_Expression = max(as.numeric(get(six1_gene_id)), na.rm = TRUE))

results <- merge(results, max_expression, by = "Cancer_Type")

# 创建带有样本计数的癌症类型标签
results <- results %>%
  mutate(Cancer_Type_Label = paste0(Cancer_Type, " (T=", Tumor_Count, ", N=", Normal_Count, ")"))

combined_data <- combined_data %>%
  left_join(results, by = c("Cancer_Type" = "Cancer_Type")) %>%
  mutate(Cancer_Type_Label = factor(Cancer_Type_Label, levels = results$Cancer_Type_Label))


# 绘制小提琴图
ggplot(combined_data, aes(x = Cancer_Type_Label, y = as.numeric(get(six1_gene_id)), fill = sample_type)) +
  geom_violin(trim = FALSE, scale = "width", color = "black", alpha = 0.7) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), color = "black", outlier.size = 0.5) +
  theme_minimal(base_size = 15) +
  labs(title = "SIX1 Expression in Normal vs Tumor Samples Across Cancer Types",
       x = "Cancer Type (Tumor Count, Normal Count)",
       y = "Expression",
       fill = "Sample Type") +
  scale_fill_manual(values = c("Tumor" = "#E69F00", "Normal" = "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "top",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  geom_text(data = results, aes(x = Cancer_Type_Label, y = Max_Expression + 1, label = Significance), inherit.aes = FALSE, size = 5, vjust = -0.5, color = "black", fontface = "bold") +
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1))




