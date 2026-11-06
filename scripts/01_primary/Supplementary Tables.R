# =========================================================================
# 自动化打包生成 Supplementary Tables.xlsx
# =========================================================================

# 1. 检测并安装所需的包 (openxlsx 专用于完美输出多Sheet的Excel文件)
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
library(openxlsx)

# 2. 设置你的工作目录 (如果已设置，此行可跳过；如果未设置，请取消注释并修改路径)
# setwd(Sys.getenv("WB_DATA_ROOT", unset = "data/raw"))

print("正在读取原始 CSV 文件...")

# 3. 根据你的文件名读取数据，并映射到新的 Table S1-S4 顺序
# 对应 Table S1: 全部差异基因
table_s1_data <- read.csv("Table_S1_All_DEGs.csv") 

# 对应 Table S2: 4个 Hub Genes
table_s2_data <- read.csv("Table_S4_Hub_Genes.csv") 

# 对应 Table S3: Cox 回归系数
table_s3_data <- read.csv("Table_S5_Cox_Coefficients.csv") 

# 对应 Table S4: 仅 TCGA 的患者风险得分分层
table_s4_data <- read.csv("Table_S6_Patient_Risk_Score_Result.csv") 

print("正在构建整合版 Excel 工作簿...")

# 4. 创建一个新的空白 Excel 工作簿 (Workbook)
wb <- createWorkbook()

# 5. 添加名为 Table S1 到 S4 的工作表 (Sheets)
addWorksheet(wb, "Table S1")
addWorksheet(wb, "Table S2")
addWorksheet(wb, "Table S3")
addWorksheet(wb, "Table S4")

# 6. 将数据分别写入对应的工作表
writeData(wb, sheet = "Table S1", x = table_s1_data)
writeData(wb, sheet = "Table S2", x = table_s2_data)
writeData(wb, sheet = "Table S3", x = table_s3_data)
writeData(wb, sheet = "Table S4", x = table_s4_data)

# 可选优化：固定首行并自动加粗（让表格看起来更具顶刊质感）
headerStyle <- createStyle(textDecoration = "bold")
addStyle(wb, sheet = "Table S1", headerStyle, rows = 1, cols = 1:ncol(table_s1_data), gridExpand = TRUE)
addStyle(wb, sheet = "Table S2", headerStyle, rows = 1, cols = 1:ncol(table_s2_data), gridExpand = TRUE)
addStyle(wb, sheet = "Table S3", headerStyle, rows = 1, cols = 1:ncol(table_s3_data), gridExpand = TRUE)
addStyle(wb, sheet = "Table S4", headerStyle, rows = 1, cols = 1:ncol(table_s4_data), gridExpand = TRUE)
freezePane(wb, sheet = "Table S1", firstRow = TRUE)
freezePane(wb, sheet = "Table S2", firstRow = TRUE)
freezePane(wb, sheet = "Table S3", firstRow = TRUE)
freezePane(wb, sheet = "Table S4", firstRow = TRUE)

# 7. 保存并输出文件
output_filename <- "Supplementary_Tables.xlsx"
saveWorkbook(wb, output_filename, overwrite = TRUE)

print(paste("恭喜！多 Sheet 的整合表格", output_filename, "已经成功生成！"))