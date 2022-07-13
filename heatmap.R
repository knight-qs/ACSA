# 将联配文件转化为phylip格式
library(devtools)
library(phylotools)
dat <- read.fasta("ali_500pr.fa")

dat2phylip(dat, outfile = "500pr.phy")

# Notepad++打开，alt+shift多光标批量删除，将phylip文件开头去除干净

# 绘制热图
num = 500 # 家族蛋白数量
a_c_count <- read.csv('query_results.csv', sep = ',', header = T) # 'X.'实际为'-'
sites <- read.csv('heatmap_sites.csv', sep = ',', header = F)
a_c_count1 <- a_c_count
row.names(a_c_count1) <- sites$V1
a_c_count2 <- t(a_c_count1)
a_c_count3 <- a_c_count2/num
a_c_count4 <- a_c_count3[-dim(a_c_count3)[1],] # 去除'X.'行

library(pheatmap)

a_c_count_heatmap <- pheatmap(a_c_count4,
                        cluster_cols = F,
                        cluster_rows = T,
                        cellheight = 10,
                        cellwidth = 10,
                        treeheight_row = 40,
                        legend = T,
                        fontsize_row = 5,
                        fontsize_col = 5,
                        angle_col = 45,
                        show_rownames = T,
                        show_colnames = T,
                        color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(10),
                        filename = 'a_c_count_heatmap.pdf' #调整好尺寸后再加这行保存文件
)
