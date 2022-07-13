# ACSA

1. 在NCBI比对蛋白质序列，获取蛋白质家族，如bla_500pr.fasta文件

2. 使用在线工具，联配各蛋白质序列，如ali_500pr.fa文件

3. 使用R脚本中的第一部分，将联配文件ali_500pr.fa文件转化为phylip格式，如500pr.phy

4. 使用文本编辑器多光标编辑功能，如notepad++，将各行开头删除干净，如protein.txt

5. 在bla_500pr.fasta文件中复制粘贴自己的蛋白质序列，生成新文件，如query_protein.txt

6. 运行ACSA.py，筛选家族蛋白中保守的残基进化位点

7. 结果文件amino_acid_existed.csv，所有位点(以目标蛋白为基准)氨基酸种类、频数统计

8. 结果文件amino_acid_filted.csv, 筛选目标蛋白对应的氨基酸频数并非最大的位点

9. 结果文件query_results.csv，筛选出的氨基酸种类及对应频数，下游R脚本绘制热图所需数据文件

10. 结果文件heatmap_sites.csv，筛选出的蛋白质位点，下游R脚本绘制热图所需数据文件
