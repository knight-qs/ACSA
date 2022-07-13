#############################################################################################################
# 蛋白氨基酸残基进化分析脚本: 统计所有位置的氨基酸种类、数量，筛选家族蛋白中保守的残基进化位点，导出绘制热图所需数据
#############################################################################################################

from typing import Counter
import csv
import os

path_amino_acid_existed = 'D:\\R\\python2\\amino_acid_existed.csv'
path_amino_acid_filted = 'D:\\R\\python2\\amino_acid_filted.csv'
path_query_results = 'D:\\R\\python2\\query_results.csv'
path_heatmap_sites = 'D:\\R\\python2\\heatmap_sites.csv'
path_query_protein = 'D:\\R\\python2\\query_protein.txt'
path_protein = 'D:\\R\\python2\\protein.txt'

if os.path.exists(path_amino_acid_existed):
    os.remove("amino_acid_existed.csv")
if os.path.exists(path_amino_acid_filted):
    os.remove("amino_acid_filted.csv")
if os.path.exists(path_query_results):
    os.remove("query_results.csv")
if os.path.exists(path_heatmap_sites):
    os.remove("heatmap_sites.csv")

def csv_site_writer(site):
    with open(path_amino_acid_existed, 'a') as existed_csvfile:
        site = str(site) + '\n'     #这里加了str()，因为这里的query_sites由数字组成，而之前的脚本读入后query_sites由字符串组成
        existed_csvfile.writelines(site)

# 读列表，列表内均为元组
def csv_existed_writer(existed_a_c):
    with open(path_amino_acid_existed, 'a') as existed_csvfile:
        existed = csv.writer(existed_csvfile)
        existed.writerow(existed_a_c)

def csv_filte_writer(flit):
    with open(path_amino_acid_filted, 'a') as filted_csvfile:
        flit = flit + '\n'
        filted_csvfile.writelines(flit)

# 读列表，此列表内元素均为字典
def csv_writer(query_result):
    dict = query_result
    with open(path_query_results, 'a', newline = '', encoding = 'utf-8') as n:
            writer = csv.DictWriter(n, fieldnames = amino_acid_all)  
            writer.writeheader()  # 写入列名
            for p in dict:
                writer.writerow(p)  # 写入数据

def csv_heatmap_sites_writer(heatmap_site):
    with open(path_heatmap_sites, 'a') as heatmap_sites:
        heatmap_site = heatmap_site + '\n'
        heatmap_sites.write(heatmap_site)

amino_acid = []
with open(path_query_protein) as a_c:
    char = a_c.read(1)
    while char:
        if char != '\n':
            amino_acid.append(char)
            char = a_c.read(1)
        else:
            char = a_c.read(1)

query_sites = []
for z in range(len(amino_acid)):
    query_sites.append(z)
del query_sites[0]      # 蛋白质位置索引从1开始，而非0
query_sites.append(len(amino_acid))     # range()起始为0，最后一位为len(amino_acid)-1，所以末尾需要补全

query_info = {}.fromkeys(query_sites)
d = 0
for e in query_sites:
    query_info[e] = amino_acid[d]
    d += 1

amino_acid_all = ['Z', 'X', 'H', 'R', 'W', 'K', 'E', 'D', 'Y', 'F', 'L', 'N', 'P', 'C', 'I', 'Q', 'T', 'V', 'A', 'M', 'S', 'G', '-']

#######################################################
#   开始统计、筛选
#######################################################
every_protein = []
for f in query_sites:
    site = f      #这里去掉int()，因为这里的query_sites由数字组成，而之前的脚本读入后query_sites由字符串组成
    aim_amino_acid = query_info[f]
    protein = []
    sum = []
    g = 0
    h = 0
    with open(path_protein, encoding = "utf-8") as p:
        for i in p.readlines():
            protein.append(i)
        
        # 判断相应氨基酸在phylip文件中的实际索引位置
        while g < site:
            h += 1
            if protein[0][h] != '-':
                g += 1

        if protein[0][h] == aim_amino_acid:
            # 出现氨基酸频数
            print('Amino acid site',site,'is OK')
            for j in protein:
                sum.append(j[h])
            result = Counter(sum)
            print('site', site, ':', result, '\n')  # 需要手动复制粘贴到文本文件，不如items()输出方便

            # 字典排序
            tuple = zip(result.values(),result.keys())      # 返回列表，列表内各元素为元组，key和value位置对调
            a_acid_sort = list(sorted(tuple, reverse = True))       # 按各元组第一位（这里是value）排序
            csv_site_writer(f)
            csv_existed_writer(a_acid_sort)

            # 筛选突变位点
            # 标准：被blast的蛋白质的对应位置氨基酸种类已知，若联配后该位置该种类氨基酸数量并非第一
            # 则认为该位置为突变位点，并返回该位点，及该位点上数量更多的氨基酸种类、相应数量
            filte = []
            for s in a_acid_sort:
                if aim_amino_acid != s[1]:
                    filte.append(s)
                else:
                    break
            if filte != []:
                trans = ''
                for t in filte:
                    t0 = t[0]
                    t1 = t[1]
                    trans = trans + str(t1) + ':' + str(t0) + ' '
                trans = 'site' + str(f) + '    ' + 'here is ' + aim_amino_acid + ':' + str(s[0]) + '    ' + trans
                csv_filte_writer(trans)     # 上一行的f加了str()，因为这里的query_sites由数字组成，而之前的脚本读入后query_sites由字符串组成

                # 全部种类氨基酸频数
                amino_acids_whole = {}.fromkeys(amino_acid_all)
                for k in amino_acid_all:
                    amino_acids_whole[k] = 0
                amino_acids_whole.update(result)
                every_protein.append(amino_acids_whole)
                csv_heatmap_sites_writer(str(f))    # 收集热图对应的氨基酸位置

        else:
            print('ERROR')
            print('Now amino acid site is', g ,', and phylip site is', h)
            print('This amino acid is:', protein[0][h], ', not your aim:', aim_amino_acid)

# 汇总，准备导入R，绘制热图
csv_writer(every_protein)
