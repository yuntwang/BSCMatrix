#!/usr/bin/env python
import sys
import gzip

bc1 = sys.argv[1]
bc2 = sys.argv[2]
bc3 = sys.argv[3]
barcode_10x = sys.argv[4]

########################读取序列文件
with open(bc1,'r') as fa:
    bc1_list = []
    for line in fa:
        # 去除末尾换行符
        line = line.replace('\n','')
        if line.startswith(">"):
             line = line.replace(">","")
             bc1_list.append(line)

with open(bc2,'r') as fa:
    bc2_list = []
    for line in fa:
        # 去除末尾换行符
        line = line.replace('\n','')
        if line.startswith(">"):
             line = line.replace(">","")
             bc2_list.append(line)

with open(bc3,'r') as fa:
    bc3_list = []
    for line in fa:
        # 去除末尾换行符
        line = line.replace('\n','')
        if line.startswith(">"):
             line = line.replace(">","")
             bc3_list.append(line)

########################读取10x barcode文件
with gzip.open(barcode_10x,'rt') as fa:
    barcode = []
    for line in fa:
        # 去除末尾换行符
        line = line.replace('\n','')
        barcode.append(line)

########################bc1 bc2 bc3 组合
bc_combine = [i + "-" + j + "-" + k for i in bc1_list for j in bc2_list for k in bc3_list]

with open("bmk_10x_barcode.xls","w") as file:
    for i,j in enumerate(bc_combine):
         file.write(j + "\t" + barcode[i] + "\n")







