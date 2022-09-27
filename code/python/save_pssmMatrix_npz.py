# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 21:31:59 2022

对 PSIBLAST 生成的单个 pssmMatrix.txt 提取其中想要的 L*20 PSSM矩阵 
（通常只要前面的20个氨基酸的打分值，前面的20列）

结果：allpssm - 列表， 每个值对应一个pssm矩阵
      ylab - 类标

（对应于 MATLAB中，save_pssmMatrix_mat.m）

@author: wenjie
"""

import numpy as np
import os

pssmdir = '../Dataset/PDB1075_pssmMatrix_uniref50/'
listfile = os.listdir(pssmdir)

allpssm = []

# 获取子目录下的文件数目(序列数目)
file_num = len(listfile)

for k in range(file_num):
    pssm_fn = 'PSSMmatrix{}.txt'.format(str(k+1))  # pssmMatrix1.txt 文件名 
    fn = os.path.join(pssmdir,pssm_fn)
    
    with open(fn, 'r') as pssmfile:  
        list_pssm = []
        count = 0
        for eachline in pssmfile:
            count += 1
            if count <= 3:   # 跳过前3行
                continue
            if not len(eachline.strip()):   # 出现 空行 表明文件尾，跳出
                break
            oneline = eachline[10:-1]       # 读取 每行数据（跳过前10行氨基酸行号）     
            oneline_pssm = [int(x) for k,x in enumerate(oneline.split()) 
                         if x !='' and k < 20]   # 取非空的数值，前20个
            list_pssm.append(oneline_pssm) 
        
        onepssm = np.array(list_pssm)
    allpssm.append(onepssm)
   
    
label = np.concatenate((np.ones(525),np.zeros(550)))

savfn=r'..\feat_npz\PDB1075_uniref50_pssmMatrix.npz'
np.savez(savfn,allpssm=allpssm,ylab=label)
    