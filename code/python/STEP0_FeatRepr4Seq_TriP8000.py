#!/usr/bin/env python
# coding: utf-8

# ## 基于氨基酸残基序列的**```特征表示```**
# 
# #### 该模块用于计算给定蛋白质序列的
# - 氨基酸组成，AminoAcidsComposition：AAC=20；
# - 二肽组成，DipeptideComposition：DIP=400
# - 3-聚体（三肽）组成，tri-peptides： =8000   注：高维稀疏，效果不好！
# 共可得 8420 个特征(描述符)
# 
# 参考文献:
# 
# [1]: Reczko, M. and Bohr, H. (1994) The DEF data base of sequence based protein
# fold class predictions. Nucleic Acids Res, 22, 3616-3619.
# 
# [2]: Hua, S. and Sun, Z. (2001) Support vector machine approach for protein
# subcellular localization prediction. Bioinformatics, 17, 721-728.
# 
# [3]:Grassmann, J., Reczko, M., Suhai, S. and Edler, L. (1999) Protein fold class
# prediction: new methods of statistical classification.

# In[1]:


import numpy as np
import time 
import os

time_start = time.time()


### 练习测试：计算三肽构成成份

def CalculateTriPeptideComposition(ProteinSequence):
    """
    ########################################################################
    功能：计算给定蛋白质序列的三肽组成（Tripeptide composition）。
    输入：蛋白质序列（氨基酸残基序列）。
    输出：包含 8000 个三肽组成的字典。
    ########################################################################
    """
    LengthSequence = len(ProteinSequence)
    AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    Result = {}
    for i in AALetter:
        for j in AALetter:
            for k in AALetter:
                Tripeptide = i + j + k
                Result[Tripeptide] = round(ProteinSequence.count(Tripeptide) / (LengthSequence - 2) * 100, 3)
    return Result


def fasta2dict(inf):
    # 按行读取序列
    # 输入fasta文件，返回名称，序列
    global name
    dict = {}
    for line in inf:
        line = line.strip()
        if line.startswith('>'):
            name = line[1:]
            dict[name] = ''
        else:
            dict[name] += line
    return dict



# 测试算法是否OK？
# protein = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWTADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASGTGTWTTWTTWT"

# dict_tmp = CalculateTriPeptideComposition(protein)

# sum(dict_tmp.values())


# PDB1075数据集的特征工程：X-1075*8000; y-1075
# 共有1075条蛋白序列，前面525条为DNA结合蛋白，后面550条为非NA结合蛋白
label = np.concatenate((np.ones(525),np.zeros(550)))

seq_dir = '../Dataset/PDB1075_seq_FASTA'

# 获取子目录下的文件数目(序列数目)
file_num = len(os.listdir(seq_dir))
print(file_num)
trip8000 = np.zeros((file_num,20*20*20))

for k in range(file_num):
    seq_fn = 'seq_FASTA{}.fasta'.format(str(k+1))  # 序列文件名 
    fn = os.path.join(seq_dir,seq_fn)
    inf = open(fn)
    dict_kv = fasta2dict(inf)
    tmp_seq = list(dict_kv.values())[0]    # 读取第k条蛋白序列
    TriP_8000 = CalculateTriPeptideComposition(tmp_seq)  # 计算序列的三肽组成
    trip8000[k,:] = list(TriP_8000.values())     
    print('当前是第{}条序列，已经完成三肽组成计算'.format(k))   
      
time_end = time.time()
print('CPU计算时间：{}'.format(time_end - time_start))

# 计算特征保存成文件
savfn='../feat_npz/PDB1075_trip8000.npz'
np.savez(savfn, X = trip8000, y = label)
