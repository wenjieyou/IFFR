# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 23:41:05 2022

《数学建模与创新》案例讲解用

特征表示： gapdpc (参考) 和 gapdis 和 gapcov (KBS)
参数取值：ngap=0..L-1; kmer=1..L
@author: wenjie
"""


import scipy.io as sio
import numpy as np

import time

start = time.time()

# 读入 matlab 数据 cell 格式
data = sio.loadmat('..\pssm_mat\PDB1075_uniref50_pssmMatrix.mat');
allpssm=data['pssm']    # 读取字典对象，字典是另一种可变容器模型
ylab = data['label']
num_seq = allpssm.shape[0]

lenseq = []
for i in range(num_seq):
    pssmmatrix=allpssm[i][0]  
    lenseq.append(pssmmatrix.shape[0])
minlen = np.min(lenseq)
print('\n 全部{}序列中，最短序列长度为：{}'.format(num_seq,minlen))     
       

def gapdpc(pssm, ngap):      # 原始文献定义：基于相关系数
    L, kk = np.shape(pssm)
    # pssm = 1.0/(1+np.exp(-pssm))   # sigmoid变换：1/(1+e^-x)
    if ngap >= L:
        print("\nERROR!!!, Gap value cannot exceed L-1!\n")
        return
    # pssm矩阵的前块和后块，之间距离为gap
    col_FrontBlock = pssm[:L-ngap]
    col_BackBlock = pssm[ngap:]
    gap_dpc = np.zeros((kk,kk),dtype=float)
    for k in range(kk):
        gap_dpc[k] = np.mean((col_BackBlock.T*col_FrontBlock[:,k]).T,0)
    return gap_dpc   # 返回方阵


def gapcov(pssm, ngap):      # 我们KBS论文中的定义：基于IFFR算法
    L, n = pssm.shape    
    pssm = 1.0/(1+np.exp(-pssm))   # sigmoid变换：1/(1+e^-x)
    if L - ngap < 1:  # 序列长度m - 跳空距离 要大于1
        print("\nERROR!!!, Gap={} <= L-1={}!".format(ngap,L-1))
        return
    if ngap == 0:
        gapPSSM = pssm
    else:         
        gapPSSM = np.zeros((L-ngap,n),dtype=float)
        for j in range(L-ngap):
            gapPSSM[j] = pssm[[j,j+ngap]].mean(axis=0)    # 抽取指定两行
    gap_cov = np.dot(gapPSSM.T,gapPSSM)
    return gap_cov   # 返回方阵


n_gap = 50

# gapcov_4d = np.zeros((num_seq,50,20,20),dtype=np.float)
gapdpc_4d = np.zeros((num_seq,n_gap,20,20),dtype=np.float)


for i in range(num_seq):
    print('第{}条序列'.format(i))
    pssmmatrix=allpssm[i][0]  

    for ngap in range(n_gap):
        # 调用 gapdpc, gapcov
        gapXpssm = gapdpc(pssmmatrix,ngap)        
        gapdpc_4d[i,ngap,:] = gapXpssm
    

end = time.time()

print('CPU运行时间：%s 秒'%(end-start))

datX = gapdpc_4d   # mercov_4d  gapcov_4d


# 保存 npz格式文件
savfn=r'..\feat_npz\PDB1075_uniref50_gapdpc.npz'
np.savez(savfn,datX=datX,ylab=ylab)
    
