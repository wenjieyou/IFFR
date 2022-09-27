# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 23:41:05 2022

《数学建模与创新》案例讲解用

特征表示： pssm-dpc(400) 和 pssm-aac(20) 和  pssm-cov(210)

@author: wenjie
"""


import numpy as np
import time

start = time.time()

# 读入 pssmMatrix 数据 npz 格式
data = np.load('..\pssm_mat\PDB1075_uniref50_pssmMatrix.npz', allow_pickle=True);

allpssm=data['allpssm']    
ylab = data['ylab']
num_seq = allpssm.shape[0]

lenseq = []
for i in range(num_seq):
    pssmmatrix=allpssm[i]  
    lenseq.append(pssmmatrix.shape[0])
minlen = np.min(lenseq)
print('\n 全部{}序列中，最短序列长度为：{}'.format(num_seq,minlen))     
       

def pssmdpc(pssm):      # 原始文献定义：基于相关系数，400维
    L, kk = np.shape(pssm)
    pssm = 1.0/(1+np.exp(-pssm))   # sigmoid变换：1/(1+e^-x)    
    # pssm矩阵的前块和后块，之间距离为gap
    col_FrontBlock = pssm[:L-1]
    col_BackBlock = pssm[1:]
    pssm_dpc = np.zeros((kk,kk),dtype=float)
    for k in range(kk):
        pssm_dpc[k] = np.mean((col_BackBlock.T*col_FrontBlock[:,k]).T,0)
    pssmdpc400 = pssm_dpc.flatten()
    return pssmdpc400   # 返回方阵


def pssmaac(pssm):      # 20维
    L, kk = np.shape(pssm)
    pssm = 1.0/(1+np.exp(-pssm))   # sigmoid变换：1/(1+e^-x)
    pssmaac20 = np.mean(pssm, axis=0)  
    return pssmaac20   # 返回20维向量


def pssmcov(pssm):      # 参考KBS论文中的定义：基于IFFR算法
    L, n = pssm.shape    
    pssm = 1.0/(1+np.exp(-pssm))   # sigmoid变换：1/(1+e^-x)
    pssm_cov = np.dot(pssm.T,pssm)
    
    tmp_triu = np.triu(pssm_cov)
    tmp_vec400 = tmp_triu.flatten()
    one_triu = np.triu(np.ones((20,20)))
    one_vec400 = one_triu.flatten().astype(bool)    
    pssmcov210 = tmp_vec400[one_vec400]      
    return pssmcov210   # 返回210维向量

# allpssmdpc = np.zeros((num_seq,20*20),dtype=np.float)
allpssmaac = np.zeros((num_seq,20),dtype=np.float)

for i in range(num_seq):
    print('第{}条序列'.format(i))
    pssmmatrix=allpssm[i]
    
    # 调用 pssmdpc, pssmcov
    Xpssm = pssmaac(pssmmatrix)        
    allpssmaac[i,:] = Xpssm
    
end = time.time()

print('CPU运行时间：%s 秒'%(end-start))

datX = allpssmaac  

# 保存 npz格式文件
savfn=r'..\feat_npz\PDB1075_uniref50_pssmaac.npz'
np.savez(savfn,datX=datX,ylab=ylab)
    
