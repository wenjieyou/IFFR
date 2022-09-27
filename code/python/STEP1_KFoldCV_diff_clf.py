# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 22:54:25 2022

对PDB1075序列数据集，基于序列信息的特征表示方法，生成氨基酸组成(AAC20),二肽组成(DIP400)
和三肽组成(TriP8000),分别使用不同的分类器进行分类预测。

@author: wenjie
"""


import numpy as np
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
import time

start = time.time()

seed = 2

alldata=np.load('../feat_npz/PDB1075_uniref50_pssmaac.npz')
# print(list(alldata))
datX = alldata['datX']
ylab = alldata['ylab']
print(datX.shape,ylab.shape)

num_seq = datX.shape[0]

    
# 实例化：分类器
# clf = SVC(kernel='linear',C=1.0,probability=True)
# clf = KNeighborsClassifier(n_neighbors=5)
# clf = DecisionTreeClassifier(criterion='gini',max_depth=5)
clf = RandomForestClassifier(criterion='gini',n_estimators=500)

# 使用KFoldCV划分数据集. shuffle：每次划分是否进行洗牌
skf = StratifiedKFold(n_splits=10,random_state=seed, shuffle=True)

list_mcc = []
list_acc = []
list_auc = []
list_sen = []
list_spe = []

for train_index,test_index in skf.split(datX,ylab):
    train_X, test_X = datX[train_index,:], datX[test_index,:]
    train_y, test_y = ylab[train_index], ylab[test_index]
  
    
    # 数据标准化：先对训练集标准化再将规则用于测试集
    # StandardScaler使得经过处理的数据符合标准正态分布，即均值为0，标准差为1
    # stdScale = StandardScaler().fit(train_X)  # 生成标准化规则
    # train_X = stdScale.transform(train_X)     # 将规则应用于训练集
    # test_X = stdScale.transform(test_X)       # 将规则应用于测试集
    
    # 训练分类器
    clf.fit(train_X, train_y)
    ypred = clf.predict(test_X)
    yprob = clf.predict_proba(test_X)[:,1]
    # print('真实类={},\n预测类={}'.format(test_y,ypred))
       
    acc = metrics.accuracy_score(test_y,ypred) # 计算ACC
    mcc = metrics.matthews_corrcoef(test_y,ypred) # 计算MCC相关系数
    auc = metrics.roc_auc_score(test_y,yprob)  # 计算AUC
    conf_mat = metrics.confusion_matrix(test_y, ypred)   #计算混淆矩阵
    TP = conf_mat[1,1]
    TN = conf_mat[0,0]
    FP = conf_mat[0,1]
    FN = conf_mat[1,0]
    sen = TP / (float(TP+FN))
    spe = TN / (float(TN+FP))
    print('mcc={}'.format(mcc))

    list_mcc.append(mcc)
    list_acc.append(acc)
    list_auc.append(auc)
    list_sen.append(sen)
    list_spe.append(spe)

print('\n\n随机状态种子={}'.format(seed))
print('ACC的总平均值={}'.format(np.mean(list_acc)*100))
print('MCC的总平均值={}'.format(np.mean(list_mcc)))
print('AUC的总平均值={}'.format(np.mean(list_auc)))
    
# ## 保存结果
# np.savez('Result.txt',list_mcc=list_mcc,list_auc=list_auc)
    
end = time.time()

print('CPU运行时间：%s 秒'%(end-start))