# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 22:20:46 2021

直接利用管道函数 (pipeline) 计算分类性能(30次的10折叉)
注：特征表示时用sigmoid，后续分类无标准化

@author: wenjie
"""

import numpy as np
from sklearn.svm import SVC
from sklearn.naive_bayes import MultinomialNB, GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
import time

start = time.time()


alldata=np.load('../feat_npz/PDB1075_uniref50_pssmdpc.npz')
# print(list(alldata))
datX = alldata['datX']
ylab = alldata['ylab']
print(datX.shape,ylab.shape)

num_seq = datX.shape[0]
ylab = ylab.ravel()

pipe = Pipeline(steps=[
    ('scaler', StandardScaler()),
    
    # ('selector',SelectKBest(f_classif, k=400)),  # 单变量特征选择
    # ('selector',SelectFromModel(SVC(kernel='linear',C=1.0),max_features=400)),
    # ('selector',RFECV(svc,step=0.1,cv=StratifiedKFold(2),scoring='accuracy')),
    # ('selector',RFE(svc,n_features_to_select=400,step=0.1)), 
      
    # ('predictor',KNeighborsClassifier(n_neighbors=3))
    # ('predictor',GaussianNB())   # MultinomialNB(alpha=1.2)
    # ('predictor',DecisionTreeClassifier(criterion='gini',max_depth=5))
    ('predictor',SVC(kernel='linear',probability=True))
    # ('predictor',RandomForestClassifier(criterion='entropy',n_estimators=1000))
    # ('predictor',XGBClassifier(use_label_encoder=False,eval_metric='logloss'))
    ])

alrs_mcc = []           # 10次10foldCV的结果: all_rand_state_MCC
alrs_acc = []
alrs_auc = []
alrs_sen = []
alrs_spe = []

for rs in range(10):    
    
    # 使用KFoldCV划分数据集. shuffle：每次划分是否进行洗牌
    skf = StratifiedKFold(n_splits=10,random_state=rs, shuffle=True)
    
    list_mcc = []
    list_acc = []
    list_auc = []
    list_sen = []
    list_spe = []
    
    for train_index,test_index in skf.split(datX,ylab):
        train_X, test_X = datX[train_index,:], datX[test_index,:]
        train_y, test_y = ylab[train_index], ylab[test_index]      
        
        # 浅层学习模型：训练分类器
        pipe.fit(train_X, train_y)
        ypred = pipe.predict(test_X)
        yprob = pipe.predict_proba(test_X)[:,1]
        # print('真实类={},\n预测类={}'.format(test_y,ypred))  
        
        acc = metrics.accuracy_score(test_y,ypred) # 计算ACC
        mcc = metrics.matthews_corrcoef(test_y,ypred) # 计算MCC相关系数
        auc = metrics.roc_auc_score(test_y,yprob)  # 计算AUC
        conf_mat = metrics.confusion_matrix(test_y, ypred)
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
    
    print('RandState={},MCC的总平均值=:{}\n'.format(rs, np.mean(list_mcc)))
    
    alrs_mcc.append(list_mcc)
    alrs_acc.append(list_acc)
    alrs_auc.append(list_auc)
    alrs_sen.append(list_sen)
    alrs_spe.append(list_spe)
    
## 保存结果
# output_fn = 'Result_PhosP_pssm20_POSSUM_ksepbigram_sig1_time30_10foldcv_ngap{}_STD_SVC.npz'.format(ngap)
# np.savez(output_fn,mcc=np.array(alrs_mcc),acc=np.array(alrs_acc),
#          auc=np.array(alrs_auc),sen=np.array(alrs_sen),spe=np.array(alrs_spe))
    
end = time.time()
print('\n所有MCC的平均值：{}'.format(np.array(alrs_mcc).mean(1)))
print('\n最终MCC平均值值：{}'.format(np.mean(np.array(alrs_mcc).mean(1))))
print('CPU运行时间：%s 秒'%(end-start))
