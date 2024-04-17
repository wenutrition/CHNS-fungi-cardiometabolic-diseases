import pandas as pd
import lightgbm as lgb
from sklearn.model_selection import  train_test_split
from sklearn import metrics 
shap.initjs()
import numpy as np
from sklearn.model_selection import cross_val_score  

data = pd.read_stata('D:/真菌分析/data_pre/metabolism_predict_cluster34.dta') 
data = pd.read_stata('D:/真菌分析/data_pre/metabolism_predict_cluster24.dta') 
data = pd.read_stata('D:/真菌分析/data_pre/metabolism_predict_cluster14.dta') 
data = pd.read_stata('D:/真菌分析/data_pre/metabolism_predict_cluster4others.dta') 

data.index=data['SampleID']
train_labels=data['outcome']
data.drop(['outcome'],axis=1,inplace=True)
data.drop(['SampleID'],axis=1,inplace=True)
predictions=[]
labels=[]
SampleID=[]

evals_result={}
params = {
        
    'boosting_type': 'gbdt',  
    'objective': 'binary',
    'learning_rate': 0.05,       # Learning rate, controls size of a gradient descent step
    'metric': 'auc',  # Area under ROC curve as the evaulation metric

              }

kf=StratifiedKFold(n_splits=10,shuffle=True,random_state=123)  
kf
x,y=data,pd.DataFrame(train_labels)  #将数据dataframe化，后面进行iloc 等选项
for i,(train_index,valid_index) in enumerate(kf.split(x,y)): #这里需要输入y
    print("第",i+1,"次")
    x_train,y_train=x.iloc[train_index],y.iloc[train_index]
    x_valid,y_valid=x.iloc[valid_index],y.iloc[valid_index]  #取出数据
    lgb_train = lgb.Dataset(x_train, y_train,silent=True)
    lgb_eval  = lgb.Dataset(x_valid, y_valid, reference=lgb_train,silent=True)
    gbm = lgb.train(params, lgb_train, num_boost_round=10000, valid_sets=[lgb_train, lgb_eval],verbose_eval=100,
                        early_stopping_rounds=200,evals_result=evals_result)
    vaild_preds = gbm.predict(x_valid, num_iteration=gbm.best_iteration)
    X=x_train
    predictions.append(vaild_preds)
    labels.append(y_valid['outcome'])
    SampleID.append(y_valid.index)

predictions_validation=np.array([item for sublist in predictions for item in sublist])
labels_validation=np.array([item for sublist in labels for item in sublist])
sampleid_validation=np.array([item for sublist in SampleID for item in sublist])
auc_test_cv = metrics.roc_auc_score(labels_validation,predictions_validation)
print('模型在测试集上的效果是{:.5f}。'.format(auc_test_cv)) 

predict=pd.DataFrame()
predict['predict']=predictions_validation
predict['label']=labels_validation
predict['SampleID']=sampleid_validation
predict.to_excel('D:/真菌分析/Final results/predict_update/metabolis_fungi_discovery_auc.xlsx')
