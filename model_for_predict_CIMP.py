import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.metrics import *
from sklearn.linear_model import LogisticRegression
import os, os.path
import argparse

def process_data(num, folder='8_feature'):
    df  = pd.read_csv('data/group_info_for_34samples.csv', delimiter=',')
    labels = {}
    for i, row in df.iterrows():
        labels[row['sample_name']] = 1 if row['group_id'] == 1 else 0
    #print(labels)

    X = pd.read_csv('%s/%dsamples_on_features_proteincoding_FPKM.txt'%(folder, num), delimiter='\t')
    print(X.shape)
    cells = X.columns[1:]
    features = X['feature'].values
    
    X = X.T
    X.columns = features

    X.to_csv('%s/%d.csv'%(folder, num))

def read_data(num, folder):
    df  = pd.read_csv('data/group_info_for_%dsamples.csv'%(num), delimiter=',')
    labels = {}
    for i, row in df.iterrows():
        labels[row['sample_name']] = 1 if row['group_id'] == 1 else 0
    #print(labels)

    X = pd.read_csv('%s/%d.csv'%(folder, num), delimiter=',')
    print(X.shape)
    X['label'] = X['cell'].apply(lambda x: labels[x])
    
    return X

def preprocess(X, means=None, stds=None):
    means1, stds1 = [], []
    for i, col in enumerate(X.columns):
        if means is None:
            mean, std = X[col].mean(), X[col].std()
            means1.append(mean)
            stds1.append(std)
        else:
            mean, std = means[i], stds[i]
        if std != 0:
            X[col] = (X[col]-mean)/std
    return X, means1, stds1
        

def train(args):
    df_train = read_data(34, args.dir)
    df_test = read_data(255, args.dir)

    X_train, y_train = df_train[[col for col in df_train.columns if col not in ['cell', 'label']]], df_train['label'].values
    X_test, y_test = df_test[[col for col in df_test.columns if col not in ['cell', 'label']]], df_test['label'].values

    #model = 'random_forest' # logistic, logistic_l1

    if args.model != 'rf':
        X_train, means, stds = preprocess(X_train)
        X_test, _, _ = preprocess(X_test, means, stds)
        X_train.to_csv(os.path.join(args.dir, '34_process.csv'))
    

    #print(y_test)
    if args.model == 'rf':
        clf = RandomForestClassifier(max_depth=2, random_state=0)
    elif args.model == 'lr':
        clf = LogisticRegression(random_state=0, penalty='l1', C=1, solver='saga')

    clf.fit(X_train, y_train)
    #y_pred = clf.predict(X_test)
    y_pred_prob = clf.predict_proba(X_test)[:, 1]

    fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred_prob)
    
    print('auc', metrics.auc(fpr, tpr))
    
    

    roc_auc = auc(fpr, tpr)
    plot(fpr, tpr, roc_auc, args)

    if args.model == 'rf':
        fig, ax = plt.subplots(figsize=(7, 5))
        #print(len(clf.feature_importances_), np.max(clf.feature_importances_), np.min(clf.feature_importances_))
        sorted_idx = np.argsort(-clf.feature_importances_)
        #print(sorted_idx, X_train.columns[sorted_idx][:20], clf.feature_importances_[sorted_idx][:20])
        with open('%s/feature_importance_rf.txt'%(args.dir), 'w') as f:
            f.write('column\tfeature importance\n')
            for col, fi in zip(X_train.columns[sorted_idx], clf.feature_importances_[sorted_idx]):
                f.write('%s\t%f\n'%(col, fi))
        plt.barh(X_train.columns[sorted_idx][:20], clf.feature_importances_[sorted_idx][:20])
        plt.xlabel("Random Forest Feature Importance")
        #ax.bar(range(len(clf.feature_importances_)), clf.feature_importances_)
        #ax.set_title("Feature Importances")
        #plt.savefig('importance.png')
        plt.savefig(os.path.join(args.dir, 'rf_importance.png'))
    else:
        #print('coef', clf.coef_)
        with open(os.path.join(args.dir, 'lr_coef.txt'), 'w') as f:
            f.write('feature\tcoef\n')
            for col, coef in zip(X_train.columns, clf.coef_[0]):
                f.write('%s\t%f\n'%(col, coef))

def plot(fpr, tpr, roc_auc, args):
    plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
            lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.savefig(os.path.join(args.dir, args.model+'_roc_auc.png'))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Classifying scATAC data')

    parser.add_argument('--dir', default='2m_feature', type=str, help='data directory')
    parser.add_argument('--threshold', default=0.95, type=float, help='the threshold value for discriminating pos/neg samples')
    parser.add_argument('--model', default='rf', type=str, help='the model: rf, lr')

    args = parser.parse_args()
    
    #process_data(34, folder='2m_feature')
    #process_data(255, folder='2m_feature')
    train(args)
    #train('8_feature')
    