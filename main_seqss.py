import pandas as pd
import numpy as np
import os
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
import sys

#
protease = int(sys.argv[1])
another = int(sys.argv[2])
# protease = 0
# another = 1
folderlist = [r'D:\study\protease\data\inputs\multiple\GluC.csv',r'D:\study\protease\data\inputs\multiple\MMP9.csv']
df = pd.read_csv(folderlist[protease])
df_train, df_test = train_test_split(df,test_size=0.2,random_state=42)

negative_data = pd.read_csv(r'D:\study\protease\data\inputs\multiple\neg_nor_seqss.csv')
negdf = negative_data.sample(len(df_train))
neg_diff = negative_data.merge(negdf, how='outer', indicator=True)
neg_diff = neg_diff.loc[neg_diff['_merge'] != 'both']
neg_diff = neg_diff.iloc[:,:5]


train = pd.concat([df_train,negdf])
### merge GO PPI info
roots = [r'D:\study\protease\data\inhouse_data\GlucGO',r'D:\study\protease\data\inhouse_data\MMP9GO']
GOs = ['Gluc_all_protein_GO.csv','MMP9_all_protein_GO.csv']
PPIs = ['Gluc_PPI.csv','MMP9_PPI.csv']

godata = pd.read_csv(os.path.join(roots[protease],GOs[protease]))
ppidata = pd.read_csv(os.path.join(roots[protease],PPIs[protease]))

merged_df = pd.merge(train, godata, on='accid', how='inner')
merged_df = pd.merge(merged_df, ppidata, on='accid', how='inner')

merged_df = merged_df.sample(frac=1,random_state=55)
### get X and y
peptides = merged_df['peptide'].values
y  = merged_df['label'].values

alphabet_file = "https://raw.githubusercontent.com/brunoalvarez89/data/master/algorithms_in_bioinformatics/part_3/alphabet"
alphabet = np.loadtxt(alphabet_file, dtype=str)
blosum_file = "https://raw.githubusercontent.com/brunoalvarez89/data/master/algorithms_in_bioinformatics/part_3/blosum50"
_blosum50 = np.loadtxt(blosum_file, dtype=float).reshape((24, -1)).T
blosum50 = {}
for i, letter_1 in enumerate(alphabet):
    blosum50[letter_1] = {}
    for j, letter_2 in enumerate(alphabet):
        blosum50[letter_1][letter_2] = _blosum50[i, j] / 5.0

def encode(peptides, encoding_scheme, alphabet):
    encoded_peptides = []
    for peptide in peptides:
        encoded_peptide = []
        for peptide_letter in peptide:
            for alphabet_letter in alphabet:
                encoded_peptide.append(encoding_scheme[peptide_letter][alphabet_letter])
        # add a 1 (bias)
        # encoded_peptide.append(1)
        # store peptide
        encoded_peptides.append(encoded_peptide)
    return np.array(encoded_peptides)

seqs = encode(peptides, blosum50, alphabet)
#######structure encoding
ss = "SHL"
pepdss = merged_df['dss'].values
dss = np.zeros((len(pepdss), len(pepdss[0])*len(ss)), dtype=int)
for i, pep in enumerate(pepdss):
    for j, aa in enumerate(pep):
        dss[i, j*len(ss) + ss.index(aa)] = 1


tempdf = merged_df[['sasa', 'CC_GO1','CC_GO2', 'CC_GO3', 'BP_GO1', 'BP_GO2', 'BP_GO3', 'MF_GO1', 'MF_GO2','MF_GO3', 'PPI_1', 'PPI_2']]
# tempdf = merged_df[['sasa']]
X = np.hstack((seqs,dss,tempdf.values))
# X = seqs

# Create a Random Forest classifier on training data
fold = 5
## optimal n_estimators
param_test1 = {"n_estimators":range(1,101,10)}
gsearch1 = GridSearchCV(estimator=RandomForestClassifier(),param_grid=param_test1,
                        scoring='roc_auc',cv=fold)
gsearch1.fit(X,y)
print(gsearch1.best_params_)
print("best accuracy:%f" % gsearch1.best_score_)
optimal_estimators = gsearch1.best_params_['n_estimators']

## optimal max_features
param_test2 = {"max_features":range(1,21,1)}
gsearch1 = GridSearchCV(estimator=RandomForestClassifier(n_estimators=optimal_estimators,
                        random_state=10),
                        param_grid = param_test2,scoring='roc_auc',cv=fold)
gsearch1.fit(X,y)
print(gsearch1.best_params_)
print('best accuracy:%f' % gsearch1.best_score_)
optimal_features = gsearch1.best_params_['max_features']

# get optimal rf clf
clf_rf = RandomForestClassifier(n_estimators=optimal_estimators,max_features=optimal_features,
                             oob_score=True,random_state=10)

clf_rf.fit(X,y)
rf_acc = round(clf_rf.oob_score_,4)
print("rf_train_accuracy: ", rf_acc)

# save
saveroot=r'D:\study\protease\data\inputs\multiple\results'
pickle.dump(clf_rf, open(os.path.join(saveroot,'model_'+str(protease+1)), "wb"))
pickle.dump(df_test, open(os.path.join(saveroot,'test_'+str(protease+1)), "wb"))
pickle.dump(neg_diff, open(os.path.join(saveroot,'neg_'+str(protease+1)), "wb"))