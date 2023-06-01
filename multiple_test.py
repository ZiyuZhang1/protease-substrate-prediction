import pandas as pd
import numpy as np
import pickle
import os
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report
import warnings
warnings.filterwarnings("ignore")


# neg1 = pickle.load(open(r'D:\study\protease\data\inputs\multiple\results\neg_1', "rb"))
# neg2 = pickle.load(open(r'D:\study\protease\data\inputs\multiple\results\neg_2', "rb"))
# neg = pd.merge(neg1,neg2,how='inner')
# pickle.dump(neg, open(r'D:\study\protease\data\inputs\multiple\results\neg', "wb"))
print('multiple_classification')
test1 = pickle.load(open(r'D:\study\protease\data\inputs\multiple\results\test_1', "rb"))
test2 = pickle.load(open(r'D:\study\protease\data\inputs\multiple\results\test_2', "rb"))
negdf = pickle.load(open(r'D:\study\protease\data\inputs\multiple\results\neg', "rb"))
neg = negdf.sample(n=len(test2)*10,random_state=42)

# set1 = set(test1['peptide'].tolist())&set(test2['peptide'].tolist())
### get test data
all_test = pd.concat([test1, test2, neg])
all_test1 = all_test.copy()
all_test2 = all_test.copy()
all_test1.iloc[len(test1):(len(test1)+len(test2)), 2] = 0
all_test2.iloc[0:len(test1), 2] = 0
### merge GO PPI info
roots = [r'D:\study\protease\data\inhouse_data\GlucGO',r'D:\study\protease\data\inhouse_data\MMP9GO']
GOs = ['Gluc_all_protein_GO.csv','MMP9_all_protein_GO.csv']
PPIs = ['Gluc_PPI.csv','MMP9_PPI.csv']
protease = 0
godata = pd.read_csv(os.path.join(roots[protease],GOs[protease]))
ppidata = pd.read_csv(os.path.join(roots[protease],PPIs[protease]))
merged_df1 = pd.merge(all_test1, godata, on='accid', how='inner')
merged_df1 = pd.merge(merged_df1, ppidata, on='accid', how='inner')

protease = 1
godata = pd.read_csv(os.path.join(roots[protease],GOs[protease]))
ppidata = pd.read_csv(os.path.join(roots[protease],PPIs[protease]))
merged_df2 = pd.merge(all_test2, godata, on='accid', how='inner')
merged_df2 = pd.merge(merged_df2, ppidata, on='accid', how='inner')

# if merged_df1['peptide'].tolist() == merged_df2['peptide'].tolist():
#     a = 'yes'

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
#######load model
M1 = pickle.load(open(r'D:\study\protease\data\inputs\multiple\results\model_1', "rb"))
M2 = pickle.load(open(r'D:\study\protease\data\inputs\multiple\results\model_2', "rb"))

peptides1 = merged_df1['peptide'].values
y1  = merged_df1['label'].values
seqs1 = encode(peptides1, blosum50, alphabet)
#######structure encoding
ss = "SHL"
pepdss = merged_df1['dss'].values
dss1 = np.zeros((len(pepdss), len(pepdss[0])*len(ss)), dtype=int)
for i, pep in enumerate(pepdss):
    for j, aa in enumerate(pep):
        dss1[i, j*len(ss) + ss.index(aa)] = 1

tempdf1 = merged_df1[['sasa', 'CC_GO1','CC_GO2', 'CC_GO3', 'BP_GO1', 'BP_GO2', 'BP_GO3', 'MF_GO1', 'MF_GO2','MF_GO3', 'PPI_1', 'PPI_2']]
# tempdf = merged_df[['sasa']]
X1 = np.hstack((seqs1,dss1,tempdf1.values))
y_proba1 = M1.predict_proba(X1)[:, 1]


peptides2 = merged_df2['peptide'].values
y2  = merged_df2['label'].values
seqs2 = encode(peptides2, blosum50, alphabet)
#######structure encoding
ss = "SHL"
pepdss = merged_df2['dss'].values
dss2 = np.zeros((len(pepdss), len(pepdss[0])*len(ss)), dtype=int)
for i, pep in enumerate(pepdss):
    for j, aa in enumerate(pep):
        dss2[i, j*len(ss) + ss.index(aa)] = 1

tempdf2 = merged_df2[['sasa', 'CC_GO1','CC_GO2', 'CC_GO3', 'BP_GO1', 'BP_GO2', 'BP_GO3', 'MF_GO1', 'MF_GO2','MF_GO3', 'PPI_1', 'PPI_2']]
# tempdf = merged_df[['sasa']]
X2 = np.hstack((seqs2,dss2,tempdf2.values))
y_proba2 = M2.predict_proba(X2)[:, 1]

results_df = pd.DataFrame()
results_df = merged_df1[['peptide','accid']]
def get_label(row):
    if row['lab1'] == 1:
        return 1
    elif row['lab2'] == 1:
        return 2
    else:
        return 0
def get_prab(row):
    if row['prab1'] > row['prab2'] and row['prab1'] >= 0.5:
        return 1
    elif row['prab2'] > row['prab1'] and row['prab2'] >= 0.5:
        return 2
    else:
        return 0
results_df['lab1'] = y1
results_df['lab2'] = y2
results_df['label'] = results_df.apply(get_label, axis=1)

results_df['prab1'] = y_proba1
results_df['prab2'] = y_proba2
results_df['predict'] = results_df.apply(get_prab, axis=1)

acc = accuracy_score(results_df['label'], results_df['predict'])
print(acc)

y_true1, y_pred1 = y1, M1.predict(X1)
test_rf1 = classification_report(y_true1, y_pred1)
print(test_rf1)

y_true2, y_pred2 = y2, M2.predict(X2)
test_rf2 = classification_report(y_true2, y_pred2)
print(test_rf2)

results_df.to_csv('final.csv',index=False)