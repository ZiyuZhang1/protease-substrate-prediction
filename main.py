# get sequence info
## get pratease substrate list of prptides, accessID
### positive data
#### in house data
import numpy as np
import pandas as pd
import pickle
import random
from pymol import cmd
from pymol import stored
import os
from sklearn import preprocessing


all_kmer_path = 'D:/study/protease/data/processed/HTPS_16_kmer_acc'
def get_accid(pep_list):
    n=1
    accid = list()
    for pep in pep_list:
        for id in all_kmer_dict:
            if pep in all_kmer_dict[id]:
                accid.append(id)
                break
    return accid

infile = 'D:/study/protease/data/inhouse_data/MMP9_cleavages_R1.csv'

in_data = pd.read_csv(infile)
pos_data = pd.DataFrame()
pos_data['peptide'] = in_data['Cleavage_environment']
pos_data['accid'] = in_data['Accession']

#filter
filter_data = pos_data[pos_data['peptide'].str.contains('-')]
pos_data = pos_data.append(filter_data)
pos_data = pos_data.drop_duplicates(subset=['peptide','accid'],keep=False)
# labeling
lable = [1] * len(pos_data)
pos_data['label'] = lable

### negative data
# neg_size = 5*len(pos_data)
neg_size = len(pos_data)
with open(all_kmer_path, "rb") as fp:   # Unpickling
   all_kmer_dict = pickle.load(fp)

all_pdb_kmers = set()
pdb_list = list()
for i in os.listdir('D:/study/protease/data/pdb'):
    pdb_list.append(i.split('.')[0])


for id in pdb_list:
    all_pdb_kmers.update(all_kmer_dict[id])


neg_set = set(all_pdb_kmers) - set(pos_data['peptide'].tolist())
neg_list = random.sample(neg_set,neg_size)
neg_acc = get_accid(neg_list)
neg_data = pd.DataFrame()
neg_data['peptide'] = neg_list
neg_data['accid'] = neg_acc
neg_data['label'] = [0]*len(neg_data)
data = pd.concat([pos_data,neg_data])


# get structure info
pdbroot = 'D:/study/protease/data/pdb/'
dsslist = list()
sasalist = list()
for i, pep in enumerate(data['peptide'].tolist()):
    peptide_chain = pep
    pdb_path = pdbroot + data['accid'].tolist()[i] + '.pdb'
    if os.path.isfile(pdb_path):
        # refresh
        cmd.delete("all")
        # load in PDB structure
        cmd.load(pdb_path)
        # Calculate the secondary structure using DSSP
        cmd.dss()

        # ss = cmd.get_fastastr("all")
        cmd.select("peptide", f"pepseq {peptide_chain}")
        cmd.create("pep", "peptide")
        stored.ss = ""
        cmd.iterate('(pep and n. CA)', 'stored.ss = stored.ss + ("%1s"%ss)')
        sasa = cmd.get_area('pep')
        if len(stored.ss) == 16:
            dsslist.append(stored.ss)
        elif len(stored.ss) == None:
            dsslist.append(np.nan)
        else:
            dsslist.append(stored.ss[:16])
        sasalist.append(sasa)
    else:
        dsslist.append(np.nan)
        sasalist.append(np.nan)
    # Print the secondary structure information
    # break

# get raw dataframe
data['dss'] = dsslist
data['sasa'] = sasalist


savepath = 'D:/study/protease/data/inputs/' + infile.split('/')[-1].split('_')[0] + '_raw_seqss.csv'
# savepath = 'D:/study/protease/data/inputs/' + infile.split('/')[-1].split('_')[0] + '_5neg_raw_seqss.csv'
data.to_csv(savepath,index=False)

# drop NA
data = data.dropna()
# scale sasa
data['sasa'] = preprocessing.scale(np.array(data['sasa']))

savepath2 = 'D:/study/protease/data/inputs/' + infile.split('/')[-1].split('_')[0] + '_nor_seqss.csv'
# savepath2 = 'D:/study/protease/data/inputs/' + infile.split('/')[-1].split('_')[0] + '_5neg_nor_seqss.csv'
data.to_csv(savepath2,index=False)