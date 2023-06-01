import requests
import pickle
import os
# root = '/study/protease/data/align/FIX'
# save_path = '/study/protease/data/align/' + 'FIX_pep_accessnum'
# with open(root, "rb") as fp:   # Unpickling
#    positive_data = pickle.load(fp)
# peplist = []
# accnum = []
# n = 0
# for pep in positive_data:
#    if len(positive_data[pep]) == 1:
#       peplist.append(pep)
#       accnum.append(positive_data[pep][0].split('|')[1])
#    else:
#       n += 1
# accnum = list(set(accnum))

# filename = 'D:/study/protease/data/raw/HTPS_db.fasta'
# acc_list = []
# with open(filename, 'r') as f:
#    seq = ""
#    key = ""
#    for line in f.readlines():
#       # Loop over lines in file
#       if line.startswith(">"):
#          # if we get '>' it is time for a new sequence
#          if key and seq:
#             # if it wasn't the first we should print it before overwriting the variables
#             acc_list.append(key.split('|')[1])
#          # store name after '>' and reset sequence
#          key = line[1:].strip()
#          seq = ""
#       else:
#          # accumulate kmer until we hit another '>'
#          seq += line.strip()
#       pass
#    print()
#    # when we are done with all the lines, print the last sequence
#    acc_list.append(key.split('|')[1])
# f.close()
path = 'D:/study/protease/data/processed/HTPS_16_kmer_acc'
with open(path, "rb") as fp:   # Unpickling
   all_kmer = pickle.load(fp)
acc_list = list(all_kmer.keys())
for id in acc_list:
   save_path = '/study/protease/data/pdb/' + id + '.pdb'
   if os.path.isfile(save_path):
      print("File already exists!")
   else:
      url = 'https://alphafold.ebi.ac.uk/files/AF-' + id + '-F1-model_v1.pdb'
      response = requests.get(url)
      with open(save_path, 'wb') as f:
         f.write(response.content)