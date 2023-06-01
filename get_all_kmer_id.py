import pickle
filename = 'D:/study/protease/data/raw/human.fasta'
save_path = 'D:/study/protease/data/processed/HTPS_16_kmer_acc'

all_kmers = dict()
def get_kmers(seq):
    # Extract your code into a function and print header for current kmer
    # print("%s\n################################" %name)
    kmers = set()
    k = 16
    for i in range(len(seq) - k + 1):
       kmer = seq[i:i+k]
       kmers.add(kmer)
    return kmers

with open(filename, 'r') as f:
    seq = ""
    key = ""
    for line in f.readlines():
        # Loop over lines in file
        if line.startswith(">"):
            # if we get '>' it is time for a new sequence
            if key and seq:
                # if it wasn't the first we should print it before overwriting the variables
                temp = get_kmers(seq)
                key = key.split('|')[1]
                all_kmers[key] = temp
            # store name after '>' and reset sequence
            key = line[1:].strip()
            seq = ""
        else:
            # accumulate kmer until we hit another '>'
           seq += line.strip()
        pass
    # when we are done with all the lines, print the last sequence
    key = key.split('|')[1]
    all_kmers[key] = temp
    print()

# save_path = save_path + str(len(list(all_kmers)))
with open(save_path, "wb") as fp:   #Pickling
    pickle.dump(all_kmers, fp)