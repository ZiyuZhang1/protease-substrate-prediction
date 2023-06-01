import subprocess
import pickle
import pandas as pd
# # Run file.py five times and record output to results.txt
# with open("seq_results.txt", "w") as f:
#     for i in range(5):
#         subprocess.run(["python", "D:/study/protease/code/My_PP/main_seq.py", '0','1'])
#         subprocess.run(["python", "D:/study/protease/code/My_PP/main_seq.py", '1', '0'])
#
#         neg1 = pickle.load(open(r'D:\study\protease\data\inputs\multiple\results\neg_1', "rb"))
#         neg2 = pickle.load(open(r'D:\study\protease\data\inputs\multiple\results\neg_2', "rb"))
#         neg = pd.merge(neg1,neg2,how='inner')
#         pickle.dump(neg, open(r'D:\study\protease\data\inputs\multiple\results\neg', "wb"))
#
#         subprocess.run(["python", "D:/study/protease/code/My_PP/mian_seq_test.py"])
#         process = subprocess.Popen(["python", "D:/study/protease/code/My_PP/main_seq_test.py"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         stdout, stderr = process.communicate()
#         f.write(stdout.decode())


# # # Run file.py five times and record output to results.txt
# with open("bad_multi_10nge_results.txt", "w") as f:
#     for i in range(5):
#         subprocess.run(["python", "D:/study/protease/code/My_PP/main_seqss.py", '0','1'])
#         subprocess.run(["python", "D:/study/protease/code/My_PP/main_seqss.py", '1', '0'])
#
#         neg1 = pickle.load(open(r'D:\study\protease\data\inputs\multiple\results\neg_1', "rb"))
#         neg2 = pickle.load(open(r'D:\study\protease\data\inputs\multiple\results\neg_2', "rb"))
#         neg = pd.merge(neg1,neg2,how='inner')
#         pickle.dump(neg, open(r'D:\study\protease\data\inputs\multiple\results\neg', "wb"))
#         # subprocess.run(["python", "D:/study/protease/code/My_PP/mian_seqss_test.py"])
#         subprocess.run(["python", "D:/study/protease/code/My_PP/multiple_test.py"])
#         # process = subprocess.Popen(["python", "D:/study/protease/code/My_PP/main_seq_test.py"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         # stdout, stderr = process.communicate()
#         # f.write(stdout.decode())


# Run file.py five times and record output to results.txt
with open("multiple_class10neg_results.txt", "w") as f:
    for i in range(5):
        subprocess.run(["python", "D:/study/protease/code/My_PP/df_split.py", '0','1'])
        subprocess.run(["python", "D:/study/protease/code/My_PP/df_split.py", '1', '0'])

        neg1 = pickle.load(open(r'D:\study\protease\data\inputs\multiple\results\neg_1', "rb"))
        neg2 = pickle.load(open(r'D:\study\protease\data\inputs\multiple\results\neg_2', "rb"))
        neg = pd.merge(neg1,neg2,how='inner')
        pickle.dump(neg, open(r'D:\study\protease\data\inputs\multiple\results\neg', "wb"))

        process = subprocess.Popen(["python", "D:/study/protease/code/My_PP/multiple_test.py"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        f.write(stdout.decode())

