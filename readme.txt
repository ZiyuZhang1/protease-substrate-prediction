If you want to generate data:
1. Use get_all_kmer_id.py to generate all 16-kmers of human.faste. Then later you can choose negative data from this generated dictionary.
2. Use get pdb.py to download all structure pdb files.
3. Use main.py to generate sequence and structure information with given 16-kmer (positive or negative)
4. download GO and PPI information. (same method and data source from [Bell, Peter A., et al. "Integrating knowledge of protein sequence with protein function for the prediction and validation of new MALT1 substrates." Computational and Structural Biotechnology Journal 20 (2022): 4717-4732.])
Then you can get sequence, structure, and GO, PPI information for any 16-kmers.

get_all_kmer_id.py: generate a dictionary with the key as protein ID, and values as 16-kmers of the protein. The input file is human.fasta containing all human proteins. (data provided by Konstantinos Kalogeropoulos)
get_pdb.py: download protein structures from AlphaFold.
main.py: generate positive data and negative data only with sequence and structure (seq+ss) information for a single protease predictor.

lost a file to generate GO information, it might be manually downloaded
lost a file to get PPI information.
But both two processed data can be obtained from 'data/inhouse_data/GlucGO' and 'data/inhouse_data/MMP9GO'


If you want to train model and test:
1) final_mian.py: run test set on trained models and generate results five times.
	To train and test on binary task(eg. MMP9 substrate or not), using main_seqss.py and main_seqss_test.py
	To train and test on multiple classification task(eg. MMP9 substrate, GluC substrate or none of them), using df_split.py and multiple_test.py
1.1.1) main_seqss.py: train the model for single protease; save the trained model and generate test data for binary predictor.
1.1.2) main_seqss_test.py: load two trained models and test with data, generate binary classification results.
1.2.1) df_split.py: train the model for single protease; save the trained model and generate test data for multiple protease predictors.
	  positive data defined with parameter 'protease' in the df_split file.
	  negative data is generated with non-MMP2 and non-GluC data.
	  the model trained with sequence, structure, and GO information.
1.2.2) multiple_test.py: load two trained models and test with data (MMP9-pos, GluC pos, and all-neg), generate multiple classification results.


