from sklearn.model_selection import KFold
from load_data import *

cancer_samples, sample_names, cancer_types = load_cancer_samples_data(path_to_cancer_samples)


# -----------------------------Scheme 1: Randomly divide all samples-----------------------------------#
samples = cancer_samples.index.values
kf1 = KFold(n_splits=5,  
            shuffle=True,  
            random_state=0) 

cv_split = {}
idx = 0
round_idx = 25
for train_index, test_index in kf1.split(range(len(samples))):
    split_samples = {"train_sample": samples[train_index],
                     "val_sample": samples[test_index]}
    
    f_save = open(path_to_split_data_idx + "Random_split_data/fold_" + str(round_idx) + '_index.pkl', 'wb')
    pickle.dump(split_samples, f_save)
    f_save.close()
    round_idx += 1


    kf2 = KFold(n_splits=5,  
                shuffle=True, 
                random_state=0)
    for CV_sub_train_index, CV_sub_val_index in kf2.split(train_index): 
        CV_train_index = train_index[CV_sub_train_index]
        CV_val_index = train_index[CV_sub_val_index]
        split_samples = {"train_sample": samples[CV_train_index],
                         "val_sample": samples[CV_val_index]}
        f_save = open(path_to_split_data_idx + "Random_split_data/fold_" + str(idx) + '_index.pkl', 'wb')
        pickle.dump(split_samples, f_save)
        f_save.close()
        idx += 1





# ---------------------------Scheme 2: Each cancer type was randomly divided and then combined-------------------------------#
cancer_cv_split = {}
for cancer in cancer_types:
    samples = cancer_samples[cancer_samples["cancer_type"] == cancer].index.values  

    kf1 = KFold(n_splits=5,  
                shuffle=True,  
                random_state=0)  
    cancer_cv_split[cancer] = {}
    idx = 0
    round_idx = 25
    for train_index, test_index in kf1.split(range(len(samples))):
        cancer_cv_split[cancer][round_idx] = {"train_sample": samples[train_index],
                                              "val_sample": samples[test_index]}

        round_idx += 1

        kf2 = KFold(n_splits=5,  
                    shuffle=True,  
                    random_state=0)
        for CV_sub_train_index, CV_sub_val_index in kf2.split(train_index):  
            CV_train_index = train_index[CV_sub_train_index]
            CV_val_index = train_index[CV_sub_val_index]
            cancer_cv_split[cancer][idx] = {"train_sample": samples[CV_train_index],
                                            "val_sample": samples[CV_val_index]}
            idx += 1


final_CV_split = {}
for idx in range(30):
    final_CV_split[idx] = {}
    train_sample_list = [list(cancer_cv_split[cancer][idx]["train_sample"]) for cancer in cancer_types]
    test_sample_list = [list(cancer_cv_split[cancer][idx]["val_sample"]) for cancer in cancer_types]
    train_sample_list = [x for j in train_sample_list for x in j]
    test_sample_list = [x for j in test_sample_list for x in j]

    split_samples = {"train_sample": np.array(train_sample_list),
                     "val_sample": np.array(test_sample_list)}

    f_save = open(path_to_split_data_idx + "Cancer_split_data/fold_" + str(idx) + '_index.pkl', 'wb')
    pickle.dump(split_samples, f_save)
    f_save.close()
