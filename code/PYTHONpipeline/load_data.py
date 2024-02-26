import pickle
import numpy as np
import pandas as pd
from configs import *
from sklearn import preprocessing


# import expression data
def load_feature_data(filename):
    data = pd.read_csv(filename, sep=',', index_col=0)
    features = data.columns.values
    sample_names = data.index.values
    data = data.to_numpy()  
    return data, features, sample_names


# load The phenotypes' scores of the samples
def load_phenotype_data(filename):
    pheno_score = pd.read_csv(filename, sep=',', index_col=0)
    pheno_labels = pheno_score.columns.values
    sample_names = pheno_score.index.values
    return pheno_score, pheno_labels, sample_names


# import the corresponding cancer type of each sample
def load_cancer_samples_data(filename):
    cancer_samples = pd.read_csv(filename, sep=',')
    sample_names = cancer_samples.index.values  
    cancer_types = cancer_samples['cancer_type'].tolist()
    cancer_types = list({}.fromkeys(cancer_types).keys())  
    return cancer_samples, sample_names, cancer_types


# Import the splited data set and shuffle the sample labels
def load_data_for_fold(X, Y_labels, sample_names, path_to_split_data_idx, fold_idx):
    np.random.seed(0)
    f_read = open(path_to_split_data_idx + "/fold_" + str(fold_idx) + '_index.pkl', 'rb')
    split_samples = pickle.load(f_read)
    train_index = np.flatnonzero(np.in1d(sample_names, split_samples["train_sample"]))
    val_index = np.flatnonzero(np.in1d(sample_names, split_samples["val_sample"]))
    X_train = X[train_index].astype(np.float32)
    X_valid = X[val_index].astype(np.float32)

    labels_train = Y_labels.iloc[train_index, :]
    labels_valid = Y_labels.iloc[val_index, :]

    # The input data of training set and validation set are normalized respectively
    min_max_scaler = preprocessing.MinMaxScaler()
    normalized_X_train = min_max_scaler.fit_transform(X_train)
    normalized_X_valid = min_max_scaler.fit_transform(X_valid)

    # shuffle sample order
    shuffle_idx_train = np.random.permutation(range(X_train.shape[0]))
    shuffle_idx_valid = np.random.permutation(range(X_valid.shape[0]))
    X_train_shuffle = normalized_X_train[shuffle_idx_train]
    X_valid_shuffle = normalized_X_valid[shuffle_idx_valid]

    labels_train_shuffle = labels_train.iloc[shuffle_idx_train, :]
    labels_valid_shuffle = labels_valid.iloc[shuffle_idx_valid, :]

    phenotypes = labels_train.columns.values 
    y_train = {}
    y_valid = {}
    for phen in phenotypes:
        y_train[phen] = labels_train_shuffle[phen].astype(float).values
        y_valid[phen] = labels_valid_shuffle[phen].astype(float).values

    return X_train_shuffle, X_valid_shuffle, y_train, y_valid



def load_final_data(X, Y_labels, shuffle=True):
    np.random.seed(0)
    # All the features are normalized
    min_max_scaler = preprocessing.MinMaxScaler()
    normalized_X = min_max_scaler.fit_transform(X)

    if shuffle is True:
        #shuffle sample order
        shuffle_idx = np.random.permutation(range(X.shape[0]))
        X_shuffle = normalized_X[shuffle_idx]
        Y_labels_shuffle = Y_labels.iloc[shuffle_idx, :]

        phenotypes = Y_labels_shuffle.columns.values 
        Y = {}
        for phen in phenotypes:
            Y[phen] = Y_labels_shuffle[phen].astype(float).values

        return X_shuffle, Y, shuffle_idx
    else:
        phenotypes = Y_labels.columns.values  
        Y = {}
        for phen in phenotypes:
            Y[phen] = Y_labels[phen].astype(float).values
        return normalized_X, Y

# normalized the validation data
def normalized_val_data(val_feature_data, train_feature_data):
    scaled_val_exp_data = np.zeros((val_feature_data.shape[0], val_feature_data.shape[1]))
    if val_feature_data.shape[1] == train_feature_data.shape[1]:
        for i in range(train_feature_data.shape[1]):
            train_values = train_feature_data[:, i]
            val_values = val_feature_data[:, i]
            train_mean = np.mean(train_values)
            train_sd = np.std(train_values)
            val_mean = np.mean(val_values)
            val_sd = np.std(val_values)
            rescaled_value = ((val_values - val_mean) * train_sd / val_sd) + train_mean

            scaled_val_exp_data[:, i] = rescaled_value

        return scaled_val_exp_data
    else:
        print("Different feature numbers!")
        return val_feature_data
