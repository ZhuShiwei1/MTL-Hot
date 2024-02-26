import os
os.environ["CUDA_VISIBLE_DEVICES"]="-1"  
from Train_CV_model import *


data_exp_tcga, exp_gene_names, exp_sample_names = load_feature_data(path_to_exp_info)
pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)


X = data_exp_tcga
Y_labels = pheno_score


# -------------------------------1. Add shared layer--------------------------------------#
print("Add shared layer")
hyperparams = {"epochs": [200],
               "inner_activation": ["relu"],
               "hidden_sizes_shared": [[200, 100]],
               "hidden_sizes_separate": [[50, 10]],
               "dropout": [.3],
               "k_reg": [0.001, .01, 0.1],
               "learning_rate": [.0005, .001],
               "batch_size": [20, 50]}

res_save_path = MTL_model_path+"Add_shared_layer/Random_split_CV/"
path_to_final_chosen_models = res_save_path + "final_models_chosen/"
if not os.path.isdir(res_save_path):
    os.makedirs(res_save_path)
if not os.path.isdir(path_to_final_chosen_models):
    os.makedirs(path_to_final_chosen_models)


hy_dict_list = list(ParameterGrid(hyperparams))
for fold_idx in range(30):
    Train_CV_models(X, Y_labels, sample_names, fold_idx, hy_dict_list, res_save_path,
                    path_to_split_data_idx=path_to_split_data_idx+"Random_split_data/")

best_fnames_in_split = Select_best_hyperparms_from_CV(res_save_path, performance_measure="MSE")


performances = {}
for round_idx in range(5):
    performances[round_idx] = Get_performance(res_save_path, performance_measure=["MSE", "MAE", "RMSE", "R2", "corr", "pearson_corr"],
                                              foldidx_range=[round_idx+25, round_idx+26],
                                              specific_fnames=[best_fnames_in_split[round_idx]])
print(performances)
pickle.dump(performances, open(path_to_final_chosen_models + "performances.p", "wb"))


# ------------------------------- 2.Add specific layer--------------------------------------#
print("Add specific layer")

hyperparams = {"epochs": [200],
               "inner_activation": ["relu"],
               "hidden_sizes_shared": [[100]],
               "hidden_sizes_separate": [[100, 50, 10]],
               "dropout": [.3],
               "k_reg": [0.001, .01, 0.1],
               "learning_rate": [.0005, .001],
               "batch_size": [20, 50]}


res_save_path = MTL_model_path+"Add_specific_layer/Random_split_CV/"
path_to_final_chosen_models = res_save_path + "final_models_chosen/"
if not os.path.isdir(res_save_path):
    os.makedirs(res_save_path)
if not os.path.isdir(path_to_final_chosen_models):
    os.makedirs(path_to_final_chosen_models)


hy_dict_list = list(ParameterGrid(hyperparams))
for fold_idx in range(30):
    Train_CV_models(X, Y_labels, sample_names, fold_idx, hy_dict_list, res_save_path,
                    path_to_split_data_idx=path_to_split_data_idx+"Random_split_data/")


best_fnames_in_split = Select_best_hyperparms_from_CV(res_save_path, performance_measure="MSE")

performances = {}
for round_idx in range(5):
    performances[round_idx] = Get_performance(res_save_path, performance_measure=["MSE", "MAE", "RMSE", "R2", "corr", "pearson_corr"],
                                              foldidx_range=[round_idx+25, round_idx+26],
                                              specific_fnames=[best_fnames_in_split[round_idx]])
print(performances)
pickle.dump(performances, open(path_to_final_chosen_models + "performances.p", "wb"))


# ------------------------------- 3.reduce specific layer --------------------------------------#
print("Reduce specific layer")
hyperparams = {"epochs": [200],
               "inner_activation": ["relu"],
               "hidden_sizes_shared": [[100]],
               "hidden_sizes_separate": [[10]],
               "dropout": [.3],
               "k_reg": [0.001, .01, 0.1],
               "learning_rate": [.0005, .001],
               "batch_size": [20, 50]}


res_save_path = MTL_model_path+"Reduce_specific_layer/Random_split_CV/"
path_to_final_chosen_models = res_save_path + "final_models_chosen/"
if not os.path.isdir(res_save_path):
    os.makedirs(res_save_path)
if not os.path.isdir(path_to_final_chosen_models):
    os.makedirs(path_to_final_chosen_models)


hy_dict_list = list(ParameterGrid(hyperparams))
for fold_idx in range(30):
    Train_CV_models(X, Y_labels, sample_names, fold_idx, hy_dict_list, res_save_path,
                    path_to_split_data_idx=path_to_split_data_idx+"Random_split_data/")
best_fnames_in_split = Select_best_hyperparms_from_CV(res_save_path, performance_measure="MSE")

performances = {}
for round_idx in range(5):
    performances[round_idx] = Get_performance(res_save_path, performance_measure=["MSE", "MAE", "RMSE", "R2", "corr", "pearson_corr"],
                                              foldidx_range=[round_idx+25, round_idx+26],
                                              specific_fnames=[best_fnames_in_split[round_idx]])
print(performances)
pickle.dump(performances, open(path_to_final_chosen_models + "performances.p", "wb"))


# ------------------------------- 4.Split late --------------------------------------#
print("Split late")
hyperparams = {"epochs": [200],
               "inner_activation": ["relu"],
               "hidden_sizes_shared": [[100, 50]],
               "hidden_sizes_separate": [[10]],
               "dropout": [.3],
               "k_reg": [0.001, .01, 0.1],
               "learning_rate": [.0005, .001],
               "batch_size": [20, 50]}


res_save_path = MTL_model_path+"Split_late/Random_split_CV/"
path_to_final_chosen_models = res_save_path + "final_models_chosen/"
if not os.path.isdir(res_save_path):
    os.makedirs(res_save_path)
if not os.path.isdir(path_to_final_chosen_models):
    os.makedirs(path_to_final_chosen_models)

hy_dict_list = list(ParameterGrid(hyperparams))
for fold_idx in range(30):
    Train_CV_models(X, Y_labels, sample_names, fold_idx, hy_dict_list, res_save_path,
                    path_to_split_data_idx=path_to_split_data_idx+"Random_split_data/")
best_fnames_in_split = Select_best_hyperparms_from_CV(res_save_path, performance_measure="MSE")


performances = {}
for round_idx in range(5):
    performances[round_idx] = Get_performance(res_save_path, performance_measure=["MSE", "MAE", "RMSE", "R2", "corr", "pearson_corr"],
                                              foldidx_range=[round_idx+25, round_idx+26],
                                              specific_fnames=[best_fnames_in_split[round_idx]])
print(performances)
pickle.dump(performances, open(path_to_final_chosen_models + "performances.p", "wb"))
