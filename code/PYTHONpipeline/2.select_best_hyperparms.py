from Train_CV_model import *


data_exp_tcga, exp_gene_names, exp_sample_names = load_feature_data(path_to_exp_info)
pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)

X = data_exp_tcga
Y_labels = pheno_score

hy_dict_list = list(ParameterGrid(hyperparams)) 

# ------------------------- 5-fold CV in the splited data from Scheme 1(Randomly divide all samples) -------------------------------#
CV_save_path = MTL_model_path + "MTL_Hot/CV/Random_split_CV/"
path_to_final_chosen_models = CV_save_path + "final_models_chosen/"
if not os.path.isdir(CV_save_path):
    os.makedirs(CV_save_path)
if not os.path.isdir(path_to_final_chosen_models):
    os.makedirs(path_to_final_chosen_models)

#  5-fold CV
for fold_idx in range(30):
    Train_CV_models(X, Y_labels, sample_names, fold_idx, hy_dict_list, CV_save_path,
                    path_to_split_data_idx=path_to_split_data_idx + "Random_split_data/")

# Select best hyperparms for each round
best_fnames_in_split = Select_best_hyperparms_from_CV(CV_save_path, performance_measure="MSE")

print("Best hyper-parameter in each test:")
for key, value in best_fnames_in_split.items():
    print('{key}:{value}'.format(key=key, value=value))

performances = {}
for round_idx in range(5):
    performances[round_idx] = Get_performance(CV_save_path, performance_measure=["MSE", "MAE", "RMSE", "R2", "corr", "pearson_corr"],
                                              foldidx_range=[round_idx + 25, round_idx + 26],
                                              specific_fnames=[best_fnames_in_split[round_idx]]
                                              )
print(performances)
pickle.dump(performances, open(path_to_final_chosen_models + "Random_split_performances.p", "wb"))

final_best_fname = Select_best_hyperparms_for_final(CV_save_path, performance_measure="MSE")
print("Best hyper-parameter in final model:%s" % final_best_fname)


# -------------------------  5-fold CV in the splited data from Scheme 2 (Each cancer type was randomly divided and then combined)-------------------------------#
CV_save_path = MTL_model_path + "MTL_Hot/CV/Cancer_split_CV/"
path_to_final_chosen_models = CV_save_path + "final_models_chosen/"
if not os.path.isdir(CV_save_path):
    os.makedirs(CV_save_path)
if not os.path.isdir(path_to_final_chosen_models):
    os.makedirs(path_to_final_chosen_models)

# 5-fold CV
for fold_idx in range(30):
    Train_CV_models(X, Y_labels, sample_names, fold_idx, hy_dict_list, CV_save_path,
                    path_to_split_data_idx=path_to_split_data_idx + "Cancer_split_data/")

# Select best hyperparms for each round
best_fnames_in_split = Select_best_hyperparms_from_CV(CV_save_path, performance_measure="MSE")


print("Best hyper-parameter in each test:")
for key, value in best_fnames_in_split.items():
    print('{key}:{value}'.format(key=key, value=value))


performances = {}
for round_idx in range(5):
    performances[round_idx] = Get_performance(CV_save_path, performance_measure=["MSE", "MAE", "RMSE", "R2", "corr", "pearson_corr"],
                                              foldidx_range=[round_idx + 25, round_idx + 26],
                                              specific_fnames=[best_fnames_in_split[round_idx]])
print(performances)
pickle.dump(performances, open(path_to_final_chosen_models + "Cancer_split_performances.p", "wb"))
