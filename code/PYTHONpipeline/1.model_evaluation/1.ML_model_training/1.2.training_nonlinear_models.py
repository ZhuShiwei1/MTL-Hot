import os
import time
from load_data import *
from single_output_model_performance import *


data_exp_tcga, exp_gene_names, exp_sample_names = load_feature_data(path_to_exp_info)
pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)


X = data_exp_tcga
Y_labels = pheno_score

# ---------------------------------------SVM_linear------------------------------------------------------- #
method = "SVM_linear"
path_to_final_chosen_models = non_Linear_models_res_path + "/%s/" % method
if not os.path.isdir(path_to_final_chosen_models):
    os.makedirs(path_to_final_chosen_models)


Random_split_performances = {}
Random_split_best_alphas = {}
Noise_performances = {}
Noise_best_alphas = {}
for round_idx in range(5):
    print(round_idx)

    Random_split_performances[round_idx], Random_split_best_alphas[round_idx] = \
        single_output_model_performance(method, X, Y_labels, sample_names,
                                        fold_idx=round_idx + 25,
                                        path_to_split_data_idx=path_to_split_data_idx + "Random_split_data/",
                                        PCA_dim=300)
    print(Random_split_performances[round_idx])

    Noise_performances[round_idx], Noise_best_alphas[round_idx] = \
        single_output_model_performance(method, X, Y_labels, sample_names,
                                        fold_idx=round_idx + 25,
                                        path_to_split_data_idx=path_to_split_data_idx + "Random_split_data/",
                                        PCA_dim=300,
                                        add_noise=0.01)
    print(Noise_performances[round_idx])

print("SVM linear  Random_split_performances:", Random_split_performances)
pickle.dump(Random_split_performances, open(path_to_final_chosen_models + "Random_split_performances.p", "wb"))
pickle.dump(Random_split_best_alphas, open(path_to_final_chosen_models + "Random_split_best_alphas.p", "wb"))

print("SVM linear  Noise_performances:", Noise_performances)
pickle.dump(Noise_performances, open(path_to_final_chosen_models + "Noise_performances.p", "wb"))


# ---------------------------------------SVM_RBF------------------------------------------------------- #
method = "SVM_RBF"
path_to_final_chosen_models = non_Linear_models_res_path + "/%s/" % method
if not os.path.isdir(path_to_final_chosen_models):
    os.makedirs(path_to_final_chosen_models)

gamma_list = [1e-3, 0.01, 0.1, 0.5]  
Random_split_performances = {}
Random_split_best_alphas = {}
Noise_performances = {}
Noise_best_alphas = {}
for round_idx in range(5):
    print(round_idx)
    Random_split_performances[round_idx], Random_split_best_alphas[round_idx] = \
        single_output_model_performance(method, X, Y_labels, sample_names,
                                        fold_idx=round_idx + 25,
                                        path_to_split_data_idx=path_to_split_data_idx + "Random_split_data/",
                                        gamma_list=gamma_list,
                                        PCA_dim=300)

    Noise_performances[round_idx], Noise_best_alphas[round_idx] = \
        single_output_model_performance(method, X, Y_labels, sample_names,
                                        fold_idx=round_idx + 25,
                                        path_to_split_data_idx=path_to_split_data_idx + "Random_split_data/",
                                        gamma_list=None,
                                        best_gamma=Random_split_best_alphas[round_idx],
                                        PCA_dim=300,
                                        add_noise=0.01)

print("SVM RBF Random_split_performances:", Random_split_performances)
pickle.dump(Random_split_performances, open(path_to_final_chosen_models + "Random_split_performances.p", "wb"))
pickle.dump(Random_split_best_alphas, open(path_to_final_chosen_models + "Random_split_best_alphas.p", "wb"))
print("SVM RBF Noise_performances:", Noise_performances)
pickle.dump(Noise_performances, open(path_to_final_chosen_models + "Noise_performances.p", "wb"))

# -----------------------------------RandomForest-----------------------------------------#
method = "RandomForest"
path_to_final_chosen_models = non_Linear_models_res_path + "/%s/" % method
if not os.path.isdir(path_to_final_chosen_models):
    os.makedirs(path_to_final_chosen_models)

RF_params = {'n_estimators': [100, 300, 500], 'max_features': ['log2', 'sqrt', None]}
Random_split_performances = {}
Random_split_best_alphas = {}
Noise_performances = {}
Noise_best_alphas = {}
for round_idx in range(5):
    print(round_idx)
    t1 = time.time()
    Random_split_performances[round_idx], Random_split_best_alphas[round_idx] = \
        single_output_model_performance(method, X, Y_labels, sample_names,
                                        fold_idx=round_idx + 25,
                                        path_to_split_data_idx=path_to_split_data_idx + "Random_split_data/",
                                        params=RF_params, best_params=None,
                                        PCA_dim=300)
    print('\n\nModel training completed in %.1f mins.' % ((time.time() - t1) / 60))
    
    Noise_performances[round_idx], Noise_best_alphas[round_idx] = \
        single_output_model_performance(method, X, Y_labels, sample_names,
                                        fold_idx=round_idx + 25,
                                        path_to_split_data_idx=path_to_split_data_idx + "Random_split_data/",
                                        params=None, best_params=Random_split_best_alphas[round_idx],
                                        PCA_dim=300,
                                        add_noise=0.01)

print("RandomForest Random_split_performances:", Random_split_performances)
pickle.dump(Random_split_performances, open(path_to_final_chosen_models + "Random_split_performances.p", "wb"))
pickle.dump(Random_split_best_alphas, open(path_to_final_chosen_models + "Random_split_best_alphas.p", "wb"))

print("RandomForest Noise_performances:", Noise_performances)
pickle.dump(Noise_performances, open(path_to_final_chosen_models + "Noise_performances.p", "wb"))

#########################################################################
# Identify optimal hyperparameters for all samples
normalized_X_tcga, Y_train = load_final_data(X, Y_labels, shuffle=False)
exp_pca = PCA(n_components=300)
exp_pca.fit(normalized_X_tcga)  
X_tcga_transformed = exp_pca.transform(normalized_X_tcga)

method = "SVM_RBF"
path_to_final_chosen_models = non_Linear_models_res_path + "/%s/" % method
if not os.path.isdir(path_to_final_chosen_models):
    os.makedirs(path_to_final_chosen_models)
best_parameters = {}
for i in range(4):
    var = phenotypes[i]
 
    Y_train_var = Y_train[var]
    print(var)
    model, parameters = single_output_model(method, X_tcga_transformed, Y_train_var,
                                            params=[1e-3, 0.01, 0.1, 0.5])
    best_parameters[var] = parameters
pickle.dump(best_parameters, open(path_to_final_chosen_models + "best_parameters.p", "wb"))


method = "RandomForest"
RF_params = {'n_estimators': [100, 300, 500], 'max_features': ['log2', 'sqrt', None]}

path_to_final_chosen_models = non_Linear_models_res_path + "/%s/" % method
if not os.path.isdir(path_to_final_chosen_models):
    os.makedirs(path_to_final_chosen_models)
best_parameters = {}
for i in range(4):
    var = phenotypes[i]
  
    Y_train_var = Y_train[var]
    print(var)
    model, parameters = single_output_model(method, X_tcga_transformed, Y_train_var,
                                            params=RF_params)
    best_parameters[var] = parameters
pickle.dump(best_parameters, open(path_to_final_chosen_models + "best_parameters.p", "wb"))

gc.collect()