import os
from load_data import *
from single_output_model_performance import *

data_exp_tcga, exp_gene_names, exp_sample_names = load_feature_data(path_to_exp_info)
pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)

X = data_exp_tcga
Y_labels = pheno_score


# ------------------------------ ridge ------------------------------------#
method = "ridge"
path_to_final_chosen_models = Linear_models_res_path + "/%s/" % method
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
                                        alphas=np.linspace(0.001, 10, 100),  
                                        PCA_dim=300)
    print(Random_split_performances[round_idx])

    Noise_performances[round_idx], Noise_best_alphas[round_idx] = \
        single_output_model_performance(method, X, Y_labels, sample_names,
                                        fold_idx=round_idx + 25,
                                        path_to_split_data_idx=path_to_split_data_idx + "Random_split_data/",
                                        best_alpha=Random_split_best_alphas[round_idx],
                                        PCA_dim=300,
                                        add_noise=0.01)
    print(Noise_performances[round_idx])


print("Ridge Random_split_performances:", Random_split_performances)
pickle.dump(Random_split_performances, open(path_to_final_chosen_models + "Random_split_performances.p", "wb"))
pickle.dump(Random_split_best_alphas, open(path_to_final_chosen_models + "Random_split_best_alphas.p", "wb"))
print("Ridge Noise_performances:", Noise_performances)
pickle.dump(Noise_performances, open(path_to_final_chosen_models + "Noise_performances.p", "wb"))


# ------------------------------ lasso ------------------------------------#
method = "lasso"

path_to_final_chosen_models = Linear_models_res_path + "/%s/" % method
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
                                        alphas=None,  
                                        PCA_dim=300)
    print(Random_split_performances[round_idx])

    Noise_performances[round_idx], Noise_best_alphas[round_idx] = \
        single_output_model_performance(method, X, Y_labels, sample_names,
                                        fold_idx=round_idx + 25,
                                        path_to_split_data_idx=path_to_split_data_idx + "Random_split_data/",
                                        best_alpha=Random_split_best_alphas[round_idx],
                                        PCA_dim=300,
                                        add_noise=0.01)
    print(Noise_performances[round_idx])


print("Lasso Random_split_performances:", Random_split_performances)
pickle.dump(Random_split_performances, open(path_to_final_chosen_models + "Random_split_performances.p", "wb"))
pickle.dump(Random_split_best_alphas, open(path_to_final_chosen_models + "Random_split_best_alphas.p", "wb"))
print("Lasso Noise_performances:", Noise_performances)
pickle.dump(Noise_performances, open(path_to_final_chosen_models + "Noise_performances.p", "wb"))


# ------------------------------ ElasticNet ------------------------------------#
method = "ElasticNet"

path_to_final_chosen_models = Linear_models_res_path + "/%s/" % method
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
                                        alphas=None, 
                                        PCA_dim=300)
    print(Random_split_performances[round_idx])
    Noise_performances[round_idx], Noise_best_alphas[round_idx] = \
        single_output_model_performance(method, X, Y_labels, sample_names,
                                        fold_idx=round_idx + 25,
                                        path_to_split_data_idx=path_to_split_data_idx + "Random_split_data/",
                                        best_alpha=Random_split_best_alphas[round_idx],
                                        PCA_dim=300,
                                        add_noise=0.01)
    print(Noise_performances[round_idx])

print("ElasticNet Random_split_performances:", Random_split_performances)
pickle.dump(Random_split_performances, open(path_to_final_chosen_models + "Random_split_performances.p", "wb"))
pickle.dump(Random_split_best_alphas, open(path_to_final_chosen_models + "Random_split_best_alphas.p", "wb"))
print("ElasticNet Noise_performances:", Noise_performances)
pickle.dump(Noise_performances, open(path_to_final_chosen_models + "Noise_performances.p", "wb"))


#########################################################################
# Identify optimal hyperparameters for all samples
normalized_X_tcga, Y_train = load_final_data(X, Y_labels, shuffle=False)
exp_pca = PCA(n_components=300)
exp_pca.fit(normalized_X_tcga)  
X_tcga_transformed = exp_pca.transform(normalized_X_tcga)
for method in ["ridge", "lasso", "ElasticNet"]:
    path_to_final_chosen_models = Linear_models_res_path + "/%s/" % method
    if not os.path.isdir(path_to_final_chosen_models):
        os.makedirs(path_to_final_chosen_models)
    best_parameters = {}
    for i in range(4):
        var = phenotypes[i]
        Y_train_var = Y_train[var]
        print(var)

        if method == "ridge":
            model, parameters = single_output_model(method, X_tcga_transformed, Y_train_var,
                                                    alphas=np.linspace(0.001, 10, 100))
        elif method == "lasso":
            model, parameters = single_output_model(method, X_tcga_transformed, Y_train_var, alphas=None)

        elif method == "ElasticNet":
            model, parameters = single_output_model(method, X_tcga_transformed, Y_train_var, alphas=None)

        best_parameters[var] = parameters

    pickle.dump(best_parameters, open(path_to_final_chosen_models + "best_parameters.p", "wb"))

gc.collect()
