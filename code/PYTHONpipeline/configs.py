# Hyperparameters to evaluate: For each, hyperparameter, provide a LIST of options to evaluate.
# A paramter grid will be generated so all combinations of these hyperparameters will be evaluated.
hyperparams = {"epochs": [200],
               "inner_activation": ["relu"],
               "hidden_sizes_shared": [[100]],
               "hidden_sizes_separate": [[50, 10]],
               "dropout": [.3],
               "k_reg": [0.001, .01, 0.1],
               "learning_rate": [.0005, .001],
               "batch_size": [20, 50]}

# ########################################### DATA ############################################
root_path = "D:/MTLHot/"
path_to_exp_info = root_path + "data/TCGA_data/log.exp.matrix.csv"

path_to_phenotypic_info = root_path + 'data/TCGA_data/pan.cancer.phenotype.csv'

path_to_cancer_samples = root_path + "data/TCGA_data/cancer_samples.csv"

phenotypes = ["APM_out", "TCell_out", "IFNgamma_out", "PDL1_out"]

path_to_validation_data = root_path + "data/validation_data/"

path_to_split_data_idx = root_path + "data/split_dataset/"


# ##################### RESULTS ##################################

MTL_model_path = root_path + "result/model_training/MTL/"
Linear_models_res_path = root_path + "result/model_training/Linear_models/"
non_Linear_models_res_path = root_path + "result/model_training/non_Linear_models/"

model_evaluation_res_path = root_path + "result/model_evaluation/"

validate_res_path = root_path + "result/model_validation/"

path_to_embedding_cluster = root_path + "result/last_shared_layer_clustering/"

IG_save_path = root_path + "result/model_interpretation/"
ranked_features_path = root_path + "result/model_interpretation/ranked_features/"
