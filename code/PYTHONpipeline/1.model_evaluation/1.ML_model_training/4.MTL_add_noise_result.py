import scipy
from Train_CV_model import *


data_exp_tcga, exp_gene_names, exp_sample_names = load_feature_data(path_to_exp_info)
pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)

X = data_exp_tcga
Y_labels = pheno_score


# -------------------------------------After adding noise to the training set------------------------------#

CV_save_path = MTL_model_path + "MTL_Hot/CV/Random_split_CV/"

best_fnames_in_split = pickle.load(open(CV_save_path + "final_models_chosen/folds.p", "rb"))

phenotypes = ["APM_out", "TCell_out", "IFNgamma_out", "PDL1_out"]
performances = {}
for fold_idx in range(25, 30):

    X_train, X_valid, Y_train, Y_valid = load_data_for_fold(X, Y_labels, sample_names,
                                                            path_to_split_data_idx=path_to_split_data_idx+"Random_split_data/",
                                                            fold_idx=fold_idx)


    exp_pca = PCA(n_components=300)
    exp_pca.fit(X_train) 

    X_train_transformed = exp_pca.transform(X_train)
    X_valid_transformed = exp_pca.transform(X_valid)

    # Adds noise to the labels of the training set
    # Perturbed 1% of the sample labels
    for phen in phenotypes:
        labels = Y_train[phen].astype(float).tolist()
        random_samp_num = round(len(labels) * 0.01)  
        random_idx = random.sample(range(len(labels)), random_samp_num)  
        random_shuff_idx = [i for i in random_idx]
        random.shuffle(random_shuff_idx)  
        for i in range(random_samp_num):
            labels[random_idx[i]] = labels[random_shuff_idx[i]]  
        Y_train[phen] = np.array(labels)


    hy_dict = best_fnames_in_split[(fold_idx-25)]
    print(hy_dict)

    hyperparams_values = hy_dict.split("_")
    MDNN_hyperparams = {"epochs": 200,
                        "inner_activation": "relu",
                        "hidden_sizes_shared": eval(hyperparams_values[0]),
                        "hidden_sizes_separate": eval(hyperparams_values[1]),
                        "dropout": float(hyperparams_values[2]),
                        "k_reg": float(hyperparams_values[3]),
                        "learning_rate": float(hyperparams_values[4]),
                        "batch_size": int(hyperparams_values[5])}

    input_size = X_train_transformed.shape[1]
    prediction_model = get_prediction_model(MDNN_hyperparams, input_size)
    trainable_model = get_trainable_model(prediction_model, input_size)

    opt = adam_v2.Adam(learning_rate=MDNN_hyperparams["learning_rate"])  # 学习率不宜设置过大
    trainable_model.compile(optimizer=opt, loss=None)


    trainable_model.fit([X_train_transformed,
                         Y_train[phenotypes[0]], Y_train[phenotypes[1]],
                         Y_train[phenotypes[2]], Y_train[phenotypes[3]]],
                        epochs=MDNN_hyperparams["epochs"],
                        batch_size=MDNN_hyperparams["batch_size"],
                        shuffle=True,  
                        verbose=0,  
                        validation_data=([X_valid_transformed,
                                          Y_valid[phenotypes[0]], Y_valid[phenotypes[1]],
                                          Y_valid[phenotypes[2]], Y_valid[phenotypes[3]]],)
                        )
    loss_weight = [np.exp(K.get_value(log_var[0])) ** 0.5 for log_var in trainable_model.layers[-1].log_vars]

    y_pred_all = prediction_model.predict(X_valid_transformed)
    performance = {}
    for i in range(4):
        var = phenotypes[i]
        y_true = Y_valid[var]
        y_pred = y_pred_all[i]
        mse = mean_squared_error(y_true, y_pred)
        corr, p1 = scipy.stats.spearmanr(y_pred, y_true)
        pearson_corr, p2 = scipy.stats.pearsonr(y_pred.reshape(y_true.shape), y_true)
        rmse = math.sqrt(mse)
        MAE = mean_absolute_error(y_true, y_pred)
        r2 = r2_score(y_true, y_pred)
        performance[var] = {"MSE": mse,
                            "MAE": MAE,
                            "RMSE": rmse,
                            "R2": r2,
                            "corr": corr,
                            "pearson_corr": pearson_corr}
        path_to_output = MTL_model_path + "MTL_Hot/CV/Add_noise_CV/rounds_best_hyperparms_outputs/%s/" % str(fold_idx-25)
        if not os.path.isdir(path_to_output):
            os.makedirs(path_to_output)
        np.savetxt(path_to_output + var + "_predicts.txt", y_pred)

    performances[(fold_idx-25)] = performance

    pickle.dump(performances, open(MTL_model_path + "MTL_Hot/CV/Add_noise_CV/rounds_best_hyperparms_outputs/" + "performances.p", "wb"))
