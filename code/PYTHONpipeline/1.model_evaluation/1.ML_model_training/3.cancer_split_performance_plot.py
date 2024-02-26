import scipy
from matplotlib.backends.backend_pdf import PdfPages
from Train_CV_model import *


data_exp_tcga, exp_gene_names, exp_sample_names = load_feature_data(path_to_exp_info)
pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)

X = data_exp_tcga
Y_labels = pheno_score

if not os.path.isdir(model_evaluation_res_path + 'model_comparison/cancer_split/'):
    os.makedirs(model_evaluation_res_path + 'model_comparison/cancer_split/')

# -------------------------------------performance for the splited data of Scheme 2(Each cancer type was randomly divided and then combined)------------------------------#
CV_save_path = MTL_model_path + "MTL_Hot/CV/Cancer_split_CV/"

best_fnames_in_split = pickle.load(open(CV_save_path + "final_models_chosen/folds.p", "rb"))

phenotypes = ["APM_out", "TCell_out", "IFNgamma_out", "PDL1_out"]

performances = {}
for fold_idx in range(25, 30):
   
    X_train, X_valid, Y_train, Y_valid = load_data_for_fold(X, Y_labels, sample_names,
                                                            path_to_split_data_idx=path_to_split_data_idx+"Cancer_split_data/",
                                                            fold_idx=fold_idx)
    
    exp_pca = PCA(n_components=300)
    exp_pca.fit(X_train)  

    X_train_transformed = exp_pca.transform(X_train)
    X_valid_transformed = exp_pca.transform(X_valid)


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

    opt = adam_v2.Adam(learning_rate=MDNN_hyperparams["learning_rate"])  
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
    with PdfPages(model_evaluation_res_path + 'model_comparison/cancer_split/Cancer_split_round_%d_corr_plot.pdf' % (fold_idx-25)) as pdf:
        for i in range(4):
            var = phenotypes[i]
            y_true = Y_valid[var]
            y_pred = y_pred_all[i]
            plt.figure(figsize=(4, 4))

            plt.scatter(y_pred, y_true, color="blue")
            plt.title(var+', samples=%d' % len(y_pred))
            plt.ylabel('Original labels')
            plt.xlabel('Predicted labels')
            plt.grid()

            lrModel = linear_model.LinearRegression()
            lrModel.fit(y_pred, y_true)
            deter_coef = lrModel.score(y_pred, y_true)
            plt.plot(y_pred, lrModel.predict(y_pred), color='red')  
            pdf.savefig()
            plt.close()

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
            path_to_output = CV_save_path + "rounds_best_hyperparms_outputs/%s/" % str(fold_idx - 25)
            if not os.path.isdir(path_to_output):
                os.makedirs(path_to_output)
            np.savetxt(path_to_output + var + "_predicts.txt", y_pred)

    performances[(fold_idx - 25)] = performance
    pickle.dump(performances, open(CV_save_path + "rounds_best_hyperparms_outputs/" + "performances.p", "wb"))

performances = pickle.load(open(CV_save_path + "rounds_best_hyperparms_outputs/" + "performances.p", "rb"))
print(performances)

reformat_performance = {}
for var in phenotypes:
    reformat_performance[var] = []
    for fold_idx in range(5):
        reformat_performance[var].append(performances[fold_idx][var]["pearson_corr"])
reformat_performance_mean = []
for var in phenotypes:
    reformat_performance_mean.append(np.mean(reformat_performance[var]))
print(reformat_performance)


with PdfPages(model_evaluation_res_path + 'model_comparison/cancer_split/Cancer_split_corr_bar_plot.pdf') as pdf:
    
    labels = phenotypes
    plt.figure(figsize=(10, 10))
    x = np.arange(4)  # the label locations
    plt.bar(labels, reformat_performance_mean, alpha=0.7)
   
    scatter_x_labels = np.repeat(labels, 5, axis=0)
    scatter_y = np.array([list(item) for item in reformat_performance.values()])
    plt.scatter(scatter_x_labels, scatter_y, color='black', s=10)
    plt.ylim(0.86, 1)
    plt.xticks(size=16, rotation=20)  
    plt.yticks(size=16)
    plt.ylabel("corr", fontsize=20)  
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()


# -------------------------------------performance for the splited data of Scheme 1(Randomly divide all samples)------------------------------#
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
    with PdfPages(model_evaluation_res_path + 'model_comparison/cancer_split/Random_split_round_%d_corr_plot.pdf' % (fold_idx-25)) as pdf:
        for i in range(4):
            var = phenotypes[i]
            y_true = Y_valid[var]
            y_pred = y_pred_all[i]
            plt.figure(figsize=(4, 4))
            plt.scatter(y_pred, y_true, color="blue")
            plt.title(var+', samples=%d' % len(y_pred))
            plt.ylabel('Original labels')
            plt.xlabel('Predicted labels')
            plt.grid()

            lrModel = linear_model.LinearRegression()
            lrModel.fit(y_pred, y_true)
            deter_coef = lrModel.score(y_pred, y_true)
            plt.plot(y_pred, lrModel.predict(y_pred), color='red')  
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()
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
            path_to_output = CV_save_path + "rounds_best_hyperparms_outputs/%s/" % str(fold_idx-25)
            if not os.path.isdir(path_to_output):
                os.makedirs(path_to_output)
            np.savetxt(path_to_output + var + "_predicts.txt", y_pred)

    performances[(fold_idx-25)] = performance

    pickle.dump(performances, open(CV_save_path + "rounds_best_hyperparms_outputs/" + "performances.p", "wb"))


performances = pickle.load(open(CV_save_path + "rounds_best_hyperparms_outputs/" + "performances.p", "rb"))
print(performances)
