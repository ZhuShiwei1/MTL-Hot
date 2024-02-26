
import datetime
import gc
import os
import random
from scipy import stats as ss
from sklearn.model_selection import ParameterGrid
from MTL_model import *



# Perform 5-fold cross-validation and select the optimal hyperparameters
def Train_CV_models(X, Y_labels, sample_names, fold_idx, hy_dict_list, CV_save_path, path_to_split_data_idx):

    path_to_results = CV_save_path + "CV_results/"
    phenotypes = ["APM_out", "TCell_out", "IFNgamma_out", "PDL1_out"]

    # --------------------------------------------data preparation---------------------------------------------#
    # split data
    X_train, X_valid, Y_train, Y_valid = load_data_for_fold(X, Y_labels, sample_names,
                                                            path_to_split_data_idx,
                                                            fold_idx=fold_idx)

    # PCA
    exp_pca = PCA(n_components=300)
    exp_pca.fit(X_train) 
    X_train_transformed = exp_pca.transform(X_train)
    X_valid_transformed = exp_pca.transform(X_valid)

    # --------------------------------------------training model-------------------------------------------------#
    for hy_iteration in range(len(hy_dict_list)):
        print(datetime.datetime.now())
        print("FOLD:", fold_idx)
        print("HYPERPARAMETER ITERATION: %d" % hy_iteration)

        hy_dict = hy_dict_list[hy_iteration]
        title = "%s_%s_%f_%f_%f_%d" % (str(hy_dict["hidden_sizes_shared"]),
                                       str(hy_dict["hidden_sizes_separate"]),
                                       hy_dict["dropout"],
                                       hy_dict["k_reg"],
                                       hy_dict["learning_rate"],
                                       hy_dict["batch_size"])
        print(title)

        res_dest = path_to_results + "/" + title + "/"
        if not os.path.isdir(res_dest):
            os.makedirs(res_dest)

        # construct model 
        input_size = X_train_transformed.shape[1]
        prediction_model = get_prediction_model(hy_dict, input_size)
        trainable_model = get_trainable_model(prediction_model, input_size)

        opt = adam_v2.Adam(learning_rate=hy_dict["learning_rate"])  
        trainable_model.compile(optimizer=opt, loss=None)

        # training model
        trainable_model.fit([X_train_transformed,
                             Y_train[phenotypes[0]], Y_train[phenotypes[1]],
                             Y_train[phenotypes[2]], Y_train[phenotypes[3]]],
                            epochs=hy_dict["epochs"],
                            batch_size=hy_dict["batch_size"],
                            shuffle=True,  
                            verbose=0,  
                            validation_data=([X_valid_transformed,
                                              Y_valid[phenotypes[0]], Y_valid[phenotypes[1]],
                                              Y_valid[phenotypes[2]], Y_valid[phenotypes[3]]],)
                            )
        loss_weight = [np.exp(K.get_value(log_var[0])) ** 0.5 for log_var in trainable_model.layers[-1].log_vars]
        print(loss_weight)  

        # model performance
        y_pred_all = prediction_model.predict(X_valid_transformed)
        performance = {}
        for i in range(4):
            var = phenotypes[i]
            y_true = Y_valid[var]
            y_pred = y_pred_all[i]
            mse = mean_squared_error(y_true, y_pred)
            MAE = mean_absolute_error(y_true, y_pred)
            RMSE = math.sqrt(mse)
            r2 = r2_score(y_true, y_pred)
            r, p = ss.spearmanr(y_pred, y_true)
            pearson_r, p1 = stats.pearsonr(y_pred.reshape(y_true.shape), y_true)
            performance[var] = {"MSE": mse, "MAE": MAE, "RMSE": RMSE, "R2": r2, "corr": r, "pearson_corr": pearson_r}

        pickle.dump(performance, open(res_dest + '%d.p' % fold_idx, 'wb'))
        K.clear_session()
        gc.collect()


#  select the optimal combination of hyperparameters from each round
def Select_best_hyperparms_from_CV(CV_save_path, performance_measure):

    path_to_results = CV_save_path + "CV_results/"

    path_to_final_chosen_models = CV_save_path + "final_models_chosen/"
    if not os.path.isdir(path_to_final_chosen_models):
        os.makedirs(path_to_final_chosen_models)

    phenotypes = ["APM_out", "TCell_out", "IFNgamma_out", "PDL1_out"]

    performances = Get_performance(CV_save_path, performance_measure, foldidx_range=[0, 25])

    # -----------------Get the best combination of hyperparameters for the 5-fold cross-validation filter in each round ---------------------------------#
    fnames = os.listdir(path_to_results + "/")
    performances_averages = {}
    CV_rankings = {}
    rank_averages = {}
    best_fnames_in_split = {}
    for round_idx in range(5):
        performances_averages[round_idx] = {}
        for var in phenotypes:
            performances_averages[round_idx][var] = {}
            for fname in fnames:
                performances_averages[round_idx][var][fname] = np.nanmean(
                    np.array([performances[var][fname][(5 * round_idx):(5 * round_idx + 5)]]))


        CV_rankings[round_idx] = {}
        for var in phenotypes:
            fnames_rank = ss.rankdata(list(performances_averages[round_idx][var].values()))
            i = 0
            CV_rankings[round_idx][var] = {}
            for fname in performances_averages[round_idx][var].keys():
                CV_rankings[round_idx][var][fname] = fnames_rank[i]
                i += 1

        
        rank_averages[round_idx] = {}
        for fname in fnames:
            rank_averages[round_idx][fname] = np.nanmean(
                np.array([CV_rankings[round_idx][var][fname] for var in phenotypes]))  # 排秩取均值

        best_fnames_in_split[round_idx] = min(rank_averages[round_idx], key=lambda k: rank_averages[round_idx][k])

    pickle.dump(best_fnames_in_split, open(path_to_final_chosen_models + "folds.p", "wb"))

    return best_fnames_in_split


def Get_performance(CV_save_path, performance_measure, foldidx_range, specific_fnames=None):
    path_to_results = CV_save_path + "CV_results/"

    performances = {}
    for var in phenotypes:
        performances[var] = {}
        if specific_fnames == None:
            fnames = os.listdir(path_to_results + "/")
        else:
            fnames = specific_fnames

        for fname in fnames:
            if isinstance(performance_measure, str): 
                performances[var][fname] = []
            else:
                performances[var][fname] = {}

            for i in range(foldidx_range[0], foldidx_range[1]):
                per_performance = pickle.load(open((path_to_results + fname + '/%d.p' % i), 'rb'))
                if isinstance(performance_measure, str):
                    performances[var][fname].append(per_performance[var][performance_measure])
                else:
                    for measure in performance_measure:
                        performances[var][fname][measure] = []
                        performances[var][fname][measure].append(per_performance[var][measure])

    return performances


# Get the best combination of hyperparameters for final model
def Select_best_hyperparms_for_final(CV_save_path, performance_measure):
    path_to_results = CV_save_path + "CV_results/"

    path_to_final_chosen_models = CV_save_path + "final_models_chosen/"

    performances = Get_performance(CV_save_path, performance_measure, foldidx_range=[25, 30])

    fname_performances_averages = {}
    fnames = os.listdir(path_to_results + "/")
    for var in phenotypes:
        fname_performances_averages[var] = {}
        for fname in fnames:

            fname_performances_averages[var][fname] = np.nanmean(np.array([performances[var][fname][0:5]]))

    fname_rankings = {}
    for var in phenotypes:
        fnames_rank = ss.rankdata(list(fname_performances_averages[var].values()))
        i = 0
        fname_rankings[var] = {}
        for fname in fnames:
            fname_rankings[var][fname] = fnames_rank[i]
            i += 1
    fname_rank_averages = {}
    for fname in fnames:
        fname_rank_averages[fname] = np.nanmean(
            np.array([fname_rankings[var][fname] for var in phenotypes])) 

    final_best_fname = min(fname_rank_averages, key=lambda k: fname_rank_averages[k])


    pickle.dump(final_best_fname, open(path_to_final_chosen_models + "/final.p", "wb"))
    return final_best_fname

