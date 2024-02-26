import gc
import math
import random
import numpy as np
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV, Ridge, ElasticNetCV, ElasticNet, LassoCV, Lasso
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, make_scorer
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVR
from load_data import load_data_for_fold
from configs import *


def get_ridge_model(X, Y, alphas, best_alpha=None, cv=5):
    if best_alpha is None:
        ridge_model = RidgeCV(alphas=alphas, cv=cv, scoring="neg_mean_squared_error")
        ridge_model.fit(X, Y)
        alpha = ridge_model.alpha_
    else:
        alpha = best_alpha

    model = Ridge(alpha=alpha)
    model.fit(X, Y)

    return model, alpha


def get_lasso_model(X, Y, alphas=None, best_alpha=None, cv=5):
    if best_alpha is None:
        lasso_model = LassoCV(alphas=alphas, cv=cv, random_state=1)
        lasso_model.fit(X, Y)
        alpha = lasso_model.alpha_
    else:
        alpha = best_alpha

    model = Lasso(alpha=alpha)
    model.fit(X, Y)

    return model, alpha


def get_ElasticNet_model(X, Y, alphas=None, best_alpha=None, cv=5):
    if best_alpha is None:
        enet_cv_model = ElasticNetCV(l1_ratio=.5, alphas=alphas, cv=cv, random_state=0)
        enet_cv_model.fit(X, Y)  
        alpha = enet_cv_model.alpha_
    else:
        alpha = best_alpha

    model = ElasticNet(l1_ratio=.5, alpha=alpha)
    model.fit(X, Y)

    return model, alpha


def get_SVM_model(X, Y, Kernel="linear", gamma_list=None, best_gamma=None, cv=5):
    if Kernel == "linear":
        linear_svr = SVR(kernel="linear")
        linear_svr.fit(X, Y)
        return linear_svr
    elif Kernel == "RBF":
        if best_gamma is None:
            # Parameters for tuning
            scorer = make_scorer(mean_squared_error, greater_is_better=False)
            parameters = [{'kernel': ['rbf'], 'gamma': gamma_list, 'C': [1]}]  # 调gamma
            RBF_svr_model = GridSearchCV(SVR(), parameters, cv=cv, scoring=scorer)
            RBF_svr_model.fit(X, Y)
            gamma = RBF_svr_model.best_params_['gamma']
        else:
            gamma = best_gamma
        RBF_svr = SVR(kernel="rbf", gamma=gamma, C=1)
        RBF_svr.fit(X, Y)
        best_parameters = gamma
        return RBF_svr, best_parameters


def get_RandomForest_model(X, Y, RF_params=None, best_RF_params=None, cv=5):

    if best_RF_params is None:
        # Parameters for tuning
        gscv = GridSearchCV(estimator=RandomForestRegressor(random_state=42, n_jobs=-1), param_grid=RF_params,
                            scoring='neg_mean_squared_error', cv=cv)
        gscv.fit(X, Y)
        best_parameters = gscv.best_params_
    else:
        best_parameters = best_RF_params

    RandomForest_model = RandomForestRegressor(random_state=42, n_jobs=-1,
                                               n_estimators=best_parameters['n_estimators'],
                                               max_features=best_parameters['max_features'])
    RandomForest_model.fit(X, Y)
    return RandomForest_model, best_parameters


def single_output_model(method, X_train_var, Y_train_var, params=None, best_params=None):

    global model, parameters
    if method == "ridge":
        if best_params is not None:
            best_alpha_value = best_params
            model, parameters = get_ridge_model(X_train_var, Y_train_var, alphas=None,
                                                best_alpha=best_alpha_value, cv=5)
        else:
            model, parameters = get_ridge_model(X_train_var, Y_train_var, alphas=params, cv=5)

    elif method == "lasso":
        if best_params is not None:
            best_alpha_value = best_params
            model, parameters = get_lasso_model(X_train_var, Y_train_var, alphas=None,
                                                best_alpha=best_alpha_value, cv=5)
        else:
            model, parameters = get_lasso_model(X_train_var, Y_train_var, alphas=params, cv=5)

    elif method == "ElasticNet":
        if best_params is not None:
            best_alpha_value = best_params
            model, parameters = get_ElasticNet_model(X_train_var, Y_train_var, alphas=None,
                                                     best_alpha=best_alpha_value, cv=5)
        else:
            model, parameters = get_ElasticNet_model(X_train_var, Y_train_var, alphas=params, cv=5)

    elif method == "SVM_linear":
        model = get_SVM_model(X_train_var, Y_train_var)  
        parameters = None

    elif method == "SVM_RBF":
        if best_params is not None:
            best_gamma_value = best_params
            model, parameters = get_SVM_model(X_train_var, Y_train_var, Kernel="RBF", gamma_list=None,
                                              best_gamma=best_gamma_value, cv=5)
        else:
            model, parameters = get_SVM_model(X_train_var, Y_train_var, Kernel="RBF", gamma_list=params, cv=5)
    elif method == "RandomForest":
        if best_params is not None:
            model, parameters = get_RandomForest_model(X_train_var, Y_train_var, best_RF_params=best_params, cv=5)
        else:
            model, parameters = get_RandomForest_model(X_train_var, Y_train_var, RF_params=params, cv=5)

    return model, parameters


def single_output_model_performance(method, X, Y_labels, sample_names, fold_idx, path_to_split_data_idx,
                                    params=None, best_params=None, PCA_dim=None, add_noise=None):
    # ------------------------------------------- data preparation ------------------------------------------------- #
    X_train, X_valid, Y_train, Y_valid = load_data_for_fold(X, Y_labels, sample_names,
                                                            path_to_split_data_idx,
                                                            fold_idx=fold_idx)
    if PCA_dim is not None:
        exp_pca = PCA(n_components=PCA_dim)
        exp_pca.fit(X_train)  
        X_train_transformed = exp_pca.transform(X_train)
        X_valid_transformed = exp_pca.transform(X_valid)
        X_train = X_train_transformed
        X_valid = X_valid_transformed

    # 如果为训练集的标签增加噪音
    if add_noise is not None:
        # 对1%的样本标签进行扰动
        for phen in phenotypes:
            labels = Y_train[phen].astype(float).tolist()
            random_samp_num = round(len(labels) * add_noise)  
            random_idx = random.sample(range(len(labels)), random_samp_num)  
            random_shuff_idx = [i for i in random_idx]
            random.shuffle(random_shuff_idx)
            for i in range(random_samp_num):
                labels[random_idx[i]] = labels[random_shuff_idx[i]]
            Y_train[phen] = np.array(labels)

    # -------------------------------- Build and train ML models ------------------------------------------#
    performances = {}
    best_parameters = {}
    for var in list(Y_train.keys()):

        Y_train_var = Y_train[var]
        Y_valid_var = Y_valid[var]
        X_train_var = X_train
        X_valid_var = X_valid
        print(var)

        train_nan_idx = np.argwhere(np.isnan(Y_train_var))
        valid_nan_idx = np.argwhere(np.isnan(Y_valid_var))
        if len(train_nan_idx) > 0:
            Y_train_var = np.delete(Y_train_var, train_nan_idx)
            X_train_var = np.delete(X_train, train_nan_idx, 0)
        if len(valid_nan_idx) > 0:
            Y_valid_var = np.delete(Y_valid_var, valid_nan_idx)
            X_valid_var = np.delete(X_valid, valid_nan_idx, 0)

        if method in ["ridge", "lasso", "ElasticNet", "SVM_linear", "SVM_RBF", "RandomForest"]:
            model, param = single_output_model(method, X_train_var, Y_train_var, params=params, best_params=best_params[var])
            best_parameters[var] = param

            y_pred = model.predict(X_valid_var)
        else:
            y_pred = Y_valid_var
            print("It must be one of the following methods: ridge, lasso, ElasticNet, SVM_linear and SVM_RBF")

        # calculated performance

        mse = mean_squared_error(Y_valid_var, y_pred)
        corr, p = stats.spearmanr(y_pred, Y_valid_var)
        pearson_corr, p = stats.pearsonr(y_pred, Y_valid_var)
        rmse = math.sqrt(mse)
        MAE = mean_absolute_error(Y_valid_var, y_pred)
        r2 = r2_score(Y_valid_var, y_pred)
        performances[var] = {"MSE": mse,
                             "MAE": MAE,
                             "RMSE": rmse,
                             "R2": r2,
                             "corr": corr,
                             "pearson_corr": pearson_corr}
        gc.collect()  
    return performances, best_parameters

