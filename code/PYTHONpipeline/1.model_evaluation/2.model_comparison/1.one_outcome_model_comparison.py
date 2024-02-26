import os
import pickle
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from configs import *
import seaborn as sns


def reformat_models_performance(method, path_to_performance, measure):
    performance = pickle.load(open(path_to_performance + "performances.p", "rb"))
    if method == "MTL":
        reformat_performance = {}
        for var in phenotypes:
            reformat_performance[var] = []
            for fold_idx in range(5):
                keys = list(performance[fold_idx][var].keys())
                if measure not in keys:
                    fname = list(performance[fold_idx][var].keys())[0]
                    reformat_performance[var].append(performance[fold_idx][var][fname][measure][0])
                else:
                    reformat_performance[var].append(performance[fold_idx][var][measure])
    else:
        reformat_performance = {}
        for var in phenotypes:
            reformat_performance[var] = []
            for fold_idx in range(5):
                reformat_performance[var].append(performance[fold_idx][var][measure])

    reformat_performance_mean = []
    reformat_performance_sd = []
    for var in phenotypes:
        reformat_performance_mean.append(np.mean(reformat_performance[var]))
        reformat_performance_sd.append(np.std(reformat_performance[var], ddof=1))

    return reformat_performance, reformat_performance_mean, reformat_performance_sd


if not os.path.isdir(model_evaluation_res_path + 'model_comparison/one_outcome_model/'):
    os.makedirs(model_evaluation_res_path + 'model_comparison/one_outcome_model/')


# -----------------------------------------Compare the performance of multi-task learning models and ML models-------------------------------------------------#
for measure in ["RMSE", "pearson_corr"]:
    #MTL
    path_to_performance = MTL_model_path + "MTL_Hot/CV/Random_split_CV/rounds_best_hyperparms_outputs/"
    MTL_performance, MTL_mean, MTL_sd = reformat_models_performance(method="MTL",
                                                                    path_to_performance=path_to_performance,
                                                                    measure=measure)
    path_to_performance = MTL_model_path + "MTL_Hot/CV/Add_noise_CV/rounds_best_hyperparms_outputs/"
    noise_MTL_performance, noise_MTL_mean, noise_MTL_sd = reformat_models_performance(method="MTL",
                                                                                      path_to_performance=path_to_performance,
                                                                                      measure=measure)

    # linear models
    # ridge
    method = "ridge"
    path_to_performance = Linear_models_res_path + "/%s/Random_split_" % method
    ridge_performance, ridge_mean, ridge_sd = reformat_models_performance(method=method,
                                                                          path_to_performance=path_to_performance,
                                                                          measure=measure)
    path_to_performance = Linear_models_res_path + "/%s/Noise_" % method
    noise_ridge_performance, noise_ridge_mean, noise_ridge_sd = \
        reformat_models_performance(method=method,
                                    path_to_performance=path_to_performance,
                                    measure=measure)

    # lasso
    method = "lasso"
    path_to_performance = Linear_models_res_path + "/%s/Random_split_" % method
    lasso_performance, lasso_mean, lasso_sd = reformat_models_performance(method=method,
                                                                          path_to_performance=path_to_performance,
                                                                          measure=measure)
    path_to_performance = Linear_models_res_path + "/%s/Noise_" % method
    noise_lasso_performance, noise_lasso_mean, noise_lasso_sd = \
        reformat_models_performance(method=method,
                                    path_to_performance=path_to_performance,
                                    measure=measure)

    # ElasticNet
    method = "ElasticNet"
    path_to_performance = Linear_models_res_path + "/%s/Random_split_" % method
    ElasticNet_performance, ElasticNet_mean, ElasticNet_sd = reformat_models_performance(method=method,
                                                                                         path_to_performance=path_to_performance,
                                                                                         measure=measure)
    path_to_performance = Linear_models_res_path + "/%s/Noise_" % method
    noise_ElasticNet_performance, noise_ElasticNet_mean, noise_ElasticNet_sd = \
        reformat_models_performance(method=method,
                                    path_to_performance=path_to_performance,
                                    measure=measure)

    # non-linear models
    method = "SVM_linear"
    path_to_performance = non_Linear_models_res_path + "/%s/Random_split_" % method
    SVM_linear_performance, SVM_linear_mean, SVM_linear_sd = reformat_models_performance(method=method,
                                                                                         path_to_performance=path_to_performance,
                                                                                         measure=measure)
    path_to_performance = non_Linear_models_res_path + "/%s/Noise_" % method
    noise_SVM_linear_performance, noise_SVM_linear_mean, noise_SVM_linear_sd = \
        reformat_models_performance(method=method,
                                    path_to_performance=path_to_performance,
                                    measure=measure)

    method = "SVM_RBF"
    path_to_performance = non_Linear_models_res_path + "/%s/Random_split_" % method
    SVM_RBF_performance, SVM_RBF_mean, SVM_RBF_sd = reformat_models_performance(method=method,
                                                                                path_to_performance=path_to_performance,
                                                                                measure=measure)
    path_to_performance = non_Linear_models_res_path + "/%s/Noise_" % method
    noise_SVM_RBF_performance, noise_SVM_RBF_mean, noise_SVM_RBF_sd = \
        reformat_models_performance(method=method,
                                    path_to_performance=path_to_performance,
                                    measure=measure)

    method = "RandomForest"
    path_to_performance = non_Linear_models_res_path + "/%s/Random_split_" % method
    RandomForest_performance, RandomForest_mean, RandomForest_sd = reformat_models_performance(method=method,
                                                                                path_to_performance=path_to_performance,
                                                                                measure=measure)
    path_to_performance = non_Linear_models_res_path + "/%s/Noise_" % method
    noise_RandomForest_performance, noise_RandomForest_mean, noise_RandomForest_sd = \
        reformat_models_performance(method=method,
                                    path_to_performance=path_to_performance,
                                    measure=measure)



    with PdfPages(model_evaluation_res_path + 'model_comparison/one_outcome_model/' + measure + '_barplot11.pdf') as pdf:
       
        for i in range(4):
            var = phenotypes[i]
            labels = ["MTL", "ridge", "lasso", "ElasticNet", "SVM_linear", "SVM_RBF", "RandomForest"]
            plt.figure(figsize=(7, 5))
            axs = plt.gca()
            axs.spines['right'].set_color('none')
            axs.spines['top'].set_color('none')
            per_mean = [MTL_mean[i], ridge_mean[i], lasso_mean[i], ElasticNet_mean[i],
                        SVM_linear_mean[i], SVM_RBF_mean[i], RandomForest_mean[i]]
            plt.bar(labels, per_mean, alpha=0.7)
     
            scatter_x_labels = np.repeat(labels, 5, axis=0)
            scatter_y = MTL_performance[var] + ridge_performance[var] + lasso_performance[var] + ElasticNet_performance[var] + \
                SVM_linear_performance[var] + SVM_RBF_performance[var] + RandomForest_performance[var]
            plt.scatter(scatter_x_labels, scatter_y, color='black', s=10)
            plt.xticks(size=16, rotation=30)  
            plt.yticks(size=16)
            plt.ylabel(measure, fontsize=20)
            if measure == "pearson_corr":
                plt.ylim(0.78, 1)  
            plt.title(var, fontsize=20)
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()



  
    wilcox_p = {}
    for var in phenotypes:
        print(var)
        t1, MTL_ridge_wilcox_p = stats.ranksums(MTL_performance[var], ridge_performance[var])
        t2, MTL_lasso_wilcox_p = stats.ranksums(MTL_performance[var], lasso_performance[var])
        t3, MTL_ElasticNet_wilcox_p = stats.ranksums(MTL_performance[var], ElasticNet_performance[var])
        t4, SVM_linear_wilcox_p = stats.ranksums(MTL_performance[var], SVM_linear_performance[var])
        t5, SVM_RBF_wilcox_p = stats.ranksums(MTL_performance[var], SVM_RBF_performance[var])
        t6, RandomForest_wilcox_p = stats.ranksums(MTL_performance[var], RandomForest_performance[var])
        print("MTL_ridge_wilcox_p: %f" % MTL_ridge_wilcox_p)
        print("MTL_lasso_wilcox_p: %f" % MTL_lasso_wilcox_p)
        print("MTL_ElasticNet_wilcox_p: %f" % MTL_ElasticNet_wilcox_p)
        print("SVM_linear_wilcox_p: %f" % SVM_linear_wilcox_p)
        print("SVM_RBF_wilcox_p: %f" % SVM_RBF_wilcox_p)
        print("RandomForest_wilcox_p: %f" % RandomForest_wilcox_p)
        wilcox_p[var] = [MTL_ridge_wilcox_p, MTL_lasso_wilcox_p, MTL_ElasticNet_wilcox_p, SVM_linear_wilcox_p,
                         SVM_RBF_wilcox_p, RandomForest_wilcox_p]
   
    pickle.dump(wilcox_p,
                open(model_evaluation_res_path + 'model_comparison/one_outcome_model/Random_split_' + measure + "_wilcox_p.pkl", 'wb'))


 
    wilcox_p = {}
    for var in phenotypes:
        print(var)
        t1, noise_MTL_wilcox_p = stats.ranksums(noise_MTL_performance[var], MTL_performance[var])
        t2, noise_ridge_wilcox_p = stats.ranksums(noise_ridge_performance[var], ridge_performance[var])
        t3, noise_lasso_wilcox_p = stats.ranksums(noise_lasso_performance[var], lasso_performance[var])
        t4, noise_ElasticNet_wilcox_p = stats.ranksums(noise_ElasticNet_performance[var], ElasticNet_performance[var])
        t5, noise_SVM_linear_wilcox_p = stats.ranksums(noise_SVM_linear_performance[var], SVM_linear_performance[var])
        t6, noise_SVM_RBF_wilcox_p = stats.ranksums(noise_SVM_RBF_performance[var], SVM_RBF_performance[var])
        t7, noise_RandomForest_wilcox_p = stats.ranksums(noise_RandomForest_performance[var], RandomForest_performance[var])
        print("noise_MTL_wilcox_p: %f" % noise_MTL_wilcox_p)
        print("noise_ridge_wilcox_p: %f" % noise_ridge_wilcox_p)
        print("noise_lasso_wilcox_p: %f" % noise_lasso_wilcox_p)
        print("noise_ElasticNet_wilcox_p: %f" % noise_ElasticNet_wilcox_p)
        print("noise_SVM_linear_wilcox_p: %f" % noise_SVM_linear_wilcox_p)
        print("noise_SVM_RBF_wilcox_p: %f" % noise_SVM_RBF_wilcox_p)
        print("noise_RandomForest_wilcox_p: %f" % noise_RandomForest_wilcox_p)
        wilcox_p[var] = [noise_MTL_wilcox_p, noise_ridge_wilcox_p, noise_lasso_wilcox_p, noise_ElasticNet_wilcox_p,
                         noise_SVM_linear_wilcox_p, noise_SVM_RBF_wilcox_p, noise_RandomForest_wilcox_p]


    pickle.dump(wilcox_p,
                open(
                    model_evaluation_res_path + 'model_comparison/one_outcome_model/Noise_' + measure + "_wilcox_p.pkl",
                    'wb'))


    
    with PdfPages(model_evaluation_res_path + 'model_comparison/one_outcome_model/' + measure + '_boxplot11.pdf') as pdf:
        
        for i in range(4):
            var = phenotypes[i]
            labels = ["MTL", "ridge", "lasso", "ElasticNet", "SVM_linear", "SVM_BRF", "RandomForest"]
            plt.figure(figsize=(10, 5))
            x_labels = np.repeat(labels, 5, axis=0)
            y = MTL_performance[var] + ridge_performance[var] + lasso_performance[var] + ElasticNet_performance[var] + \
                SVM_linear_performance[var] + SVM_RBF_performance[var] + RandomForest_performance[var]
            df = pd.DataFrame({'Method': x_labels, 'Performance': y})
            sns.boxplot(x='Method', y='Performance', data=df, width=0.6)
            sns.stripplot(
                data=df, x="Method", y="Performance", hue=None,
                alpha=.7, zorder=1, legend=False, jitter=0.1, dodge=True
            )  

            plt.ylabel(measure)  
            if measure == "pearson_corr":
                plt.ylim(0.78, 1)
            plt.title(var)
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

 
    with PdfPages(model_evaluation_res_path + 'model_comparison/one_outcome_model/noise_' + measure + '_boxplot11.pdf') as pdf:
       
        for i in range(4):
            var = phenotypes[i]
            labels = ["MTL", "ridge", "lasso", "ElasticNet", "SVM_linear", "SVM_BRF", "RandomForest"]
            plt.figure(figsize=(10, 5))
            x_labels = np.repeat(labels, 10, axis=0)
            y = MTL_performance[var] + noise_MTL_performance[var] +\
                ridge_performance[var] + noise_ridge_performance[var] +\
                lasso_performance[var] + noise_lasso_performance[var] + \
                ElasticNet_performance[var] + noise_ElasticNet_performance[var] +\
                SVM_linear_performance[var] + noise_SVM_linear_performance[var] +\
                SVM_RBF_performance[var] + noise_SVM_RBF_performance[var] + \
                RandomForest_performance[var] + noise_RandomForest_performance[var]
            noise = np.tile(np.concatenate((np.repeat("No noise", 5), np.repeat("Add noise", 5)), axis=0), 7)
            df = pd.DataFrame({'Method': x_labels, 'Performance': y, 'Noise': noise})

            sns.boxplot(x='Method', y='Performance', hue='Noise', data=df, width=0.6)
            sns.stripplot(
                data=df, x="Method", y="Performance", hue='Noise',
                alpha=.7, zorder=1, legend=False, jitter=0.1, dodge=True
            )  
            plt.ylabel(measure)  
            plt.title(var)
            if measure == "pearson_corr":
                plt.ylim(0.78, 1)  
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()
