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


# -----------------------------------------Compare the model performance between neural networks with different layer depths-------------------------------------- #


if not os.path.isdir(model_evaluation_res_path + 'model_comparison/diff_depth_model/'):
    os.makedirs(model_evaluation_res_path + 'model_comparison/diff_depth_model/')


diff_depth = ["Add_shared_layer", "Add_specific_layer", "Reduce_specific_layer", "Split_late"]

for measure in ["MSE", "MAE", "RMSE", "R2", "corr"]:
    print(measure)

    path_to_performance = MTL_model_path + "MTL_Hot/CV/Random_split_CV/rounds_best_hyperparms_outputs/"
    MTL_performance, MTL_mean, MTL_sd = reformat_models_performance(method="MTL",
                                                                    path_to_performance=path_to_performance,
                                                                    measure=measure)
    diff_depth_performance = {}
    diff_depth_mean = {}
    diff_depth_sd = {}
    for depth_type in diff_depth:
        path_to_performance = MTL_model_path + depth_type + "/Random_split_CV/final_models_chosen/"
        diff_depth_performance[depth_type], diff_depth_mean[depth_type], diff_depth_sd[depth_type] = \
            reformat_models_performance(method="MTL",
                                        path_to_performance=path_to_performance,
                                        measure=measure)


    wilcox_p = {}
    diff_depth_wilcox_p = {}
    for var in phenotypes:
        diff_depth_wilcox_p[var] = {}
        for depth_type in diff_depth:
            t, diff_depth_wilcox_p[var][depth_type] = stats.ranksums(MTL_performance[var], diff_depth_performance[depth_type][var])
    print(diff_depth_wilcox_p)

   
    pickle.dump(diff_depth_wilcox_p, open(model_evaluation_res_path + 'model_comparison/diff_depth_model/' + measure + "_wilcox_p.pkl", 'wb'))

    
    with PdfPages(model_evaluation_res_path + 'model_comparison/diff_depth_model/' + measure + '_barplot.pdf') as pdf:
     
        colors = ['#d7191c', '#fdae61', '#abdda4', '#2b83ba']
        for i in range(4):
            var = phenotypes[i]
            labels = ["MTL"]
            labels.extend(diff_depth)
            plt.figure(figsize=(10, 10))
            x = np.arange(4)  # the label locations
            per_mean = [MTL_mean[i]]
            per_mean.extend(diff_depth_mean[depth_type][i] for depth_type in diff_depth)
            plt.bar(labels, per_mean, color=colors[i], alpha=0.7)
         
            scatter_x_labels = np.repeat(labels, 5, axis=0)
            scatter_y = MTL_performance[var]
            for depth_type in diff_depth:
                scatter_y = np.concatenate((scatter_y, np.array(diff_depth_performance[depth_type][var])), axis=0)
            plt.scatter(scatter_x_labels, scatter_y, color='black', s=10)
            plt.xticks(size=16, rotation=20) 
            plt.yticks(size=16)
            plt.ylabel(measure, fontsize=20)  
            plt.title(var, fontsize=20)
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()


    with PdfPages(model_evaluation_res_path + 'model_comparison/diff_depth_model/' + measure + '_boxplot.pdf') as pdf:

        colors = ['#d7191c', '#fdae61', '#abdda4', '#2b83ba']
        for i in range(4):
            var = phenotypes[i]
            labels = ["MTL"]
            labels.extend(diff_depth)
            plt.figure(figsize=(10, 10))
            scatter_x_labels = np.repeat(labels, 5, axis=0)
            y = MTL_performance[var]
            for depth_type in diff_depth:
                y = np.concatenate((y, np.array(diff_depth_performance[depth_type][var])), axis=0)
            df = pd.DataFrame({'MTL_Models': scatter_x_labels, 'Performance': y})
            sns.boxplot(x='MTL_Models', y='Performance', data=df, width=0.6, color='skyblue')
            sns.stripplot(data=df, x="MTL_Models", y="Performance", hue=None,
                          # palette='dark:black',
                          alpha=.5, zorder=1, legend=False, jitter=0.1, dodge=True) 
            plt.ylabel(measure) 
            plt.xticks(size=16, rotation=20)  
            plt.yticks(size=16)
            plt.title(var)
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()
