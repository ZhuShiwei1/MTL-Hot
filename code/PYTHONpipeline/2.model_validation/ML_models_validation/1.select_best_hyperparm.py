import os
import joblib
import scipy
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn import linear_model
from load_data import *
from single_output_model_performance import *

# 导入数据
data_exp_tcga, exp_gene_tcga, exp_sample_names = load_feature_data(path_to_exp_info)
pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)

exclude_data_list = ["ExpO_project", "Braun.NatMed.2021", "Mariathasan.Nature.2018"]

for dset in exclude_data_list:
    print(dset)
    validate_data_path = path_to_validation_data + "%s/processed_data/" % dset

    # retrain后的模型存储路径
    final_model_path = validate_res_path + "%s/rescaled_final_ML_model/" % dset
    if not os.path.isdir(final_model_path):
        os.makedirs(final_model_path)

    # ---------------------数据准备-----------------------------#
    # 导入验证集的表型分数
    path_to_val_phenotypic_info = validate_data_path + "cancer.phenotype.csv"
    val_pheno_score, val_pheno_labels, val_sample_names = load_phenotype_data(path_to_val_phenotypic_info)

    # 导入验证数据集的特征
    path_to_exp_info = validate_data_path + "log_exp_data.csv"
    data_exp_val, exp_gene_names, exp_sample_names = load_feature_data(path_to_exp_info)
    # 只保留验证数据集中使用的特征
    valid_exp_genes_idx = np.where(exp_gene_tcga == exp_gene_names[:, None])[1]  # 返回验证集基因在TCGA基因列表中的位置（按验证基因集的顺序返回idx）
    X_tcga = data_exp_tcga[:, valid_exp_genes_idx]
    normalized_X_tcga, Y_train = load_final_data(X_tcga, pheno_score, shuffle=False)
    # 根据训练模型的所有TCGA样本中基因的均值方差重新标化验证集的表达谱
    scaled_val_exp = normalized_val_data(val_feature_data=data_exp_val, train_feature_data=normalized_X_tcga)
    print("load data")

    # 使用PCA对数据进行降维，并使用降维后的数据作为多任务学习模型的input
    exp_pca = PCA(n_components=300)
    exp_pca.fit(normalized_X_tcga)  # 拟合降维模型
    X_val_transformed = exp_pca.transform(scaled_val_exp)
    X_tcga_transformed = exp_pca.transform(normalized_X_tcga)
    # 存储pca模型
    joblib.dump(exp_pca, final_model_path + "exp_pca.m")

    print("PCA transformed")

    #for method in ["ridge", "lasso", "ElasticNet", "SVM_linear", "SVM_RBF", "RandomForest"]:
    for method in ["RandomForest"]:
        # retrain后的模型结果存储路径
        final_model_res_path = validate_res_path + "%s/rescaled_final_ML_model/%s/" % (dset, method)
        if not os.path.isdir(final_model_res_path):
            os.makedirs(final_model_res_path)
        # 针对每个变量训练ML模型，并在验证集中进行验证
        # 针对全部的TCGA样本识别最优超参数组合
        performances = {}
        #best_parameters = {}
        per_pred_res = []
        for i in range(4):
            var = phenotypes[i]
            # 真实值
            y_true = val_pheno_score.iloc[:, i].astype(float).values
            # 针对特定的表型
            Y_train_var = Y_train[var]
            print(var)
            # ridge
            if method in ["ridge", "lasso", "ElasticNet"]:
                best_parameters = pickle.load(open(Linear_models_res_path + "/%s/" % method + "best_parameters.p", "rb"))
                model, parameters = single_output_model(method, X_tcga_transformed, Y_train_var,
                                                        best_params=best_parameters[var])
            elif method == "SVM_linear":
                model, parameters = single_output_model(method, X_tcga_transformed, Y_train_var)

            elif method in ["SVM_RBF", "RandomForest"]:
                best_parameters = pickle.load(
                    open(non_Linear_models_res_path + "/%s/" % method + "best_parameters.p", "rb"))
                model, parameters = single_output_model(method, X_tcga_transformed, Y_train_var,
                                                        best_params=best_parameters[var])
            # 模型预测
            y_pred = model.predict(X_val_transformed)
            #best_parameters[var] = parameters
            per_pred_res.append(np.hstack(y_pred))

            # 计算性能
            mse = mean_squared_error(y_true, y_pred)
            corr, p = scipy.stats.spearmanr(y_pred, y_true)
            pearson_corr, p1 = scipy.stats.pearsonr(y_pred.reshape(y_true.shape), y_true)
            rmse = math.sqrt(mse)
            MAE = mean_absolute_error(y_true, y_pred)
            r2 = r2_score(y_true, y_pred)
            performances[var] = {"MSE": mse,
                                 "MAE": MAE,
                                 "RMSE": rmse,
                                 "R2": r2,
                                 "corr": corr,
                                 "pearson_corr": pearson_corr}

            # 绘制相关点图
            with PdfPages(final_model_res_path + '%s_%s_corr_plot.pdf' % (method, var)) as pdf:
                # 计算相关性
                pearson_corr, p2 = scipy.stats.pearsonr(y_pred.reshape(y_true.shape), y_true)
                plt.figure(figsize=(4, 4))
                # 可视化预测值与真实值之间的相关性
                plt.scatter(y_pred, y_true, color="blue")
                plt.title(var + ', pearson_corr=%.3f' % pearson_corr)
                plt.ylabel('Original labels')
                plt.xlabel('Predicted labels')
                plt.grid()
                # 绘制最佳拟合曲线
                lrModel = linear_model.LinearRegression()
                lrModel.fit(y_pred.reshape(-1, 1), y_true)
                deter_coef = lrModel.score(y_pred.reshape(-1, 1), y_true)
                plt.plot(y_pred.reshape(-1, 1), lrModel.predict(y_pred.reshape(-1, 1)), color='red')  # 拟合回归直线
                pdf.savefig()  # saves the current figure into a pdf page
                plt.close()

        # 存储预测结果和最优超参数
        pred_df = pd.DataFrame(np.array(per_pred_res, dtype=float)).T
        pred_df.columns = phenotypes
        pred_df.to_csv(final_model_res_path + "val_preds.csv")

        pickle.dump(performances, open(final_model_res_path + "performances.p", "wb"))
        #pickle.dump(best_parameters, open(final_model_res_path + "best_parameters.p", "wb"))

        gc.collect()  # 回归垃圾
