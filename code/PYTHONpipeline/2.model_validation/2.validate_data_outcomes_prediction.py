
# 对额外的验证集进行预测，并取100个模型预测结果的均值作为最终的模型预测结果
import joblib
import scipy
from keras.models import load_model
from matplotlib import pyplot as plt
from scipy.stats import stats
from sklearn import linear_model
import math
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from load_data import *
from matplotlib.backends.backend_pdf import PdfPages


# 导入tcga的数据
data_exp_tcga, exp_gene_tcga, exp_sample_names = load_feature_data(path_to_exp_info)
pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)

exclude_data_list = ["ExpO_project", "Braun.NatMed.2021", "Mariathasan.Nature.2018"]
#exclude_data_list = ["Sanjeev.Nature.2018"]
for dset in exclude_data_list:
    validate_data_path = path_to_validation_data + "%s/processed_data/" % dset
    per_validate_res_path = validate_res_path + "%s/" % dset
    # retrain后的模型存储路径
    final_model_path = validate_res_path + "%s/rescaled_final_model/" % dset

    # 导入验证集的表型分数
    path_to_val_phenotypic_info = validate_data_path + "cancer.phenotype.csv"
    val_pheno_score, val_pheno_labels, val_sample_names = load_phenotype_data(path_to_val_phenotypic_info)

    # 导入验证数据集的特征
    path_to_exp_info = validate_data_path + "log_exp_data.csv"
    data_exp_val, exp_gene_names, exp_sample_names = load_feature_data(path_to_exp_info)
    # 只保留验证数据集中使用的特征
    valid_exp_genes_idx = np.where(exp_gene_tcga == exp_gene_names[:, None])[1]  # 返回验证集基因在TCGA基因列表中的位置（按验证基因集的顺序返回idx）
    X_tcga = data_exp_tcga[:, valid_exp_genes_idx]
    normalized_X_tcga, Y = load_final_data(X_tcga, pheno_score, shuffle=False)
    # 根据训练模型的所有TCGA样本中基因的均值方差重新标化验证集的表达谱
    scaled_val_exp = normalized_val_data(val_feature_data=data_exp_val, train_feature_data=normalized_X_tcga)
    print("load data")

    # 根据之前的pca模型直接对数据降维
    exp_pca = joblib.load(final_model_path + "exp_pca.m")
    X_val_transformed = exp_pca.transform(scaled_val_exp)
    X_tcga_transformed = exp_pca.transform(normalized_X_tcga)
    print("PCA transformed")

    # 使用模型直接预测结果
    all_pred_vals = []
    all_pred_tcga = []
    reps = 1
    for rep in range(reps):
        model_path = final_model_path + str(rep) + "/"
        # 导入模型
        prediction_model = load_model(model_path + "prediction_model.h5")

        # 验证集的预测值
        y_predicts = prediction_model.predict(X_val_transformed)
        per_pred_res = []

        with PdfPages(model_path + 'corr_plot.pdf') as pdf:
            for i in range(4):
                var = phenotypes[i]
                y_true = val_pheno_score.iloc[:, i].astype(float).values
                y_pred = y_predicts[i]
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
                lrModel.fit(y_pred, y_true)
                deter_coef = lrModel.score(y_pred, y_true)
                plt.plot(y_pred, lrModel.predict(y_pred), color='red')  # 拟合回归直线
                pdf.savefig()  # saves the current figure into a pdf page
                plt.close()
                #np.savetxt(model_path + var + "_predicts.txt", y_pred)
                per_pred_res.append(np.hstack(y_pred))

        pred_df = pd.DataFrame(np.array(per_pred_res, dtype=float)).T
        pred_df.columns = phenotypes
        pred_df.to_csv(final_model_path + "val_preds.csv")
        all_pred_vals.append(pred_df.values)

        # tcga的预测值
        tcga_y_predicts = prediction_model.predict(X_tcga_transformed)
        per_pred_res = []
        for i in range(4):
            var = phenotypes[i]
            y_pred = tcga_y_predicts[i]
            #np.savetxt(model_path + var + "_tcga_predicts.txt", y_pred)
            per_pred_res.append(np.hstack(y_pred))
        tcga_pred_df = pd.DataFrame(np.array(per_pred_res, dtype=float)).T
        tcga_pred_df.columns = phenotypes
        #tcga_pred_df.to_csv(final_model_path + "tcga_preds.csv")
        all_pred_tcga.append(tcga_pred_df.values)

        # 计算性能
        performances = {}
        for i in range(4):
            var = phenotypes[i]
            y_true = val_pheno_score.iloc[:, i].astype(float).values
            y_pred = y_predicts[i]
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
        pickle.dump(performances, open(model_path + "performances.p", "wb"))
