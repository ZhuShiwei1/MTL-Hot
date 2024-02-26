# 使用和验证数据具有相同特征的TCGA数据重新训练模型
import joblib
from Train_CV_model import *

# 导入最终模型的最优超参数组合
path_to_final_chosen_models = MTL_model_path + "MTL_Hot/CV/Random_split_CV/final_models_chosen/"
final_best_fname = pickle.load(open(path_to_final_chosen_models + "/final.p", "rb"))

# 将字符串转化为字典
hyperparams_values = final_best_fname.split("_")
MDNN_hyperparams = {"epochs": 200,
                    "inner_activation": "relu",
                    "hidden_sizes_shared": eval(hyperparams_values[0]),
                    "hidden_sizes_separate": eval(hyperparams_values[1]),
                    "dropout": float(hyperparams_values[2]),
                    "k_reg": float(hyperparams_values[3]),
                    "learning_rate": float(hyperparams_values[4]),
                    "batch_size": int(hyperparams_values[5])}
print(MDNN_hyperparams)
# 导入tcga的数据
data_exp_tcga, exp_gene_tcga, exp_sample_names = load_feature_data(path_to_exp_info)
pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)

"""
exclude_data_list = ["ExpO_project", "Braun.NatMed.2021", "Mariathasan.Nature.2018", "Cirenajwis.Oncotarget.2015",
                     "Cui.NPJGenomMed.2016_Riaz17", "Cui.NPJGenomMed.2016_Gide19", "Cui.NPJGenomMed.2016_Hugo16",
                     "Cui.NPJGenomMed.2016_Liu19", "Cui.NPJGenomMed.2016_PUCH", "Cui.NPJGenomMed.2016_VanAllen",
                     "combind_skcm_datasets"]
"""
exclude_data_list = ["Sanjeev.Nature.2018"]
for dset in exclude_data_list:
    print(dset)
    validate_data_path = path_to_validation_data + "%s/processed_data/" % dset
    # retrain后的模型存储路径
    final_model_path = validate_res_path + "%s/rescaled_final_model/" % dset
    if not os.path.isdir(final_model_path):
        os.makedirs(final_model_path)
    # ---------------------数据准备-----------------------------#
    # 导入验证数据集
    data_exp_valid, exp_gene_valid, exp_sample_valid = load_feature_data(validate_data_path + "log_exp_data.csv")
    print(data_exp_valid.shape)
    # 只保留验证数据集中使用的特征
    valid_exp_genes_idx = np.where(exp_gene_tcga == exp_gene_valid[:, None])[1]  # 返回验证集基因在TCGA基因列表中的位置（按验证基因集的顺序返回idx）
    X = data_exp_tcga[:, valid_exp_genes_idx]
    print(X.shape)
    Y_labels = pheno_score

    # 打乱全部样本的标签
    X_shuffle, Y_shuffle, shuffle_idx = load_final_data(X, Y_labels)

    print("load data")
    # 使用PCA对数据进行降维，并使用降维后的数据作为多任务学习模型的input
    exp_pca = PCA(n_components=300)
    exp_pca.fit(X_shuffle)  # 拟合降维模型
    X_train_transformed = exp_pca.transform(X_shuffle)
    # 存储pca模型
    joblib.dump(exp_pca, final_model_path + "exp_pca.m")
    # 存储打乱后的样本标签顺序
    np.savetxt(final_model_path + "shuffle_idx.txt", shuffle_idx)
    Y_train = Y_shuffle
    print("PCA transformed")

    # ----------------------重复训练模型100次--------------------------------#
    # 重复训练模型100次，并存储每次的模型结果
    input_size = X_train_transformed.shape[1]
    if dset in ["ExpO_project"]:
        run_times = 100
    else:
        run_times = 1

    for rep in range(run_times):
        print(rep)
        model_path = final_model_path + str(rep) + "/"
        if not os.path.isdir(model_path):
            os.makedirs(model_path)

        prediction_model = get_prediction_model(MDNN_hyperparams, input_size)
        trainable_model = get_trainable_model(prediction_model, input_size)
        opt = adam_v2.Adam(learning_rate=MDNN_hyperparams["learning_rate"])  # 学习率不宜设置过大
        trainable_model.compile(optimizer=opt, loss=None)
        trainable_model.fit([X_train_transformed,
                             Y_train[phenotypes[0]], Y_train[phenotypes[1]],
                             Y_train[phenotypes[2]], Y_train[phenotypes[3]]],
                            epochs=MDNN_hyperparams["epochs"],
                            batch_size=MDNN_hyperparams["batch_size"],
                            shuffle=True,  # 在每个epoch前重新打乱样本顺序
                            verbose=0  # 不输出训练过程
                            )
        loss_weight = [np.exp(K.get_value(log_var[0])) ** 0.5 for log_var in trainable_model.layers[-1].log_vars]
        print(loss_weight)  # 每个任务的权重
        # -------------------------------存储模型以及模型参数-------------------------------------#
        # 存储模型
        prediction_model.compile(optimizer=adam_v2.Adam(), loss=None)  # 这个用不到，但是得写，防止出警告
        prediction_model.save(model_path + "prediction_model.h5")

        # 存储全部的模型参数
        weight_list = []
        for layer in prediction_model.layers:
            weight_list.append(layer.get_weights())
        pickle.dump(weight_list, open(model_path + "prediction_model_weight.p", "wb"))
        K.clear_session()
        gc.collect()
