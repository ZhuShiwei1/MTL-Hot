# 获得验证集在最后一个共享层的预测值
import os
import gc
import joblib
from keras.models import Model, load_model
from keras import backend as K
from keras.optimizers import adam_v2
from load_data import *


def get_model_layers(model_file, num_layers):
    # 导入模型
    model = load_model(model_file)
    # define new model that cuts off the last several layers
    newmodel = Model(inputs=model.input, outputs=model.layers[num_layers - 1].output)
    # agian, need to specify these parameters, but they aren't used since we don't retrain the model
    opt = adam_v2.Adam()
    newmodel.compile(optimizer=opt, loss=None)
    return newmodel


# ---------------------------- 计算最后一个共享层预测值 -------------------------------------#
# 导入tcga的数据
data_exp_tcga, exp_gene_tcga, exp_sample_names = load_feature_data(path_to_exp_info)
pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)
#exclude_data_list = ["ExpO_project", "Braun.NatMed.2021", "Mariathasan.Nature.2018"]
exclude_data_list = ["Sanjeev.Nature.2018"]
for dset in exclude_data_list:
    print(dset)
    validate_data_path = path_to_validation_data + "%s/processed_data/" % dset
    # retrain后的模型存储路径
    final_model_path = validate_res_path + "%s/rescaled_final_model/" % dset

    # 导入验证集数据
    data_exp_valid, exp_gene_names, exp_sample_names = load_feature_data(
        validate_data_path + "log_exp_data.csv")
    # 导入TCGA的训练数据
    valid_exp_genes_idx = np.where(exp_gene_tcga == exp_gene_names[:, None])[1]  # 返回验证集基因在TCGA基因列表中的位置（按验证基因集的顺序返回idx）
    X_tcga = data_exp_tcga[:, valid_exp_genes_idx]
    normalized_X_tcga, Y = load_final_data(X_tcga, pheno_score, shuffle=False)
    # 根据训练模型的所有TCGA样本中基因的均值方差重新标化验证集的表达谱
    scaled_val_exp = normalized_val_data(val_feature_data=data_exp_valid, train_feature_data=normalized_X_tcga)

    print("load data")

    # 根据之前的pca模型直接对数据降维
    exp_pca = joblib.load(final_model_path + "exp_pca.m")
    X_val_transformed = exp_pca.transform(scaled_val_exp)
    X_tcga_transformed = exp_pca.transform(normalized_X_tcga)

    final_rep_embeddings_path = validate_res_path + "%s/final_embeddings/" % dset
    if not os.path.isdir(final_rep_embeddings_path):
        os.makedirs(final_rep_embeddings_path)

    for rep in range(1):
        path_to_model = final_model_path + "%s/prediction_model.h5" % rep
        res_path = final_rep_embeddings_path + "/%s/" % str(rep)
        if not os.path.isdir(res_path):
            os.makedirs(res_path)

        # 感兴趣的神经层
        HIDDEN_LAYER = 2
        mod_to_layer = get_model_layers(path_to_model, num_layers=HIDDEN_LAYER)
        X_val_embedding = mod_to_layer.predict(X_val_transformed)
        np.savetxt("%s/%d/val_rescale.txt" % (final_rep_embeddings_path, rep), X_val_embedding)

        X_tcga_embedding = mod_to_layer.predict(X_tcga_transformed)
        np.savetxt("%s/%d/tcga_rescale.txt" % (final_rep_embeddings_path, rep), X_tcga_embedding)

        K.clear_session()
        gc.collect()
