import os
os.environ["OMP_NUM_THREADS"] = '1'  
from DEC import *
import joblib
from load_data import *
import pandas as pd

data_exp_tcga, exp_gene_tcga, exp_sample_tcga = load_feature_data(path_to_exp_info)
pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)

# ------------------------------------重新标化每套验证集的表达谱----------------------------------------------------#
exclude_data_list = ["ExpO_project", "Braun.NatMed.2021", "Mariathasan.Nature.2018", "Cirenajwis.Oncotarget.2015",
                     "Cui.NPJGenomMed.2016_Riaz17", "Cui.NPJGenomMed.2016_Gide19", "Cui.NPJGenomMed.2016_Hugo16",
                     "Cui.NPJGenomMed.2016_Liu19", "Cui.NPJGenomMed.2016_PUCH", "Cui.NPJGenomMed.2016_VanAllen","Sanjeev.Nature.2018"]

for dset in exclude_data_list:
    print(dset)
    validate_data_path = path_to_validation_data + "%s/processed_data/" % dset
    per_validate_res_path = validate_res_path + "%s/" % dset

    path_to_exp_info = validate_data_path + "log_exp_data.csv"
    data_exp_val, exp_gene_names, exp_sample_names = load_feature_data(path_to_exp_info)
    print(exp_gene_tcga.shape)
    print(exp_gene_names.shape)

    valid_exp_genes_idx = np.where(exp_gene_tcga == exp_gene_names[:, None])[1] 
    X_tcga = data_exp_tcga[:, valid_exp_genes_idx]
    normalized_X_tcga, Y = load_final_data(X_tcga, pheno_score, shuffle=False)
    # normalized the expression profile of the validation set
    scaled_val_exp = normalized_val_data(val_feature_data=data_exp_val, train_feature_data=normalized_X_tcga)
    data1 = pd.DataFrame(columns=exp_gene_names, data=scaled_val_exp)
    data1.to_csv(path_to_validation_data + '%s/scaled_val_exp.csv' % dset)  


# 5 datasets of SKCM were imported and combined
data_exp, skcm_intersect_genes, exp_sample = load_feature_data(path_to_validation_data + "combind_skcm_datasets/processed_data/log_exp_data.csv")
dset = "Cui.NPJGenomMed.2016_Riaz17"
data1 = pd.read_csv(path_to_validation_data + '%s/scaled_val_exp.csv' % dset, sep=",", index_col=0)
intersect_genes_idx = np.where(data1.columns.values == skcm_intersect_genes[:, None])[1]  
scaled_val_exp = data1.values[0::, 0::]
combind_skcm_datasets = scaled_val_exp[:, intersect_genes_idx]
for dset in ["Cui.NPJGenomMed.2016_Gide19", "Cui.NPJGenomMed.2016_Hugo16",
             "Cui.NPJGenomMed.2016_Liu19", "Cui.NPJGenomMed.2016_PUCH"]:
    data1 = pd.read_csv(path_to_validation_data + '%s/scaled_val_exp.csv' % dset, sep=",", index_col=0)
    intersect_genes_idx = np.where(data1.columns.values == skcm_intersect_genes[:, None])[1]  
    scaled_val_exp = data1.values[0::, 0::]
    combind_skcm_datasets = np.concatenate((combind_skcm_datasets, scaled_val_exp[:, intersect_genes_idx]), axis=0)

data2 = pd.DataFrame(columns=skcm_intersect_genes, data=combind_skcm_datasets)
data2.to_csv(path_to_validation_data + 'combind_skcm_datasets/scaled_val_exp.csv')


# ----------------------------------clustering--------------------------------------------#
exclude_data_list = ["ExpO_project", "Braun.NatMed.2021", "Mariathasan.Nature.2018", "Cirenajwis.Oncotarget.2015",
                     "Cui.NPJGenomMed.2016_Riaz17", "Cui.NPJGenomMed.2016_Gide19", "Cui.NPJGenomMed.2016_Hugo16",
                     "Cui.NPJGenomMed.2016_Liu19", "Cui.NPJGenomMed.2016_PUCH", "Cui.NPJGenomMed.2016_VanAllen",
                     "combind_skcm_datasets","Sanjeev.Nature.2018"]
for dset in exclude_data_list:
    print(dset)

    data1 = pd.read_csv(path_to_validation_data + '%s/scaled_val_exp.csv' % dset, sep=",", index_col=0)
    scaled_val_exp = data1.values[0::, 0::]
    print(scaled_val_exp.shape)
    final_model_path = validate_res_path + "%s/rescaled_final_model/" % dset
    save_dir = path_to_embedding_cluster + "%s/" % dset
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)
    print("load data")

    exp_pca = joblib.load(final_model_path + "exp_pca.m")
    X_val_transformed = exp_pca.transform(scaled_val_exp)
    print("PCA transformed")

    n_clusters = 2
    rep = 0

    path_to_model = final_model_path + str(rep) + "/prediction_model.h5"
    HIDDEN_LAYER = 2
    mod_to_layer = get_model_layers(path_to_model, num_layers=HIDDEN_LAYER)

    # ------------------------------ Construct the model after adding the clustering layer-----------------------------#
    model = DEC(mod_to_layer, n_clusters,
                X_train=X_val_transformed,
                maxiterV=1000,
                update_interval=100,
                tol=1/X_val_transformed.shape[0],
                batch_size=10)
    model.save_weights(save_dir + 'DEC_model_final_%s.h5' % str(n_clusters))  

    # --------------------------------- Predict the clustering results of the sample in the validation set ----------------------------------#
    q = model.predict(X_val_transformed, verbose=0)
    p = target_distribution(q)  # update the auxiliary target distribution p

    y_pred = q.argmax(1)
    data1 = pd.DataFrame(y_pred)
    data1.to_csv(save_dir + 'y_pred_%s.csv' % str(n_clusters))
