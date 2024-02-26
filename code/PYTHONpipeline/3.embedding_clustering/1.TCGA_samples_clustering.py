import os
os.environ["OMP_NUM_THREADS"] = '1'
from DEC import *
import joblib
from load_data import *
import pandas as pd
import tensorflow as tf
gpus = tf.config.experimental.list_physical_devices('GPU')
for gpu in gpus:
    tf.config.experimental.set_memory_growth(gpu, True)


# ---------------------------------load data---------------------------------#
print("load data:------------")

data_exp_tcga, exp_gene_tcga, exp_sample_names = load_feature_data(path_to_exp_info)
pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)

X = data_exp_tcga
Y_labels = pheno_score

normalized_X, Y = load_final_data(X, Y_labels, shuffle=False)

exp_pca = joblib.load(MTL_model_path+"MTL_Hot/final_model/exp_pca.m")
X_train_transformed = exp_pca.transform(normalized_X)


HIDDEN_LAYER = 2
for rep in range(100):
    print(rep)
    for n_clusters in range(2, 7):
        save_dir = path_to_embedding_cluster + "TCGA/%d/" % rep  
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)

        # ------------------------------ Construct the model after adding the clustering layer -----------------------------#
        path_to_model = MTL_model_path + "MTL_Hot/final_model/%s/prediction_model.h5" % rep
        mod_to_layer = get_model_layers(path_to_model, num_layers=HIDDEN_LAYER)
        model = DEC(mod_to_layer, n_clusters,
                    X_train=X_train_transformed,
                    maxiterV=10000,
                    update_interval=420,
                    tol=0.001,
                    batch_size=200)

        # --------------------------------- Returns the clustering results for each sample ----------------------------------#
        q = model.predict(X_train_transformed, verbose=0)
        p = target_distribution(q)  # update the auxiliary target distribution p

        y_pred = q.argmax(1)
        data1 = pd.DataFrame(y_pred)
        data1.to_csv(save_dir + 'y_pred_%s.csv' % str(n_clusters))
        model.save_weights(save_dir + 'DEC_model_final_%s.h5' % str(n_clusters))  


# --------------------------Validate the model clustering results in a set of validation data------------------------------------#
HIDDEN_LAYER = 2
dset = "ExpO_project"
validate_data_path = path_to_validation_data + "%s/processed_data/" % dset
path_to_exp_info = validate_data_path + "log_exp_data.csv"
data_exp_val, exp_gene_names, val_sample_names = load_feature_data(path_to_exp_info)
valid_exp_genes_idx = np.where(exp_gene_tcga == exp_gene_names[:, None])[1]  
X_tcga = data_exp_tcga[:, valid_exp_genes_idx]
normalized_X_tcga, Y = load_final_data(X_tcga, pheno_score, shuffle=False)

scaled_val_exp = normalized_val_data(val_feature_data=data_exp_val, train_feature_data=normalized_X_tcga)
print("load data")

exp_pca = joblib.load(validate_res_path + "%s/rescaled_final_model/exp_pca.m" % dset)
X_val_transformed = exp_pca.transform(scaled_val_exp)
print("PCA transformed")

for rep in range(100):
    print(rep)
    final_model_path = validate_res_path + "%s/rescaled_final_model/" % dset

    for n_clusters in range(2, 7):
        save_dir = path_to_embedding_cluster + dset + "/%d/" % rep  
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)

        path_to_model = final_model_path + str(rep) + "/prediction_model.h5"

        mod_to_layer = get_model_layers(path_to_model, num_layers=HIDDEN_LAYER)

        # ------------------------------ Construct the model after adding the clustering layer-----------------------------#
        model = DEC(mod_to_layer, n_clusters,
                    X_train=X_val_transformed,
                    maxiterV=10000,
                    update_interval=100,
                    tol=0.001,
                    batch_size=100)
        model.save_weights(save_dir + 'DEC_model_final_%s.h5' % str(n_clusters))  

        # ---------------------------------  predict the clustering results of the sample in the validation set ----------------------------------#
        q = model.predict(X_val_transformed, verbose=0)
        p = target_distribution(q)  # update the auxiliary target distribution p

        y_pred = q.argmax(1)
        data1 = pd.DataFrame(y_pred)
        data1.to_csv(save_dir + 'y_pred_%s.csv' % str(n_clusters))
