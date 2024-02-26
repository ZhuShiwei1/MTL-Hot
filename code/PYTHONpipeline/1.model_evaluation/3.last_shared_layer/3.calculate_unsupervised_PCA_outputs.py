from sklearn.decomposition import PCA
from load_data import *


data_exp_tcga, exp_gene_names, exp_sample_names = load_feature_data(path_to_exp_info)

pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)
X_train, Y = load_final_data(X=data_exp_tcga, Y_labels=pheno_score, shuffle=False)

# PCA
exp_pca = PCA(n_components=100)
X_train_pca = exp_pca.fit_transform(X_train)

print(X_train_pca.shape)

path_to_embedding = MTL_model_path + "MTL_Hot/final_embeddings/"
PCA_pred_df = pd.DataFrame(np.array(X_train_pca, dtype=float))
PCA_pred_df.columns = list(range(1, 100+1))
PCA_pred_df.to_csv(path_to_embedding + "PCA_preds.csv")
