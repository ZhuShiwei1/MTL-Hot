from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.decomposition import PCA
from load_data import *
from configs import *
import gc
import json

data_exp_tcga, exp_gene_names, exp_sample_names = load_feature_data(path_to_exp_info)

pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)
X_train, Y = load_final_data(X=data_exp_tcga, Y_labels=pheno_score, shuffle=False)

pca = PCA()
X_train_pca = pca.fit_transform(X_train)
exp_var_pca = pca.explained_variance_ratio_
cum_sum_eigenvalues = np.cumsum(exp_var_pca)
print(len(cum_sum_eigenvalues))  #

print(cum_sum_eigenvalues[50])  # 0.669
print(cum_sum_eigenvalues[100])  # 0.728
print(cum_sum_eigenvalues[200])  # 0.779
print(cum_sum_eigenvalues[300])  # 0.80726
print(cum_sum_eigenvalues[400])  # 0.8262
print(cum_sum_eigenvalues[500])  # 0.84090


with PdfPages(root_path+'result/data_prepare/PCA_explained_variance.pdf') as pdf:
    plt.figure()
    plt.bar(range(0, len(exp_var_pca)), exp_var_pca, alpha=0.5, align='center', label='Individual explained variance')
    plt.step(range(0, len(cum_sum_eigenvalues)), cum_sum_eigenvalues, where='mid', label='Cumulative explained variance')
    plt.ylabel('Explained variance ratio')
    plt.xlabel('Principal component index')
    plt.legend(loc='best')
    plt.axvline(300)
    plt.tight_layout()
    pdf.savefig()
    plt.close()


gc.collect()
