# 评估最后一个共享层与各个表型的关联
import os
import gc
import joblib
from keras.models import Model, load_model
from keras import backend as K
from keras.optimizers import adam_v2
from load_data import *


def get_model_layers(model_file, num_layers):
  
    model = load_model(model_file)
    # define new model that cuts off the last several layers
    newmodel = Model(inputs=model.input, outputs=model.layers[num_layers - 1].output)
    # agian, need to specify these parameters, but they aren't used since we don't retrain the model
    opt = adam_v2.Adam()
    newmodel.compile(optimizer=opt, loss=None)
    return newmodel


# ---------------------------- Calculate the outputs of the last shared layer -------------------------------------#

data_exp_tcga, exp_gene_names, exp_sample_names = load_feature_data(path_to_exp_info)
pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)

X = data_exp_tcga
Y_labels = pheno_score
normalized_X, Y = load_final_data(X, Y_labels, shuffle=False)
print("load data")

exp_pca = joblib.load(MTL_model_path+"MTL_Hot/final_model/exp_pca.m")
X_train_transformed = exp_pca.transform(normalized_X)

path_to_models = MTL_model_path + "MTL_Hot/final_model/"
final_rep_embeddings_path = MTL_model_path + "MTL_Hot/final_embeddings/"
if not os.path.isdir(final_rep_embeddings_path):
    os.makedirs(final_rep_embeddings_path)

for rep in range(100):
    print(rep)
    path_to_model = path_to_models + "%s/prediction_model.h5" % rep
    HIDDEN_LAYER = 2
    mod_to_layer = get_model_layers(path_to_model, num_layers=HIDDEN_LAYER)
    X_transformed = mod_to_layer.predict(X_train_transformed)
    np.savetxt("%s/%s.txt" % (final_rep_embeddings_path, rep), X_transformed)

    K.clear_session()
    gc.collect()
