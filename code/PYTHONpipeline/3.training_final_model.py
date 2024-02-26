from Train_CV_model import *
import joblib

gc.collect()
path_to_final_chosen_models = MTL_model_path+"MTL_Hot/CV/Random_split_CV/final_models_chosen/"
final_best_fname = pickle.load(open(path_to_final_chosen_models + "/final.p", "rb"))

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

data_exp_tcga, exp_gene_names, exp_sample_names = load_feature_data(path_to_exp_info)
pheno_score, pheno_labels, sample_names = load_phenotype_data(path_to_phenotypic_info)

X = data_exp_tcga
Y_labels = pheno_score

X_shuffle, Y_shuffle, shuffle_idx = load_final_data(X, Y_labels)
print("load data")
# PCA
exp_pca = PCA(n_components=300)
exp_pca.fit(X_shuffle)  
X_train_transformed = exp_pca.transform(X_shuffle)
Y_train = Y_shuffle
print("PCA transformed")


if not os.path.isdir(MTL_model_path+"MTL_Hot/final_model/"):
    os.makedirs(MTL_model_path+"MTL_Hot/final_model/")

joblib.dump(exp_pca, MTL_model_path+"MTL_Hot/final_model/exp_pca.m")

np.savetxt(MTL_model_path+"MTL_Hot/final_model/shuffle_idx.txt", shuffle_idx)


# retrain the MTL model for 100 times
input_size = X_train_transformed.shape[1]
for rep in range(100):
    print(rep)
    model_path = MTL_model_path + "MTL_Hot/final_model/" + str(rep) + "/"
    if not os.path.isdir(model_path):
        os.makedirs(model_path)

    prediction_model = get_prediction_model(MDNN_hyperparams, input_size)
    trainable_model = get_trainable_model(prediction_model, input_size)
    opt = adam_v2.Adam(learning_rate=MDNN_hyperparams["learning_rate"])  
    trainable_model.compile(optimizer=opt, loss=None)
    trainable_model.fit([X_train_transformed,
                         Y_train[phenotypes[0]], Y_train[phenotypes[1]],
                         Y_train[phenotypes[2]], Y_train[phenotypes[3]]],
                        epochs=MDNN_hyperparams["epochs"],
                        batch_size=MDNN_hyperparams["batch_size"],
                        shuffle=True,  
                        verbose=0  
                        )
    loss_weight = [np.exp(K.get_value(log_var[0])) ** 0.5 for log_var in trainable_model.layers[-1].log_vars]
    print(loss_weight)  

    # save model
    prediction_model.compile(optimizer=adam_v2.Adam(), loss=None)  
    prediction_model.save(model_path+"prediction_model.h5")

    y_predicts = prediction_model.predict(X_train_transformed)
    for i in range(4):
        var = phenotypes[i]
        y_pred = y_predicts[i]
        np.savetxt(model_path + var + "_predicts.txt", y_pred)

    # save all weights in the model
    weight_list = []
    for layer in prediction_model.layers:
        weight_list.append(layer.get_weights())
    pickle.dump(weight_list, open(model_path + "prediction_model_weight.p", "wb"))
    K.clear_session()
    gc.collect()
