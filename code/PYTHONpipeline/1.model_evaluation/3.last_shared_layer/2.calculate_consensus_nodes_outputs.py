import numpy as np
from sklearn.cluster import KMeans
import pandas as pd
from scipy import stats
from configs import *


def cluster_medoid(zscores_for_cluster):
    distMatrix = np.array([[np.dot(x - y, x - y) for y in zscores_for_cluster] for x in zscores_for_cluster])
    idx = np.argmin(distMatrix.sum(axis=0))
    return idx


def get_last_layer_consensus_nodes(X_transformed_dicts, clusternum=100, rep_num=100, res_path=None):
    # 1)Combine the results of multiple training runs
    print("Combining runs...")
    non_triv_idxes = []
    runs = []
    transformed_vals_prestacked = []
    for rep in np.arange(rep_num):
        non_triv_idx = np.where(np.mean(X_transformed_dicts[rep] != 0, axis=0) > 0)[0]  
        non_triv_idxes.append(non_triv_idx)
        runs.append(np.array([rep] * len(non_triv_idx)))
        transformed_vals_prestacked.append(X_transformed_dicts[rep].T[non_triv_idx])

    transformed_vals_stacked = np.vstack(transformed_vals_prestacked)
    runs = np.hstack(runs)

    # 2)standardized
    transformed_vals_stacked_zscore = stats.zscore(transformed_vals_stacked, axis=1)

    # 3)K-means clustering
    print("Clustering nodes...")
    c = KMeans(n_clusters=clusternum).fit(transformed_vals_stacked_zscore)
    cluster_df = pd.DataFrame(np.vstack([runs, np.hstack(non_triv_idxes), c.labels_]).T,
                              columns=["run", "node_idx", "cluster"])

    # 4)Find the feature closest to the center of the cluster
    selected_medoids = []
    medoid_activations = []
    for cluster_id in range(clusternum):
        rows = np.where(cluster_df["cluster"] == cluster_id)[0]
        selected_idx = cluster_medoid(transformed_vals_stacked_zscore[rows])
        medoid_activations.append(transformed_vals_stacked_zscore[selected_idx])
        selected_medoids.append(rows[selected_idx])

    new_embedding = np.array(medoid_activations).T

    if res_path is not None:
        
        s1 = res_path + "Cluster_%i_medoids_info.csv" % clusternum
        cluster_df.iloc[selected_medoids].to_csv(s1, index=False)
        print("Saved sources of centroids to", s1)
        s2 = res_path + "Cluster_%i_embedding.txt" % clusternum
        np.savetxt(s2, new_embedding)
        print("Saved consensus embedding to", s2)

    return [new_embedding, cluster_df, c]


# ---------------------Calculate the output value of the last shared layer consistency node - when all TCGA samples are used for model training----------------------------#
clusternum = 50
rep_num = 100
res_path = MTL_model_path + "MTL_Hot/final_embeddings/"

X_transformed_dicts = {}
for rep in np.arange(rep_num):
    X_transformed = np.loadtxt("%sMTL_Hot/final_embeddings/%i.txt" % (MTL_model_path, rep))
    X_transformed_dicts[rep] = X_transformed

res = get_last_layer_consensus_nodes(X_transformed_dicts, clusternum=clusternum, rep_num=rep_num, res_path=res_path)
