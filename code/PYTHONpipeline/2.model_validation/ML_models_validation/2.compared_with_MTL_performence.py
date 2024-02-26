# 提取在三个验证集中各个算法的性能-pearson-corr
import pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from configs import *


exclude_data_list = ["ExpO_project", "Braun.NatMed.2021", "Mariathasan.Nature.2018"]

for dset in exclude_data_list:
    print(dset)
    per_validate_res_path = validate_res_path + "%s/" % dset

    MTL_performance = pickle.load(open(per_validate_res_path + "rescaled_final_model/0/performances.p", "rb"))
    pearson_corr = {}
    for var in phenotypes:
        pearson_corr[var] = []
        pearson_corr[var].append(MTL_performance[var]['pearson_corr'])
        for method in ["ridge", "lasso", "ElasticNet", "SVM_linear", "SVM_RBF", "RandomForest"]:
            per_ML_performance = pickle.load(open(per_validate_res_path + "rescaled_final_ML_model/%s/performances.p" % method, "rb"))
            pearson_corr[var].append(per_ML_performance[var]['pearson_corr'])

    # 绘制折线图
    with PdfPages(per_validate_res_path + 'models_compare_corr_line_plot.pdf') as pdf:
        plt.title(dset)
        x_axix = ["MTL", "ridge", "lasso", "ElasticNet", "SVM_linear", "SVM_RBF", "RandomForest"]
        plt.grid()
        plt.plot(x_axix, pearson_corr[phenotypes[0]], color='green', label=phenotypes[0])
        plt.plot(x_axix, pearson_corr[phenotypes[1]], color='red', label=phenotypes[1])
        plt.plot(x_axix, pearson_corr[phenotypes[2]], color='skyblue', label=phenotypes[2])
        plt.plot(x_axix, pearson_corr[phenotypes[3]], color='blue', label=phenotypes[3])
        plt.legend()  # 显示图例
        plt.xlabel('')
        plt.ylabel('pearson correlation coefficient')
        pdf.savefig()  # saves the current figure into a pdf page
        plt.show()
        plt.close()



