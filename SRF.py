import sys
import getopt
import numpy as np
import math
from sklearn.cluster import SpectralClustering
from sklearn import linear_model

def get_args():
    geneFile   = ""  # -g --genefile
    methyFile  = ""  # -m --methyfile
    mirnaFile  = ""  # -r --mirnafile
    clusterNum = 0   # -c --clusternum

    w1 = 0
    w2 = 0
    w3 = 0

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hg:m:r:c:x:y:z:",
                                   ["help", "genefile=", "methyfile=", "mirnafile=", "clusternum=", "weight1=",
                                    "weight2=", "weight3="])
    except getopt.GetoptError:
        print(
            "Error: SRF.py -g <Gene File> -m <Methylation File> -r <miRNA File> -c <Number of Clusters> -x <Weight of Gene> -y <Weight of Methylation> -z <Weight of microRNA>")
        print()
        print(
            "   or: SRF.py --genefile <Gene File> --methyfile <Methylation File> --mirnafile <miRNA File> --clusternum <Number of Clusters> --weight1 <Weight of Gene> --weight2 <Weight of Methylation> --weight3 <Weight of microRNA>")
        sys.exit()

    for opt, arg in opts:
        if (opt not in ["-h", "--help", "-g", "--genefile", "-m", "--methyfile", "-r", "--mirnafile", "-c",
                        "--clusternum", "-x", "--weight1", "-y", "--weight2", "-z", "--weight3"]):
            print("Error:Souce Data, the number of clusters and the weight of each data types are necessary.")
            sys.exit()

    for opt, arg in opts:
        if (opt in ("-h", "--help")):
            print(
                "Try: SRF.py -g <Gene File> -m <Methylation File> -r <miRNA File> -c <Number of Clusters> -x <Weight of Gene> -y <Weight of Methylation> -z <Weight of microRNA>")
            print()
            print(
                " or: SRF.py --genefile <Gene File> --methyfile <Methylation File> --mirnafile <miRNA File> --clusternum <Number of Clusters> --weight1 <Weight of Gene> --weight2 <Weight of Methylation> --weight3 <Weight of microRNA>")
            sys.exit()
        elif (opt in ("-g", "--genefile")):
            geneFile = arg
            continue
        elif (opt in ("-m", "--methyFile")):
            methyFile = arg
            continue
        elif (opt in ("-r", "--mirnaFile")):
            mirnaFile = arg
            continue
        elif (opt in ("-c", "--clusternum")):
            clusterNum = int(arg)
            continue
        elif (opt in ("-x", "--weight1")):
            w1 = float(arg)
            continue
        elif (opt in ("-y", "--weight2")):
            w2 = float(arg)
            continue
        elif (opt in ("-z", "--weight3")):
            w3 = float(arg)
            continue
    return geneFile, methyFile, mirnaFile, clusterNum, w1, w2, w3

def SRF():

    # region hyperparameter & arguments setting

    f_hyperpara = 3
    t_hyperpara = 1
    geneFile, methyFile, mirnaFile, clusterNum, w1, w2, w3 = get_args()

    print("Calculating...")

    # endregion

    # region Read Data

    coefficient    = [w1, w2, w3]

    sourcedataFile = [geneFile, methyFile, mirnaFile]

    data_list1 = []
    data_list2 = []
    data_list3 = []
    data_list  = [data_list1, data_list2, data_list3]

    isFirstLine = True

    patientID   = []

    for i in range(3):
        with open(sourcedataFile[i]) as sourcedata:
            if (i != 0):
                next(sourcedata)
            for line in sourcedata:
                if (isFirstLine):
                    patientID   = line.split()
                    isFirstLine = False
                else:
                    row = line.split()[1:]
                    data_list[i].append(list(map(float, row)))

    Gene  = np.array(data_list[0])
    Methy = np.array(data_list[1])
    Mirna = np.array(data_list[2])

    # endregion

    # region Pearson

    pearson_Gene  = np.corrcoef(Gene , rowvar=False)
    pearson_Methy = np.corrcoef(Methy, rowvar=False)
    pearson_Mirna = np.corrcoef(Mirna, rowvar=False)

    n_patient     = len(pearson_Gene)

    # endregion

    # region Fisher transformation

    pearson_Gene_  = np.zeros(np.shape(pearson_Gene) , dtype=float)
    pearson_Methy_ = np.zeros(np.shape(pearson_Methy), dtype=float)
    pearson_Mirna_ = np.zeros(np.shape(pearson_Mirna), dtype=float)

    for i in range(n_patient):
        for j in range(i + 1, n_patient):
            if (pearson_Gene[i][j] != 1.0):
                pearson_Gene_[i][j] = pearson_Gene_[j][i] = 0.5 * (math.log((1.0 + pearson_Gene[i][j]) / (1.0 - pearson_Gene[i][j])))
            else:
                pearson_Gene_[i][j] = pearson_Gene_[j][i] = 0.5 * (math.log((1.0 + pearson_Gene[i][j]) / (1.0 - pearson_Gene[i][j] + 0.001)))

            if (pearson_Methy[i][j] != 1.0):
                pearson_Methy_[i][j] = pearson_Methy_[j][i] = 0.5 * (math.log((1.0 + pearson_Methy[i][j]) / (1.0 - pearson_Methy[i][j])))
            else:
                pearson_Methy_[i][j] = pearson_Methy_[j][i] = 0.5 * (math.log((1.0 + pearson_Methy[i][j]) / (1.0 - pearson_Methy[i][j] + 0.001)))

            if (pearson_Mirna[i][j] != 1.0):
                pearson_Mirna_[i][j] = pearson_Mirna_[j][i] = 0.5 * (math.log((1.0 + pearson_Mirna[i][j]) / (1.0 - pearson_Mirna[i][j])))
            else:
                pearson_Mirna_[i][j] = pearson_Mirna_[j][i] = 0.5 * (math.log((1.0 + pearson_Mirna[i][j]) / (1.0 - pearson_Mirna[i][j] + 0.001)))

    # endregion

    # region LinearRegression & Calculate SimilarityMatrix
    rst_coef   = np.zeros([n_patient, n_patient], dtype=list)

    rst_model  = linear_model.LinearRegression()

    Logit_P_rg = np.ones([n_patient, n_patient], dtype=float)
    Logit_P_rm = np.ones([n_patient, n_patient], dtype=float)
    Logit_P_rt = np.ones([n_patient, n_patient], dtype=float)

    ProbMat_rg = np.zeros([n_patient, n_patient], dtype=float)
    ProbMat_rm = np.zeros([n_patient, n_patient], dtype=float)
    ProbMat_rt = np.zeros([n_patient, n_patient], dtype=float)

    SimilarityMatrix = np.ones([n_patient, n_patient], dtype=float)

    for i in range(n_patient):
        for j in range(i + 1, n_patient):
            xi_gene_  = np.delete(pearson_Gene_[i] , [i, j], axis=0)
            xi_methy_ = np.delete(pearson_Methy_[i], [i, j], axis=0)
            xi_mirna_ = np.delete(pearson_Mirna_[i], [i, j], axis=0)

            xj_gene_  = np.delete(pearson_Gene_[j] , [i, j], axis=0)
            xj_methy_ = np.delete(pearson_Methy_[j], [i, j], axis=0)
            xj_mirna_ = np.delete(pearson_Mirna_[j], [i, j], axis=0)

            x_gene_   = np.append(xi_gene_ , xj_gene_)
            x_methy_  = np.append(xi_methy_, xj_methy_)
            x_mirna_  = np.append(xi_mirna_, xj_mirna_)

            # region LinearRegression: beat0 + beta1*(rm,rt,rg) + beta2*(rt,rg,rm) = (rg,rm,rt)

            # region Calculate beta

            X_rm_rt_rg_ = np.concatenate((x_methy_, x_mirna_, x_gene_))
            X_rt_rg_rm_ = np.concatenate((x_mirna_, x_gene_ , x_methy_))
            Y_rg_rm_rt_ = np.concatenate((x_gene_ , x_methy_, x_mirna_))

            X = np.array([X_rm_rt_rg_, X_rt_rg_rm_])
            X = np.transpose(X)
            Y = Y_rg_rm_rt_

            rst_result = rst_model.fit(X, Y)

            beta0 = rst_result.intercept_
            beta1 = rst_result.coef_[0]
            beta2 = rst_result.coef_[1]

            rst_coef[i][j] = [beta0, beta1, beta2]

            # endregion

            # region rg: Logit(P) = β0 + β1*rm + β2*rt = ln(P/(1-P))

            Logit_P_rg[i][j] = Logit_P_rg[j][i] \
                             = f_hyperpara * (beta0 + beta1 * pearson_Methy_[i][j] + beta2 * pearson_Mirna_[i][j])

            P_rg = math.exp(Logit_P_rg[i][j]) / (1 + math.exp(Logit_P_rg[i][j]))

            ProbMat_rg[i][j] = ProbMat_rg[j][i] = 1 / (1 + math.exp(math.pow(P_rg, -t_hyperpara)))

            # endregion

            # region rm: Logit(P) = β0 + β1*rt + β2*rg = ln(P/(1-P))

            Logit_P_rm[i][j] = Logit_P_rm[j][i] \
                             = f_hyperpara * (beta0 + beta1 * pearson_Mirna_[i][j] + beta2 * pearson_Gene_[i][j])

            P_rm = math.exp(Logit_P_rm[i][j]) / (1 + math.exp(Logit_P_rm[i][j]))

            ProbMat_rm[i][j] = ProbMat_rm[j][i] = 1 / (1 + math.exp(math.pow(P_rm, -t_hyperpara)))

            # endregion

            # region rt: Logit(P) = β0 + β1*rg + β2*rm = ln(P/(1-P))

            Logit_P_rt[i][j] = Logit_P_rt[j][i] \
                             = f_hyperpara * ( beta0 + beta1 * pearson_Gene_[i][j] + beta2 * pearson_Methy_[i][j])

            P_rt = math.exp(Logit_P_rt[i][j]) / (1 + math.exp(Logit_P_rt[i][j]))

            ProbMat_rt[i][j] = ProbMat_rt[j][i] = 1 / (1 + math.exp(math.pow(P_rt, -t_hyperpara)))

            # endregion

            # endregion

            # region SimilarityMatrix(c1*rg+c2*rm+c3*rt)

            value = (coefficient[0] * ProbMat_rg[i][j]
                     + coefficient[1] * ProbMat_rm[i][j]
                     + coefficient[2] * ProbMat_rt[i][j])
            SimilarityMatrix[i][j] = SimilarityMatrix[j][i] = value

            # endregion

    # endregion

    # region Predict Labels and Write File

    sc     = SpectralClustering(n_clusters=clusterNum, affinity='precomputed', n_jobs=-1)
    labels = sc.fit_predict(SimilarityMatrix)

    fileName = "Labels_" + str(coefficient[0]) + "_" + str(coefficient[1]) + "_" + str(coefficient[2])
    with open(fileName + ".txt", "w") as file:
        for id, lbl in zip(patientID, labels):
            file.write(id + "\t" + str(lbl) + "\n")

    print("Done!")

    # endregion

if __name__ == '__main__':
    SRF()