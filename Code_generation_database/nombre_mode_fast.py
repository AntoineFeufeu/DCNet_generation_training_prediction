import numpy as np
import shutil
import os
import argparse
import subprocess
from numpy import loadtxt, savetxt

## This codes aims to create folders in samples with the right number of mode in the embedding and segmentation. This is the fast one method.
## It not does a new calcul of the embedding or segmentations : It requires the mode_6 already creates and is the version which is output and active (embedding and segmentaton directly on the folders have 6 modes), and it copies the embedding and segmentation image and errase the modes that are not wanted
## So it only work for 5 or less modes
def read_args(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("--fichier_emplacement_code",
                        default="./Code_generation_database/",
                        help="where your code of simulations SPECFEM/CPS is")
    parser.add_argument("--stockage_ML_data",
                        default="./ML_DATA/",
                        help="where your originial database is")
    parser.add_argument("--nb_modes",
                        type=int,
                        default="4",
                        help="Number of modes you want to have in your data base")
    parser.add_argument("--ideb",
                        type=int,
                        default="1",
                        help="First sample")
    parser.add_argument("--ifin",
                        type=int,
                        default="21",
                        help="Last sample (not inclued)")
    args = parser.parse_args()
    return args

def dossier_nb_mode(nb_modes,stockage_data_name,stockage_data_name_mode) :
    if not os.path.isdir(stockage_data_name_mode):
        os.makedirs(stockage_data_name_mode)
    else :
        shutil.rmtree(stockage_data_name_mode)
        os.makedirs(stockage_data_name_mode)

    shutil.copyfile(stockage_data_name+"SDISPR.TXT",stockage_data_name_mode+"SDISPR.TXT")
    shutil.copyfile(stockage_data_name+"dispersion_matrix.csv",stockage_data_name_mode+"dispersion_matrix.csv")
    with open(stockage_data_name_mode+"SDISPR.TXT", "r+") as f:
        content = f.read()

        position = content.find("RAYLEIGH WAVE      MODE #  "+str(nb_modes))
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            f.seek(position-20)
            f.truncate()

def seg_emb (nb_modes,stockage_data_name) :
    stockage_data_name_mode = stockage_data_name+"mode_"+str(nb_modes)+"/"
    stockage_data_name_mode_6 = stockage_data_name+"mode_"+str(6)+"/"
    shutil.copyfile(stockage_data_name_mode_6+"embedding.csv",stockage_data_name_mode+"embedding.csv")
    shutil.copyfile(stockage_data_name_mode_6+"segmentation.csv",stockage_data_name_mode+"segmentation.csv")
    emb=loadtxt(stockage_data_name_mode+"embedding.csv",delimiter=',')
    l_emb = [] 
    for i in range (len(emb)) :
        for j in range (len(emb[0])) :
            if emb[i,j] > nb_modes :
                emb[i,j] = 0.0
                l_emb.append((i,j))
    seg=loadtxt(stockage_data_name_mode+"segmentation.csv",delimiter=',')
    for ind in l_emb :
        seg[ind[0],ind[1]] = 0.0
    # plt.imshow(emb)
    # plt.show()
    savetxt(stockage_data_name_mode+'segmentation.csv',seg,delimiter=',')
    savetxt(stockage_data_name_mode+'embedding.csv',emb,delimiter=',')

def main (args) :
    for i in range (args.ideb,args.ifin) :
        stockage_data_name = args.stockage_ML_data+str(i).zfill(4)+"/"
        stockage_data_name_mode = stockage_data_name+"mode_"+str(args.nb_modes)+"/"
        print(i)
        dossier_nb_mode(nb_modes=args.nb_modes,stockage_data_name=stockage_data_name,stockage_data_name_mode=stockage_data_name_mode)
        seg_emb(args.nb_modes,stockage_data_name=stockage_data_name)
        if os.path.exists(stockage_data_name_mode+"dispersion_matrix.csv"):
            os.remove(stockage_data_name_mode+"dispersion_matrix.csv")


if __name__ == '__main__':
    args =  read_args()
    main(args)
