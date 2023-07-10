import numpy as np
import shutil
import os
import argparse
from numpy import loadtxt
## This code allows us to activate embedding and segmentation with the right number of modes. Only works when the folder with the number of modes wanted has been generated
## Or , if you use it with the number of mode already active, it also creates the folder with this number of modes
def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stockage_ML_data",
                        default="./ML_DATA/",
                        help="where your originial database is")
    parser.add_argument("--nb_modes_voulu",
                        type=int,
                        default="4",
                        help="number of modes wanted")
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

def nb_mode_voulu(nb_modes_voulu,stockage_ML_data,indice) :
    stockage_data_name = stockage_ML_data + str(indice).zfill(4) + "/"
    print(indice)
    if os.path.exists(stockage_data_name+"embedding.csv"):
        nb_modes_actuels = int(np.max(loadtxt(stockage_data_name+"embedding.csv",delimiter=',')))
        stockage_data_name_node_actuel = stockage_data_name +"mode_"+str(nb_modes_actuels)+"/"
        if not os.path.isdir(stockage_data_name_node_actuel):
            os.makedirs(stockage_data_name_node_actuel)
    stockage_data_name_node_voulu = stockage_data_name +"mode_"+str(nb_modes_voulu)+"/"
    
    # if not os.path.isdir(stockage_data_name_node_voulu): 
    #     print("Ce nombre de mode n'a pas été généré")
    #     sys.exit()
        
    
    if os.path.exists(stockage_data_name+"embedding.csv") :
        if os.path.exists(stockage_data_name+"SDISPR.TXT"):
            shutil.copyfile(stockage_data_name+"SDISPR.TXT",stockage_data_name_node_actuel+"SDISPR.TXT")
        if os.path.exists(stockage_data_name+"embedding.csv"):
            shutil.copyfile(stockage_data_name+"embedding.csv",stockage_data_name_node_actuel+"embedding.csv")
        if os.path.exists(stockage_data_name+"embedding.jpg"):
            shutil.copyfile(stockage_data_name+"embedding.jpg",stockage_data_name_node_actuel+"embedding.jpg")
        if os.path.exists(stockage_data_name+"segmentation.csv"):
            shutil.copyfile(stockage_data_name+"segmentation.csv",stockage_data_name_node_actuel+"segmentation.csv")
        if os.path.exists(stockage_data_name+"segmentation.jpg"):
            shutil.copyfile(stockage_data_name+"segmentation.jpg",stockage_data_name_node_actuel+"segmentation.jpg")
        if os.path.exists(stockage_data_name+"superposition_labélisée.png"):
            shutil.copyfile(stockage_data_name+"superposition_labélisée.png",stockage_data_name_node_actuel+"superposition_labélisée.png")
        if os.path.exists(stockage_data_name+"superposition.jpg"):
            shutil.copyfile(stockage_data_name+"superposition.jpg",stockage_data_name_node_actuel+"superposition.jpg")

    
    if os.path.exists(stockage_data_name_node_voulu+"SDISPR.TXT"):
        shutil.copyfile(stockage_data_name_node_voulu+"SDISPR.TXT",stockage_data_name+"SDISPR.TXT")
    if os.path.exists(stockage_data_name_node_voulu+"embedding.csv"):
        shutil.copyfile(stockage_data_name_node_voulu+"embedding.csv",stockage_data_name+"embedding.csv")
    if os.path.exists(stockage_data_name_node_voulu+"embedding.jpg"):
        shutil.copyfile(stockage_data_name_node_voulu+"embedding.jpg",stockage_data_name+"embedding.jpg")
    if os.path.exists(stockage_data_name_node_voulu+"segmentation.csv"):
        shutil.copyfile(stockage_data_name_node_voulu+"segmentation.csv",stockage_data_name+"segmentation.csv")
    if os.path.exists(stockage_data_name_node_voulu+"segmentation.jpg"):
        shutil.copyfile(stockage_data_name_node_voulu+"segmentation.jpg",stockage_data_name+"segmentation.jpg")
    if os.path.exists(stockage_data_name_node_voulu+"superposition_labélisée.png"):
        shutil.copyfile(stockage_data_name_node_voulu+"superposition_labélisée.png",stockage_data_name+"superposition_labélisée.png")
    if os.path.exists(stockage_data_name_node_voulu+"superposition.jpg"):
        shutil.copyfile(stockage_data_name_node_voulu+"superposition.jpg",stockage_data_name+"superposition.jpg")

def main (args) :
    for i in range (args.ideb,args.ifin) :
        
        nb_mode_voulu(nb_modes_voulu=args.nb_modes_voulu,stockage_ML_data=args.stockage_ML_data,indice=i)


if __name__ == '__main__':
    args =  read_args()
    main(args)
