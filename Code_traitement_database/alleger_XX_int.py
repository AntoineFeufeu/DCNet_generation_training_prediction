import numpy as np
import os
import argparse

## Le but de ce code est d'enlever dans les fichiers OUTPUT de specfem les recepteurs selon x, inutiles pour créer l'image de dispersion

def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stockage_ML_data",
                        default="./ML_DATA/")
    parser.add_argument("--ideb",
                        type=int,
                        default="9011")
    parser.add_argument("--ifin",
                        type=int,
                        default="9046")
    parser.add_argument("--nb_receivers", 
                        type=int,
                        default="96")
    args = parser.parse_args()
    return args

def enlever_xx(stockage_ML_data,indice,nb_receivers) :
    stockage_output = stockage_ML_data + str(indice).zfill(4) + "/OUTPUT_FILES/"
    print(indice)
    for j in range (1,nb_receivers+1) :
        stockage_outpout_j = stockage_output + "AA.S"+str(j).zfill(4)+".BXX.semv"
        if os.path.exists(stockage_outpout_j):
            os.remove(stockage_outpout_j)
    
def main (args) :
    for i in range (args.ideb,args.ifin) :
        
        enlever_xx(stockage_ML_data=args.stockage_ML_data,indice=i,nb_receivers=args.nb_receivers)

if __name__ == '__main__':
    args =  read_args()
    main(args)
