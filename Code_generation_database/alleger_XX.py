import numpy as np
import os
import argparse

## Le but de ce code est d'enlever dans les fichiers OUTPUT de specfem les recepteurs selon x, inutiles pour créer l'image de dispersion

def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stockage_data_name",
                        default="./ML_DATA/")
    parser.add_argument("--nb_rec", 
                        type=int,
                        default="96")
    args = parser.parse_args()
    return args

def enlever_xx(stockage_data_name,nb_rec) :
    stockage_output = stockage_data_name + "/OUTPUT_FILES/"
    for j in range (1,nb_rec+1) :
        stockage_outpout_j = stockage_output + "AA.S"+str(j).zfill(4)+".BXX.semv"
        if os.path.exists(stockage_outpout_j):
            os.remove(stockage_outpout_j)
    
def main (args) :     
    enlever_xx(stockage_data_name=args.stockage_data_name,nb_rec=args.nb_rec)

if __name__ == '__main__':
    args =  read_args()
    main(args)
