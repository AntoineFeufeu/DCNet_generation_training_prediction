import numpy as np
import shutil
import os
import argparse
from numpy import loadtxt

# This code seperates the data base in two : it keep the original one with all needed to train the model or making predictions, but it moves in 
# a Stockage_output data base all the folders OUTPUT_FILES, save to create new dispersion image with the results of the receivers
# Allows to have a smaller data base, so easier to manipulate
def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stockage_ML_data",
                        default="./ML_DATA/",
                        help="where your originial database is")
    parser.add_argument("--stockage_output",
                        default="./OUTPUT_DATA/",
                        help="where your output database is")
    parser.add_argument("--ideb",
                        type=int,
                        default="10050",
                        help="First sample")
    parser.add_argument("--ifin",
                        type=int,
                        default="10051",
                        help="Last sample (not inclued)")
    args = parser.parse_args()
    return args

def transfert(stockage_output,stockage_ML_data,indice) :
    stockage_data_name = stockage_ML_data + str(indice).zfill(4) + "/OUTPUT_FILES/"
    stockage_output_name = stockage_output + str(indice).zfill(4) + "/"
    
    if os.path.isdir(stockage_data_name) :
        print(indice)
        if os.path.isdir(stockage_output_name) :
            shutil.move(stockage_data_name,stockage_output_name)
        else :
            os.mkdir(stockage_output_name)
            shutil.move(stockage_data_name,stockage_output_name)

def main (args) :
    for i in range (args.ideb,args.ifin) :
        
        transfert(stockage_output=args.stockage_output,stockage_ML_data=args.stockage_ML_data,indice=i)


if __name__ == '__main__':
    args =  read_args()
    main(args)
