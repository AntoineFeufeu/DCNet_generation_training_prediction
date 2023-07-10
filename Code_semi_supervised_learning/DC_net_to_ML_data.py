import numpy as np
import shutil
import os
import argparse
## This code aims to transfert predictions in the data set, to make semi-supervised transfert learning
def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stockage_ML_data",
                        default="./ML_DATA/",
                        help="where your originial database is")
    parser.add_argument("--stockage_output",
                        default="./DCNet_output/",
                        help="where your predictions are stocked by your model")
    parser.add_argument("--ideb",
                        type=int,
                        default="1",
                        help="Beginning of the predictions you want to transfert")
    parser.add_argument("--ifin",
                        type=int,
                        default="3",
                        help="End (not inclued) of the predictions you want to transfert")
    args = parser.parse_args()
    return args

def transfert(stockage_output,stockage_ML_data,indice) :
    stockage_data_name = stockage_ML_data + str(indice).zfill(4) + "/"
    print(indice)
    stockage_output_name = stockage_output + str(indice).zfill(4) + "/"

    if os.path.isdir(stockage_output_name) :
         shutil.move(stockage_output_name,stockage_ML_data)
    # if os.path.exists(stockage_output_name+"embedding.csv"):
    #     shutil.move(stockage_output_name+"embedding.csv",stockage_data_name+"embedding.csv")
    # if os.path.exists(stockage_output_name+"segmentation.csv"):
    #     shutil.move(stockage_output_name+"segmentation.csv",stockage_data_name+"segmentation.csv")

def main (args) :
    for i in range (args.ideb,args.ifin) :
        
        transfert(stockage_output=args.stockage_output,stockage_ML_data=args.stockage_ML_data,indice=i)


if __name__ == '__main__':
    args =  read_args()
    main(args)
