import numpy as np
import shutil
import os
import argparse
import subprocess

## This codes aims to create folders in samples with the right number of mode in the embedding and segmentation. This is the vlassic method.
## It  does a new calcul of the embedding or segmentations : It NOT requires the mode_6 
## So it only work for 6 or less modes
def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fichier_emplacement_code",
                        default="./Code_generation_database/",
                        help="where your code of simulations SPECFEM/CPS is")
    parser.add_argument("--stockage_data",
                        default="./ML_DATA/",
                        help="where your originial database is")
    parser.add_argument("--nb_modes",
                        type=int,
                        default="3",
                        help="Number of modes you want to have in your data base")
    parser.add_argument("--ideb",
                        type=int,
                        default="1",
                        help="First sample")
    parser.add_argument("--ifin",
                        type=int,
                        default="3",
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
    
    
def superposition (stockage_data_name_mode,vpmin,vpmax,fmin,fmax,image_dim_x,image_dim_y,fichier_emplacement_code) :
    arguments = ["--vpmin="+str(vpmin),"--vpmax="+str(vpmax),"--fmin="+str(fmin),"--fmax="+str(fmax),"--image_dim_x="+str(image_dim_x),"--image_dim_y="+str(image_dim_y),"--stockage_data_name="+stockage_data_name_mode]      
    subprocess.run(["python", fichier_emplacement_code+"superposition_CPS_specfem.py"] + arguments)

def seg_emb (nb_modes,stockage_data_name_mode,vpmin,vpmax,fmin,fmax,image_dim_x,image_dim_y,fichier_emplacement_code) :
    arguments = ["--modes_number="+str(nb_modes),"--vpmin="+str(vpmin),"--vpmax="+str(vpmax),"--fmin="+str(fmin),"--fmax="+str(fmax),"--image_dim_x="+str(image_dim_x),"--image_dim_y="+str(image_dim_y),"--stockage_data_name="+stockage_data_name_mode]      
    subprocess.run(["python", fichier_emplacement_code+"segmentation-embedding.py"] + arguments)
    

def main (args) :
    for i in range (args.ideb,args.ifin) :
        
        #Algorithm which read the log file and create a dictionnary with the parameters
        with open(args.stockage_data+str(i).zfill(4)+"/log", 'r') as file:
            content = file.read()
        lines = content.split('\n')
        variable_dict = {}
        for line in lines:
            parts = line.split('=')
            if len(parts) != 2:
                continue
            name = parts[0].strip()
            value = parts[1].strip()
            variable_dict[name] = value
        stockage_data_name = args.stockage_data+str(i).zfill(4)+"/"
        stockage_data_name_mode = stockage_data_name+"mode_"+str(args.nb_modes)+"/"
        
        fmin = variable_dict['fmin']
        fmax = variable_dict['fmax']
        vpmin = variable_dict['vpmin']
        vpmax = variable_dict['vpmax']
        image_dim_x = variable_dict['image_dim_x']
        image_dim_y = variable_dict['image_dim_y']

        if not os.path.isdir(stockage_data_name_mode) :
            print(i)
            dossier_nb_mode(nb_modes=args.nb_modes,stockage_data_name=stockage_data_name,stockage_data_name_mode=stockage_data_name_mode)
            superposition(stockage_data_name_mode,vpmin,vpmax,fmin,fmax,image_dim_x,image_dim_y,args.fichier_emplacement_code)
            seg_emb(args.nb_modes,stockage_data_name_mode,vpmin,vpmax,fmin,fmax,image_dim_x,image_dim_y,args.fichier_emplacement_code)
            if os.path.exists(stockage_data_name_mode+"dispersion_matrix.csv"):
                os.remove(stockage_data_name_mode+"dispersion_matrix.csv")


if __name__ == '__main__':
    args =  read_args()
    main(args)
