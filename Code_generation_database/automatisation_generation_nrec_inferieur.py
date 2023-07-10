import numpy as np
import shutil
import os
import argparse
import subprocess
import random as rd
## Le but de ce code est de créer de nouvelle image de dispersion, à partir de moins de récépteurs que le nombre maximal, fixé sur la base de données à 96.

def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fichier_emplacement_code",
                        default="./Code_generation_database/",
                        help="where your code of simulations SPECFEM/CPS is")
    parser.add_argument("--stockage_data",
                        default="./ML_DATA/",
                        help="where your originial database is")
    parser.add_argument("--stockage_output",
                        default="./OUTPUT_DATA/",
                        help="where your output of specfem are stocked")
    parser.add_argument("--nb_modes",
                        type=int,
                        default="6",
                        help="Number of modes you want to have in your data base")
    parser.add_argument("--ideb",
                        type=int,
                        default="1",
                        help="The first sample that you want to make another dispersion matrix with less receivers")
    parser.add_argument("--ifin",
                        type=int,
                        default="2",
                        help="The last sample (not include) that you want to make another dispersion matrix with less receivers")
    parser.add_argument("--new_ideb",
                        type=int,
                        default="1000",
                        help="The first number add in your database for the new samples (so your new samples will be from new_ideb to new_ideb + (ifin-ideb) )")
    parser.add_argument("--nrec",
                        type=int,
                        default="19",
                        help="Number of receivers you want for your new samples")
    parser.add_argument("--dx",
                        type=float,
                        default="1",
                        help="Original distance beetween the receivers")
    
    args = parser.parse_args()
    return args

# This function creates the new emplacement for the new sample
def new_dossier(fichier_emplacement_code,i,stockage_data_name,stockage_output,new_stockage_data_name,nrec) :
    if not os.path.isdir(new_stockage_data_name):
        os.makedirs(new_stockage_data_name)
    else :
        shutil.rmtree(new_stockage_data_name)
        os.makedirs(new_stockage_data_name)
    stockage_data_name_mode = stockage_data_name+"mode_"+str(args.nb_modes)+"/"
    new_stockage_data_name_mode = new_stockage_data_name+"mode_"+str(args.nb_modes)+"/"    
    if not os.path.isdir(stockage_data_name_mode) :
        arguments = ["--fichier_emplacement_code="+fichier_emplacement_code,"--nb_modes="+str(args.nb_modes),"--ideb="+str(i),"--ifin="+str(i+1)]      
        subprocess.run(["python", fichier_emplacement_code+"nombre_mode.py"] + arguments)
        arguments = ["--nb_modes_voulu="+str(args.nb_modes),"--ideb="+str(i),"--ifin="+str(i+1)]      
        subprocess.run(["python", fichier_emplacement_code+"sortir_bon_mode.py"] + arguments)
        
    shutil.copyfile(stockage_data_name+"log",new_stockage_data_name+"log")
    shutil.copyfile(stockage_data_name+"SDISPR.TXT",new_stockage_data_name+"SDISPR.TXT")
    shutil.copyfile(stockage_data_name+"fv_boundary.csv",new_stockage_data_name+"fv_boundary.csv")
    shutil.copytree(stockage_output+"OUTPUT_FILES/",new_stockage_data_name+"OUTPUT_FILES/")
    shutil.copytree(stockage_data_name_mode,new_stockage_data_name_mode)
    if os.path.exists(new_stockage_data_name_mode+"superposition_labélisée.png") :
        os.remove(new_stockage_data_name_mode+"superposition_labélisée.png")
    if os.path.exists(new_stockage_data_name_mode+"dispersion_labélisée.png") :
        os.remove(new_stockage_data_name_mode+"dispersion_labélisée.png")
    if os.path.exists(new_stockage_data_name_mode+"superposition.jpg") :
        os.remove(new_stockage_data_name_mode+"superposition.jpg")
    with open(new_stockage_data_name+"log", 'a') as f:
        f.write('nrec = '+str(nrec))
    
# Making the dispersion matrix with the right number of receivers
def dispersion (nrec,dx,stockage_data_name,vpmin,vpmax,fmin,fmax,dt,image_dim_x,image_dim_y,nstep_output_receiver,fichier_emplacement_code,i) :
    arguments = ["--nrec="+str(nrec),"--dx="+str(dx),"--vpmin="+str(vpmin),"--vpmax="+str(vpmax),"--fmin="+str(fmin),"--fmax="+str(fmax),"--dt="+str(dt),"--image_dim_x="+str(image_dim_x),"--image_dim_y="+str(image_dim_y),"--nstep_output_receiver="+str(nstep_output_receiver),"--stockage_data_name="+stockage_data_name]      
    subprocess.run(["python", fichier_emplacement_code+"dispersion.py"] + arguments)

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
        new_stockage_data_name = args.stockage_data+str(i-args.ideb+args.new_ideb).zfill(4)+"/"
        stockage_output = args.stockage_output+str(i).zfill(4)+"/"
        fmin = float(variable_dict['fmin']) + 20
        fmax = variable_dict['fmax']
        vpmin = variable_dict['vpmin']
        vpmax = variable_dict['vpmax']
        image_dim_x = variable_dict['image_dim_x']
        image_dim_y = variable_dict['image_dim_y']
        dx = args.dx
        dt = variable_dict['dt']
        nstep_output_receiver = variable_dict['nstep_output_receiver']
        nrec = args.nrec

        print(i-args.ideb+args.new_ideb)
        new_dossier(fichier_emplacement_code=args.fichier_emplacement_code,i=i,stockage_data_name=stockage_data_name,stockage_output=stockage_output,new_stockage_data_name=new_stockage_data_name,nrec=nrec)
        dispersion (nrec,dx,new_stockage_data_name,vpmin,vpmax,fmin,fmax,dt,image_dim_x,image_dim_y,nstep_output_receiver,args.fichier_emplacement_code,i)
        shutil.rmtree(new_stockage_data_name+"OUTPUT_FILES/")

if __name__ == '__main__':
    args =  read_args()
    main(args)
