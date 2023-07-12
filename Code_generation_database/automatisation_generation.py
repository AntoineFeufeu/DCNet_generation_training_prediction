import numpy as np
import shutil
import os
import argparse
import matplotlib.pyplot as plt
import subprocess
import inspect
import random as rd

## Here is the main code for making a simulation. It uses other codes, to creates all parametric files, simulates specfem2d and CPS, creates the dispersion image, segmentation and embedding and the log file.


def read_args():
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--dt",
                        default=3e-5, #3e-5
                        type=float,
                        help="time step")
    parser.add_argument("--T",
                        default=0.4, 
                        type=float,
                        help="duration of the simulation")     
    parser.add_argument("--nstep_output_receiver",
                        default=10,
                        type=int,
                        help="Number of steps between two measures for the receivers")   
         
    parser.add_argument("--nb_layers",
                        default=28,
                        type = int,
                        help="Number of layers you want, max 28 for the default model")
    parser.add_argument("--nb_mesh_per_layer",
                        default=1,
                        type = int,
                        help="Number of mesh per layers you want (1 per meter is good)")
    
    parser.add_argument("--default_model",
                        default=True,
                        type = bool,
                        help="True=default model, False=custom model. Default one is the simplest way to make a soil : creates a soil with rho predefined, vp from 600m/s to 4000m/s, vs=vp/sqrt(3), thickness layer i = thickness**i ")    
    
    ########################################## PARAMETERS FOR DEFAULT MODEL #########################################
    parser.add_argument("--thickness",
                        default=1.0,
                        type=float,
                        help="ONLY FOR DEFAULT MODEL : thickness of the layers (1.0 -> 1.05)")
    parser.add_argument("--vp_middle",
                        default=3532,
                        type = float,
                        help="ONLY FOR DEFAULT MODEL : vp of the break (3100->3800)")
    parser.add_argument("--layer_middle",
                        default=12,
                        type = int,
                        help="ONLY FOR DEFAULT MODEL : layer of the break (11->15)")
    
    #################################################################################################################

    #### PARAMETERS FOR CUSTOM MODEL. From the highest layer to the lowest layer. Respect the length = nb_layers ####
    parser.add_argument("--custom_thickness",
                        default=[3,5,6,4,5,4,5.5,5,5],
                        nargs='+',
                        type=float,
                        help="ONLY FOR CUSTOM MODEL : thickness of the layers in a list")
    parser.add_argument("--custom_vp",
                        default=[600,1000,1500,1700,1900,2200,2300,2700,3000],
                        nargs='+',
                        type=float,
                        help="ONLY FOR CUSTOM MODEL : vp of the layers in a list")
    parser.add_argument("--custom_vs",
                        default=[400,850,1250,1400,1600,1900,1950,2200,2400],
                        nargs='+',
                        type=float,
                        help="ONLY FOR CUSTOM MODEL : vs of the layers in a list")
    parser.add_argument("--custom_rho",
                        default=[2000,2200,2250,2300,2400,2450,2500,2550,2650],
                        nargs='+',
                        type=float,
                        help="ONLY FOR CUSTOM MODEL : rho of the layers in a list")
    #################################################################################################################
    
    parser.add_argument("--f0",
                        default=80.,
                        type = float,
                        help="most present frequency (Hz)")
    parser.add_argument("--burst_band_width",
                        default=0.,
                        type = float,
                        help="bandwidth of the source")
    parser.add_argument("--nxmax",
                        default=90, 
                        type = int,
                        help="number of horizontal meshes")
    parser.add_argument("--fmax",
                        default=240,
                        type = float,
                        help="fmax")
    parser.add_argument("--fmin",
                        default=20,
                        type = float,
                        help="fmin")
    parser.add_argument("--vpmax",
                        default=2000,
                        type = float,
                        help="vpmax")
    parser.add_argument("--vpmin",
                        default=100,
                        type = float,
                        help="vpmin")
    parser.add_argument("--dx",
                        default= 1,
                        type=float,
                        help=" distance between receivers")
    parser.add_argument("--nrec",
                        default=96,
                        type = int,
                        help="nb receivers")
    parser.add_argument("--image_dim_x",
                        default=256,
                        type = int,
                        help="image_dim_x")
    parser.add_argument("--image_dim_y",
                        default=512,
                        type = int,
                        help="image_dim_y")
    
    parser.add_argument("--topographie",
                        default=True,
                        type = bool,
                        help="Do you want to input a topgraphie in your model are not")
    parser.add_argument("--linear_decay",
                        default=True,
                        type = bool,
                        help="Do you want to input a smooth topgraphie")
    parser.add_argument("--water_cavity",
                        default=True,
                        type = bool,
                        help="do you want to have a water cavity or not")
    parser.add_argument("--percent_solid",
                        default=0.0,
                        type = float,
                        help="proportion_of_solid_material_in_karst : %_gaz = 1 - %_sol - %_liq / MAX 1 MIN 0")
    parser.add_argument("--percent_liquid",
                        default=1,
                        type = float,
                        help="proportion_of_liquid_material_in_karst : %_gaz = 1 - %_sol - %_liq / MAX 1 MIN 0")
    parser.add_argument("--nxbeg_water",
                        default=20,
                        type = int,
                        help="number of horizontal meshes at the beginning of the water cavity")
    parser.add_argument("--nxend_water",
                        default=45,
                        type = int,
                        help="number of horizontal meshes at the end of the water cavity")
    parser.add_argument("--layerbeg_water",
                        default=2,
                        type = int,
                        help="layer at the beginning of the water cavity")
    parser.add_argument("--layerend_water",
                        default=4,
                        type = int,
                        help="layer at the end of the water cavity")
    
    parser.add_argument("--fichier_donnees_specfem",
                        default="./Code_generation_database/Test_src/",
                        help="Place where the initial pare_file/interface/SOURCE is stocked")
    parser.add_argument("--fichier_emplacement_programme_CPS",
                        default="./PROGRAMS.330/",
                        help="Place where the initial XDmod is stocked")
    parser.add_argument("--fichier_INPUT_initial",
                        default="./specfem2d/DATA/")
    parser.add_argument("--fichier_OUTPUT_initial",
                        default="./specfem2d/OUTPUT_FILES/")
    parser.add_argument("--fichier_emplacement_code",
                        default="./Code_generation_database/")
    parser.add_argument("--fichier_emplacement_traitement_code",
                        default="./Code_traitement_database/")
    parser.add_argument("--path_CPS",
                        default="/PROGRAMS.330/bin")
    parser.add_argument("--Stockage_data",
                    default="./ML_DATA/")
    
    
    parser.add_argument("--name",
                    default="test_custom_topo_eau",
                    help="Name of the sample")
    args = parser.parse_args()
    return args
############################################### FIN INITIALISATION ####################################################################


##This code creates the folder in the stockage_data directory, and the DATA folder in specfem2d, to make the simulation            
def creation_espaces (stockage_data_name,fichier_INPUT_initial) :
    
    if not os.path.isdir(stockage_data_name):
        os.makedirs(stockage_data_name)
    else :
        shutil.rmtree(stockage_data_name)
        os.makedirs(stockage_data_name)
        
    if not os.path.isdir(fichier_INPUT_initial):
        os.makedirs(fichier_INPUT_initial)
    else :
        shutil.rmtree(fichier_INPUT_initial)
        os.makedirs(fichier_INPUT_initial)

## This code creates all the paramteric files needed for CPS and specfem2d
def ecriture (name,fichier_emplacement_code,nb_layers,default_model,thickness,vp_middle,layer_middle,nb_mesh_per_layer,nxmax,f0,burst_band_width,fmax,dt,T,nstep_output_receiver,fichier_INPUT_initial,fichier_donnees_specfem,delta,linear_decay,water_cavity,nxbeg_water,nxend_water,layerbeg_water,layerend_water,percent_solid,percent_liquid,custom_thickness,custom_vp,custom_vs,custom_rho) :
    if default_model : 
        arguments = ["--name="+str(name),"--nb_layers="+str(nb_layers),"--default_model="+str("1" if default_model else ""),"--thickness="+str(thickness), "--vp_middle="+str(vp_middle), "--layer_middle="+str(layer_middle), "--nb_mesh_per_layer="+str(nb_mesh_per_layer),"--nxmax="+str(nxmax), "--f0="+str(f0),"--burst_band_width="+str(burst_band_width),"--dt="+str(dt),"--T="+str(T),"--nstep_output_receiver="+str(nstep_output_receiver),"--fichier_INPUT_initial="+fichier_INPUT_initial,"--fichier_donnees_specfem="+fichier_donnees_specfem,"--water_cavity="+str("1" if water_cavity else ""),"--linear_decay="+str("1" if linear_decay else ""),"--nxbeg_water="+str(nxbeg_water),"--nxend_water="+str(nxend_water),"--layerbeg_water="+str(layerbeg_water),"--layerend_water="+str(layerend_water),"--percent_solid="+str(percent_solid),"--percent_liquid="+str(percent_liquid)]
        subprocess.run(["python", fichier_emplacement_code+"ecriture_totale.py", "--delta", *map(str,delta)] + arguments)
    else :
        arguments = ["--name="+str(name),"--nb_layers="+str(nb_layers),"--default_model="+str("1" if default_model else ""),"--nb_mesh_per_layer="+str(nb_mesh_per_layer),"--nxmax="+str(nxmax), "--f0="+str(f0),"--burst_band_width="+str(burst_band_width),"--dt="+str(dt),"--T="+str(T),"--nstep_output_receiver="+str(nstep_output_receiver),"--fichier_INPUT_initial="+fichier_INPUT_initial,"--fichier_donnees_specfem="+fichier_donnees_specfem,"--water_cavity="+str("1" if water_cavity else ""),"--linear_decay="+str("1" if linear_decay else ""),"--nxbeg_water="+str(nxbeg_water),"--nxend_water="+str(nxend_water),"--layerbeg_water="+str(layerbeg_water),"--layerend_water="+str(layerend_water),"--percent_solid="+str(percent_solid),"--percent_liquid="+str(percent_liquid)]+['--custom_thickness']+[str(item) for item in custom_thickness]+['--custom_vp']+[str(item) for item in custom_vp]+['--custom_vs']+[str(item) for item in custom_vs]+['--custom_rho']+[str(item) for item in custom_rho]                                                                                                              
        subprocess.run(["python", fichier_emplacement_code+"ecriture_totale.py", "--delta", *map(str,delta)] + arguments)

## This code runs CPS and Specfem simulations
def run_CPS_and_SPECFEM (stockage_data_name,path_CPS,fichier_emplacement_programme_CPS,fichier_OUTPUT_initial,fichier_INPUT_initial) :
    print(";" + os.getcwd() + path_CPS)
    os.environ["PATH"] += ";" + os.getcwd() + path_CPS
    print(os.environ["PATH"])
    subprocess.call(['bash', fichier_emplacement_programme_CPS+'XDmod_28-layers.dep'])
    os.environ["PATH"] = os.environ["PATH"].replace(";" + os.getcwd() + path_CPS, "")
    shutil.move("SDISPR.TXT", stockage_data_name+"SDISPR.TXT")
    subprocess.call(['bash','./Code_generation_database/run_specfem.sh'])
    shutil.move(fichier_OUTPUT_initial,stockage_data_name)
    shutil.move( fichier_emplacement_programme_CPS+'XDmod_28-layers.dep',stockage_data_name)
    shutil.rmtree(fichier_INPUT_initial)
    os.remove('./2layers_orig_new.model.d')
    os.remove('./2layers_orig_new.rc.dat')
    os.remove('./2layers_orig_new.txt')
    os.remove('./dfile')
    os.remove('./sdisp96.dat')
    os.remove('./sdisp96.ray')
    os.remove('./SDISPR.PLT')

## Here is the generation of the dispersion image
def dispersion (nrec,dx,stockage_data_name,vpmin,vpmax,fmin,fmax,dt,image_dim_x,image_dim_y,nstep_output_receiver,fichier_emplacement_code) :
    arguments = ["--nrec="+str(nrec),"--dx="+str(dx),"--vpmin="+str(vpmin),"--vpmax="+str(vpmax),"--fmin="+str(fmin),"--fmax="+str(fmax),"--dt="+str(dt),"--image_dim_x="+str(image_dim_x),"--image_dim_y="+str(image_dim_y),"--nstep_output_receiver="+str(nstep_output_receiver),"--stockage_data_name="+stockage_data_name]      
    subprocess.run(["python", fichier_emplacement_code+"dispersion.py"] + arguments)
## Here is the generation of the superposition images
def superposition (stockage_data_name,vpmin,vpmax,fmin,fmax,dt,image_dim_x,image_dim_y,fichier_emplacement_code) :
    arguments = ["--vpmin="+str(vpmin),"--vpmax="+str(vpmax),"--fmin="+str(fmin),"--fmax="+str(fmax),"--image_dim_x="+str(image_dim_x),"--image_dim_y="+str(image_dim_y),"--stockage_data_name="+stockage_data_name]      
    subprocess.run(["python", fichier_emplacement_code+"superposition_CPS_specfem.py"] + arguments)
## Here is the generation of the segmentation and embedding
def seg_emb (stockage_data_name,vpmin,vpmax,fmin,fmax,dt,image_dim_x,image_dim_y,fichier_emplacement_code) :
    arguments = ["--vpmin="+str(vpmin),"--vpmax="+str(vpmax),"--fmin="+str(fmin),"--fmax="+str(fmax),"--image_dim_x="+str(image_dim_x),"--image_dim_y="+str(image_dim_y),"--stockage_data_name="+stockage_data_name]      
    subprocess.run(["python", fichier_emplacement_code+"segmentation-embedding.py"] + arguments)
## Here is the generation of the log file, where all the parameters used are written
def log (stockage_data_name,thickness,nb_layers,vp_middle,layer_middle,dx,f0,burst_band_width,vpmin,vpmax,fmin,fmax,dt,T,image_dim_x,image_dim_y,nstep_output_receiver,topographie,delta,water_cavity,nxbeg_water,nxend_water,layerbeg_water,layerend_water,percent_solid,percent_liquid) :
    path_read = stockage_data_name+"log"
    var_dict = locals()
    if os.path.isdir(path_read):
        os.remove(path_read)
    argspec = inspect.getfullargspec(log)

    arg_names = argspec.args
    with open(path_read, "w") as f:
        for i in arg_names :
            if i != 'stockage_data_name' :
                f.write(str(i) + " = "+str(var_dict[i])+"\n")
                
## Here is the code who delete XX receivers output, which are not used.
def alleger_xx(fichier_emplacement_code,stockage_data_name) :
    arguments = ["--stockage_data_name="+stockage_data_name]      
    subprocess.run(["python", fichier_emplacement_code+"alleger_XX.py"] + arguments)

def main(args):
    global delta
    ## Adding topographic feature or not
    if not args.topographie :
        delta = [0,0]
    else :
        # Choose how to generate the delta  for topography here
        delta= [0.]
        for j in range (13) :
            delta.append(abs(delta[-1]+2*(rd.random()-0.5)))
    delta = [-1.3501874999998336, -1.3071874999998272, -1.3501874999998336, -1.3621874999998909, -1.4311874999998508, -1.3841874999998254, -1.4741874999998572, -1.3991874999998117, -1.4551874999998518, -1.508187499999849, -1.5381874999998217, -1.5871874999999136, -1.6521874999998545, -1.7641874999998208, -1.8471874999999045, -2.0591874999998936, -1.897187499999859, -2.14318749999984, -2.2941874999999072, -2.1861874999998463, -2.0251874999999018, -2.414187499999912, -2.4301874999998745, -2.4381874999999127, -2.1911874999998417, -2.13018749999992, -1.8431874999998854, -1.5821874999999181, -1.4601874999998472, -1.1911874999998417, -1.081187499999828, -0.7711874999998827, -0.61618749999991, -0.5291874999999209, -0.35018749999983356, -0.2101874999998472, 0.07381250000014461, 0.3448125000001028, 0.5938125000001264, 0.6278125000001182, 0.6778125000001864, 0.7888125000001764, 0.9978125000001228, 1.2718125000001237, 1.557812500000182, 1.7748125000001664, 1.9838125000001128, 2.0668125000000828, 2.171812500000101, 2.195812500000102, 2.428812500000163, 2.49481250000008, 2.6908125000001064, 2.6478125000001, 2.637812500000109, 2.6838125000001583, 2.588812500000131, 2.775812500000143, 2.7748125000001664, 2.8558125000001837, 2.984812500000089, 3.1018125000001646, 3.2168125000001737, 3.24481250000008]
    stockage_data_name = args.Stockage_data + args.name + "/"
    creation_espaces (stockage_data_name,args.fichier_INPUT_initial)
    ecriture (args.name,args.fichier_emplacement_code,args.nb_layers,args.default_model,args.thickness,args.vp_middle,args.layer_middle,args.nb_mesh_per_layer,args.nxmax,args.f0,args.burst_band_width,args.fmax,args.dt,args.T,args.nstep_output_receiver,args.fichier_INPUT_initial,args.fichier_donnees_specfem,delta,args.linear_decay,args.water_cavity,args.nxbeg_water,args.nxend_water,args.layerbeg_water,args.layerend_water,percent_solid=args.percent_solid,percent_liquid=args.percent_liquid,custom_thickness=args.custom_thickness,custom_vp=args.custom_vp,custom_vs=args.custom_vs,custom_rho=args.custom_rho)
    run_CPS_and_SPECFEM (stockage_data_name,args.path_CPS,args.fichier_emplacement_programme_CPS,args.fichier_OUTPUT_initial,args.fichier_INPUT_initial)
    dispersion(args.nrec,args.dx,stockage_data_name,args.vpmin,args.vpmax,args.fmin,args.fmax,args.dt,args.image_dim_x,args.image_dim_y,args.nstep_output_receiver,args.fichier_emplacement_code)
    superposition(stockage_data_name,args.vpmin,args.vpmax,args.fmin,args.fmax,args.dt,args.image_dim_x,args.image_dim_y,args.fichier_emplacement_code)
    seg_emb(stockage_data_name,args.vpmin,args.vpmax,args.fmin,args.fmax,args.dt,args.image_dim_x,args.image_dim_y,args.fichier_emplacement_code)
    log (stockage_data_name,args.thickness,args.nb_layers,args.vp_middle,args.layer_middle,args.dx,args.f0,args.burst_band_width,args.vpmin,args.vpmax,args.fmin,args.fmax,args.dt,args.T,args.image_dim_x,args.image_dim_y,args.nstep_output_receiver,args.topographie,delta,args.water_cavity,args.nxbeg_water,args.nxend_water,args.layerbeg_water,args.layerend_water,args.percent_solid,args.percent_liquid)
    #alleger_xx(fichier_emplacement_code=args.fichier_emplacement_code,stockage_data_name=stockage_data_name)
if __name__ == '__main__':
    args =  read_args()
    main(args)
