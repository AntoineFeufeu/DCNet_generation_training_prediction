import numpy as np
import shutil
import os
import argparse
import matplotlib.pyplot as plt
import random as rd

## This code aims to create all parametric files for CPS and specfem simulations
def read_args():
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--dt",
                        default=3e-5,
                        type=float,
                        help="time step")
    parser.add_argument("--T",
                        default=0.4,
                        type=float,
                        help="duration of the simulation")     
    parser.add_argument("--nstep_output_receiver",
                        default=1,
                        type=int,
                        help="Number of steps between two measures for the receivers")            
    parser.add_argument("--thickness",
                        default=1.01,
                        type=float,
                        help="thickness of the layers (1.01 -> 1.05)")
    parser.add_argument("--nb_layers",
                        default=28,
                        type = int,
                        help="Number of layers you want, max 28")
    parser.add_argument("--nb_mesh_per_layer",
                        default=1,
                        type = int,
                        help="Number of mesh per layers you want")
    parser.add_argument("--default_model",
                        default=True,
                        type = bool,
                        help="True=default model, False=custom model. Default one is the simplest way to make a soil : creates a soil with rho predefined, vp from 600m/s to 4000m/s, vs=vp/sqrt(3), thickness layer i = thickness**i ")    
    parser.add_argument("--vp_middle",
                        default=3100.,
                        type = float,
                        help="vp of the break (3100->3800)")
    parser.add_argument("--layer_middle",
                        default=11,
                        type = int,
                        help="layer of the break (11->15)")
    parser.add_argument("--custom_thickness",
                        default=[3,5,6,4,5],
                        nargs='+',
                        type=float,
                        help="ONLY FOR CUSTOM MODEL : thickness of the layers in a list")
    parser.add_argument("--custom_vp",
                        default=[600,1000,1500,1700,1900],
                        nargs='+',
                        type=float,
                        help="ONLY FOR CUSTOM MODEL : vp of the layers in a list")
    parser.add_argument("--custom_vs",
                        default=[400,850,1250,1400,1600],
                        nargs='+',
                        type=float,
                        help="ONLY FOR CUSTOM MODEL : vs of the layers in a list")
    parser.add_argument("--custom_rho",
                        default=[2000,2200,2250,2300,2400],
                        nargs='+',
                        type=float,
                        help="ONLY FOR CUSTOM MODEL : rho of the layers in a list")
    parser.add_argument("--f0",
                        default=80.,
                        type = float,
                        help="most present frequency")
    parser.add_argument("--burst_band_width",
                        default=0.,
                        type = float,
                        help="bandwidth of the source")
    parser.add_argument("--xdeb",
                        default=0,
                        type = float,
                        help="horizontal beginning of the simulation window")
    parser.add_argument("--xfin",
                        default=120,
                        type = float,
                        help="horizontal end of the simulation window")
    parser.add_argument("--nxmax",
                        default=90,
                        type = int,
                        help="number of horizontal meshes")
    parser.add_argument("--delta",
                        default=[0.,0.],
                        nargs = '+',
                        type = float,
                        help="topography model")
    parser.add_argument("--water_cavity",
                        default=False,
                        type = bool,
                        help="do you want to have a water cavity or not")
    parser.add_argument("--linear_decay",
                        default=True,
                        type = bool,
                        help="Do you want to input a smooth topgraphie")
    parser.add_argument("--percent_solid",
                        default=65,
                        type = float,
                        help="%_of_solid_material_in_karst : %_gaz = 1 - %_sol - %_liq")
    parser.add_argument("--percent_liquid",
                        default=30,
                        type = float,
                        help="%_of_liquid_material_in_karst : %_gaz = 1 - %_sol - %_liq ")
    parser.add_argument("--nxbeg_water",
                        default=70,
                        type = int,
                        help="number of horizontal meshes at the beginning of the water cavity")
    parser.add_argument("--nxend_water",
                        default=80,
                        type = int,
                        help="number of horizontal meshes at the end of the water cavity")
    parser.add_argument("--layerbeg_water",
                        default=10,
                        type = int,
                        help="layer at the beginning of the water cavity")
    parser.add_argument("--layerend_water",
                        default=12,
                        type = int,
                        help="layer at the end of the water cavity")
    parser.add_argument("--fichier_donnees_specfem",
                        default="./Code_generation_database/Test_src/",
                        help="Place where the initial pare_file/interface/SOURCE is stocked")
    parser.add_argument("--fichier_emplacement_programme_CPS",
                        default="./PROGRAMS.330/",
                        help="Place where the initial XDmod is stocked")
    parser.add_argument("--fichier_INPUT_initial",
                        default="./specfem2d/DATA/",
                        help="Folder where you want to stock the files generated")
    parser.add_argument("--name",
                        default="0001",
                        help="Folder name of the sample")
    args = parser.parse_args()
    return args


###################################### MODIFICATION PAR_FILE ####################################################################

def ecriture_Par_file (fichier_donnees_specfem,fichier_INPUT_initial,T,dt,nstep_output_receiver,nb_layers,nxmax,material_number,rho,vp,vs,nxmin,nzmin,nzmax,water_cavity,nxbeg_water,nxend_water,layerbeg_water,layerend_water,percent_solid,percent_liquid) :
    fichier_parfile_initial = fichier_donnees_specfem+"Par_file"
    fichier_parfile_complete = fichier_INPUT_initial+"Par_file"
    shutil.copyfile(fichier_parfile_initial, fichier_parfile_complete)
    with open(fichier_parfile_complete, "r+") as f:
        content = f.read()

        position = content.find("# total number of time steps")
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            f.seek(position + len("# total number of time steps"))
            f.write("\nNSTEP                           = "+str(int(T/dt)))
        position = content.find("# duration of a time step")
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            f.seek(position + len("# duration of a time step"))
            f.write("\nDT                              = "+str(dt))

        position = content.find("# defaults to 1, which means no down-sampling")
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            f.seek(position + len("# defaults to 1, which means no down-sampling"))
            f.write("\nNTSTEP_BETWEEN_OUTPUT_SAMPLE    = "+str(nstep_output_receiver))
            
        position = content.find("# number of model materials")
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            f.seek(position + len("# number of model materials"))
            f.write("\n nbmodels                        = "+str(nb_layers+1)+" ")
        position = content.find("#       utils/attenuation/conversion_from_Qkappa_Qmu_to_Qp_Qs_from_Dahlen_Tromp_959_960.f90.")
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            percent_gaz = 1 - percent_solid - percent_solid
            rho_sol = rho[int((layerend_water+layerbeg_water)/2)]
            vp_sol = vp[int((layerend_water+layerbeg_water)/2)]
            vs_sol = vs[int((layerend_water+layerbeg_water)/2)]
            rho_karst = rho_sol*percent_solid+1000*percent_liquid+1.292*percent_gaz
            f.seek(position + len("#       utils/attenuation/conversion_from_Qkappa_Qmu_to_Qp_Qs_from_Dahlen_Tromp_959_960.f90."))
            for i in range (nb_layers):
                f.write("\n "+str(material_number[i])+" "+"1"+" "+str(rho[i])+" "+str(vp[i])+" "+str(vs[i])+" 0 0 9999 9999 0 0 0 0 0 0")
            f.write("\n "+str(material_number[-1]+1)+" "+"1"+" "+str(rho_karst)+" "+str((percent_solid*vp_sol*rho_sol+percent_liquid*1500*1000+percent_gaz*340*1.292)/rho_karst)+" "+str(percent_solid*vs_sol*rho_sol/rho_karst)+" 0 0 9999 9999 0 0 0 0 0 0")
            f.write("\n")
        position = content.find("# RENTRER nx")
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            f.seek(position + len("# RENTRER nx"))
            f.write("\n nx                              = "+str(nxmax)+" ")
        position = content.find("# define the different regions of the model in the (nx,nz) spectral-element mesh")
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            f.seek(position + len("# define the different regions of the model in the (nx,nz) spectral-element mesh"))
            if not water_cavity :
                f.write("\n nbregions                       = "+str(nb_layers+1)+" ")
            elif nxmin <= nxbeg_water-1 and nxend_water+1 <= nxmax :
                f.write("\n nbregions                       = "+str(nb_layers+1+2*(layerend_water-layerbeg_water+1))+" ")
            elif nxmin <= nxbeg_water-1 or nxend_water+1 <= nxmax :
                f.write("\n nbregions                       = "+str(nb_layers+1+1*(layerend_water-layerbeg_water+1))+" ")
            else : 
                f.write("\n nbregions                       = "+str(nb_layers+1))
        position = content.find("# format of each line: nxmin nxmax nzmin nzmax material_number")
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            f.seek(position + len("# format of each line: nxmin nxmax nzmin nzmax material_number"))
            f.write("\n "+str(nxmin)+" "+str(nxmax)+" "+str(1)+" "+str(2)+" "+str(1))
            if not water_cavity :
                for i in range (nb_layers) :
                    f.write("\n "+str(nxmin)+" "+str(nxmax)+" "+str(nzmin[i]+2)+" "+str(nzmax[i]+2)+" "+str(material_number[i]))
                f.write("\n")
            else :
                for i in range (nb_layers) :
                    if (not i >= layerbeg_water) or (not i <= layerend_water) :
                        f.write("\n "+str(nxmin)+" "+str(nxmax)+" "+str(nzmin[i]+2)+" "+str(nzmax[i]+2)+" "+str(material_number[i]))
                    else :
                        if nxmin <= nxbeg_water-1 :
                            f.write("\n "+str(nxmin)+" "+str(nxbeg_water-1)+" "+str(nzmin[i]+2)+" "+str(nzmax[i]+2)+" "+str(material_number[i]))
                        f.write("\n "+str(nxbeg_water)+" "+str(nxend_water)+" "+str(nzmin[i]+2)+" "+str(nzmax[i]+2)+" "+str(material_number[-1]+1))
                        if nxend_water+1 <= nxmax :
                            f.write("\n "+str(nxend_water+1)+" "+str(nxmax)+" "+str(nzmin[i]+2)+" "+str(nzmax[i]+2)+" "+str(material_number[i]))
                f.write("\n")

###################################### FIN MODIFICATION PAR_FILE ###############################################################



###################################### MODIFICATION INTERFACE ###############################################################
def ecriture_interfaces (fichier_donnees_specfem,linear_decay,fichier_INPUT_initial,profondeur_interface,nb_layers,xdeb,xfin,nb_mesh_per_layer,delta) :
    delta2 = list(np.zeros(len(delta)))
    #delta2=list(np.linspace(0,-3,len(delta)))
    n_interface_diff = 26
    fichier_interface_initial = fichier_donnees_specfem+"interfaces.dat"
    fichier_interface_complete = fichier_INPUT_initial+"interfaces.dat"
    shutil.copyfile(fichier_interface_initial, fichier_interface_complete)
    n = len(delta)
    with open(fichier_interface_complete, "r+") as f:
        content = f.read()
        position = content.find("# number of interfaces")
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            f.seek(position + len("# number of interfaces"))
            f.write("\n "+str(nb_layers+2)+" ")
            
        position = content.find("# for each interface below, we give the number of points and then x,z for each point\n#")
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            f.seek(position + len("# for each interface below, we give the number of points and then x,z for each point\n#"))
            f.write("\n#\n# interface number "+str(1)+"\n# ")
            f.write("\n 2")
            f.write("\n "+str(xdeb)+" "+str(round(profondeur_interface[0]-2,11)))
            f.write("\n "+str(xfin)+" "+str(round(profondeur_interface[0]-2,11)))
            for i in range (1,nb_layers+2):
                f.write("\n#\n# interface number "+str(i+1)+"\n# ")
                f.write("\n "+str(n))
                if linear_decay :
                    for j in range(n) :
                        if i <= n_interface_diff :
                            f.write("\n "+str((xfin-xdeb)*j/(n-1)+xdeb)+" "+str(round(profondeur_interface[i-1]+(delta[j]*i)/(nb_layers+1)+(delta2[j]*i)/(n_interface_diff),11)))
                        else :
                            f.write("\n "+str((xfin-xdeb)*j/(n-1)+xdeb)+" "+str(round(profondeur_interface[i-1]+(delta[j]*i)/(nb_layers+1),11)))
                else : 
                    for j in range(n) :
                        f.write("\n "+str((xfin-xdeb)*j/(n-1)+xdeb)+" "+str(round(profondeur_interface[i-1]+(delta[j])/(nb_layers+2-i),11)))
        position = content.find("# for each layer, we give the number of spectral elements in the vertical direction\n#")
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            f.seek(position + len("# for each layer, we give the number of spectral elements in the vertical direction\n#"))
            f.write("\n#\n# layer number "+str(1)+"\n# ")
            f.write("\n "+str(2))
            for i in range (1,nb_layers+1):
                f.write("\n#\n# layer number "+str(i+1)+"\n# ")
                f.write("\n "+str(nb_mesh_per_layer))

###################################### FIN MODIFICATION INTERFACE ###############################################################




###################################### MODIFICATION SOURCE ###############################################################
def ecriture_SOURCE (fichier_donnees_specfem,fichier_INPUT_initial,f0,burst_band_width,xdeb,xfin,delta,name) :
    fichier_source_initial = fichier_donnees_specfem+"SOURCE"
    fichier_source_complete = fichier_INPUT_initial+"SOURCE"
    n = len(delta)
    dx = (xfin-xdeb)/(n-1)
    ndx = int(9//dx)
    r = 9%dx
    #######
    zsource = delta[ndx] + (delta[ndx+1]-delta[ndx])/dx * r
    shutil.copyfile(fichier_source_initial, fichier_source_complete)
    
    with open(fichier_source_complete, "r+") as f:
        content = f.read()
        position = content.find("# Zs a completer")
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            f.seek(position + len("# Zs a completer"))
            f.write("\nzs                              = "+str(zsource))
        position = content.find("#BANDWIDTH TO CHANGE")
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            f.seek(position + len("#BANDWIDTH TO CHANGE"))
            f.write("\nburst_band_width                = "+str(burst_band_width)+" ")
        position = content.find("#f0 TO CHANGE")
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            f.seek(position + len("#f0 TO CHANGE"))
            f.write("\nf0                              = "+str(f0)+" ")

###################################### FIN MODIFICATION SOURCE ###############################################################

def ecriture_XDmod (fichier_donnees_specfem,fichier_INPUT_initial,profondeur_interface,nb_layers,rho,vp,vs,fichier_emplacement_programme_CPS,water_cavity,layerbeg_water,layerend_water,percent_solid,percent_liquid) :
    fichier_CPS_initial = fichier_donnees_specfem+"XDmod_28-layers.dep"
    fichier_CPS_util = fichier_emplacement_programme_CPS+"XDmod_28-layers.dep"
    shutil.copyfile(fichier_CPS_initial,fichier_CPS_util)
    
    with open(fichier_CPS_util, "r+") as f:
        content = f.read()
        position = content.find("HR      VP      VS     RHO QP  QS  ETAP ETAS FREFP FREFS")
        if position == -1:
            print("Le texte à modifier n'a pas été trouvé.")
        else:
            f.seek(position + len("HR      VP      VS     RHO QP  QS  ETAP ETAS FREFP FREFS"))
            if not water_cavity :
                for i in range (nb_layers-1,-1,-1) :
                    f.write("\n"+str((profondeur_interface[i+1]-profondeur_interface[i])/1000)+" "+str(vp[i]/1000)+" "+str(vs[i]/1000)+" "+str(rho[i]/1000000)+" "+" 0 0 0 0 1 1")
            else :
                percent_gaz = 1 - percent_solid - percent_solid
                rho_sol = rho[int((layerend_water+layerbeg_water)/2)]
                vp_sol = vp[int((layerend_water+layerbeg_water)/2)]
                vs_sol = vs[int((layerend_water+layerbeg_water)/2)]
                rho_karst = rho_sol*percent_solid+1000*percent_liquid+1.292*percent_gaz
                for i in range (nb_layers-1,-1,-1) :
                    if (not (i) >= layerbeg_water) or ((not (i) <= layerend_water)) :
                        f.write("\n"+str((profondeur_interface[i+1]-profondeur_interface[i])/1000)+" "+str(vp[i]/1000)+" "+str(vs[i]/1000)+" "+str(rho[i]/1000000)+" "+" 0 0 0 0 1 1")
                    else :
                        f.write("\n"+str((profondeur_interface[i+1]-profondeur_interface[i])/1000)+" "+str(((percent_solid*vp_sol*rho_sol+percent_liquid*1500*1000+percent_gaz*340*1.292)/rho_karst)/1000)+" "+str((percent_solid*vs_sol*rho_sol/rho_karst)/1000)+" "+str(rho_karst/1000000)+" "+" 0 0 0 0 1 1")
            f.write("\n0.00000000      1.20000      0.690000      2.00000   0.000000    0.000000       0.0000000       0.00000000       1.00000000       1.00000000")

    ## Q v(f ) = Qv(f / fv) ^ η v 
     
    
def creation_parameters (nb_layers,thickness,nb_mesh_per_layer,layer_middle,vp_middle) :
    global rho
    global vp
    global vs
    global material_number
    global nxmin
    global nzmin
    global nzmax
    global epaisseur
    global profondeur_interface
    material_number = [i for i in range(1,nb_layers+1)]
    rho = [2449.0, 2456., 2473., 2488., 2487., 2475., 2464., 2457., 2451., 2445., 2435., 2424., 2411., 2393., 2377., 2361., 2343., 2322, 2295., 2253., 2191., 2085., 1958., 1865., 1764., 1659., 1580., 1558.]
    nxmin = 1
    nzmin = []
    nzmax=[]
    vp=[]
    vs=[]
    epaisseur = 0
    for i in range(nb_layers) :
        epaisseur+= thickness**(i+1)
        nzmin.append(1+nb_mesh_per_layer*i)
        nzmax.append(nb_mesh_per_layer*(i+1))
    profondeur_interface = [-epaisseur]
    for i in range(nb_layers) :
        profondeur_interface.append(profondeur_interface[-1]+thickness**(nb_layers-i))
    for i in range(28) :
        if ((i+1)<layer_middle) :
            vp.append(4000-((4000-vp_middle)/(layer_middle-1))*i)
        elif ((i+1) == layer_middle) :
            vp.append(vp_middle)
        else :
            vp.append(vp_middle-((vp_middle-600)/(28-layer_middle))*(i-layer_middle+1))
        vs.append(vp[-1]/np.sqrt(3))
    vp = (np.array(vp))[(28-nb_layers):]
    vs = (np.array(vs))[(28-nb_layers):]
    rho = (np.array(rho))[(28-nb_layers):]

def creation_parameters_custom (custom_thickness,custom_vp,custom_vs,custom_rho,nb_layers,nb_mesh_per_layer) :
    global rho
    global vp
    global vs
    global material_number
    global nxmin
    global nzmin
    global nzmax
    global epaisseur
    global profondeur_interface
    material_number = [i for i in range(1,nb_layers+1)]
    nxmin = 1
    nzmin = []
    nzmax=[]
    epaisseur = 0
    for i in range(nb_layers) :
        epaisseur+= custom_thickness[i]
        nzmin.append(1+nb_mesh_per_layer*i)
        nzmax.append(nb_mesh_per_layer*(i+1))
    profondeur_interface = [-epaisseur]
    for i in range(nb_layers) :
        profondeur_interface.append(profondeur_interface[-1]+custom_thickness[nb_layers-i-1])
    vp=custom_vp
    vp.reverse()
    vs=custom_vs
    vs.reverse()
    rho=custom_rho
    rho.reverse()
    print(vp,vs,rho,profondeur_interface)
    
def main(args):
    print(args.custom_vp,args.custom_vs,args.custom_rho,args.custom_thickness,args.default_model)
    if args.default_model :
        creation_parameters (nb_layers=args.nb_layers,thickness=args.thickness,nb_mesh_per_layer=args.nb_mesh_per_layer,layer_middle=args.layer_middle,vp_middle=args.vp_middle)
    else :
        creation_parameters_custom (custom_thickness=args.custom_thickness,custom_vp=args.custom_vp,custom_vs=args.custom_vs,custom_rho=args.custom_rho,nb_layers=args.nb_layers,nb_mesh_per_layer=args.nb_mesh_per_layer)
    ecriture_Par_file(fichier_donnees_specfem=args.fichier_donnees_specfem,fichier_INPUT_initial=args.fichier_INPUT_initial,T=args.T,dt=args.dt,nstep_output_receiver=args.nstep_output_receiver,nb_layers=args.nb_layers,nxmax=args.nxmax,material_number=material_number,rho=rho,vp=vp,vs=vs,nxmin=nxmin,nzmin=nzmin,nzmax=nzmax,water_cavity=args.water_cavity,nxbeg_water=args.nxbeg_water,nxend_water=args.nxend_water,layerbeg_water=args.layerbeg_water,layerend_water=args.layerend_water,percent_solid=args.percent_solid,percent_liquid=args.percent_liquid)
    ecriture_interfaces(fichier_donnees_specfem=args.fichier_donnees_specfem,linear_decay=args.linear_decay,fichier_INPUT_initial=args.fichier_INPUT_initial,profondeur_interface=profondeur_interface,nb_layers=args.nb_layers,xdeb=args.xdeb,xfin=args.xfin,nb_mesh_per_layer=args.nb_mesh_per_layer,delta=args.delta)
    ecriture_SOURCE(fichier_donnees_specfem=args.fichier_donnees_specfem,fichier_INPUT_initial=args.fichier_INPUT_initial,f0=args.f0,burst_band_width=args.burst_band_width,xdeb=args.xdeb,xfin=args.xfin,delta=args.delta,name=args.name)
    ecriture_XDmod(fichier_donnees_specfem=args.fichier_donnees_specfem,fichier_INPUT_initial=args.fichier_INPUT_initial,profondeur_interface=profondeur_interface,nb_layers=args.nb_layers,rho=rho,vp=vp,vs=vs,fichier_emplacement_programme_CPS = args.fichier_emplacement_programme_CPS,water_cavity=args.water_cavity,layerbeg_water=args.layerbeg_water,layerend_water=args.layerend_water,percent_solid=args.percent_solid,percent_liquid=args.percent_liquid)
    
if __name__ == '__main__':
   
    args =  read_args()
    main(args)
