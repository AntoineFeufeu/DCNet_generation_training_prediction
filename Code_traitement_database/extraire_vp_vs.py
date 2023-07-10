import numpy as np
import shutil
import os
import argparse
import matplotlib.pyplot as plt
import random as rd

## This code makes an illustration on how vp and vs are in our soil
def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stockage_ML_data",
                        default="./ML_DATA/")
    parser.add_argument("--ideb",
                        type=int,
                        default="10001")
    parser.add_argument("--ifin",
                        type=int,
                        default="10002")
    args = parser.parse_args()
    return args
  
def creation_rho_vp_vs_nxmin_nzmin_nzmax_profondeurinterface (nb_layers,thickness,nb_mesh_per_layer,layer_middle,vp_middle) :
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

def main (args) :
    for i in range (args.ideb,args.ifin) :
        
        #Algorithm which read the log file and create a dictionnary with the parameters
        with open(args.stockage_ML_data+str(i).zfill(4)+"/log", 'r') as file:
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

        stockage_data_name = args.stockage_ML_data+str(i).zfill(4)+"/"
        
        thickness = float(variable_dict['thickness'])
        layer_middle = int (variable_dict['layer_middle'])
        vp_middle = int (variable_dict['vp_middle'])
        water_cavity = ( (variable_dict['water_cavity']) == 'True')
        layerbeg_water = int (variable_dict['layerbeg_water'])
        layerend_water = int( variable_dict['layerend_water'])
        percent_solid = float (variable_dict['percent_solid'])
        percent_liquid = float (variable_dict['percent_liquid'])
        #nb_layers = int (variable_dict['nb_layers'])
        creation_rho_vp_vs_nxmin_nzmin_nzmax_profondeurinterface (nb_layers=28,thickness=thickness,nb_mesh_per_layer=1,layer_middle=layer_middle,vp_middle=vp_middle)
        
        if water_cavity :
            percent_gaz = 1 - percent_solid - percent_solid
            rho_sol = rho[int((layerend_water+layerbeg_water)/2)]
            vp_sol = vp[int((layerend_water+layerbeg_water)/2)]
            vs_sol = vs[int((layerend_water+layerbeg_water)/2)]
            rho_karst = rho_sol*percent_solid+1000*percent_liquid+1.292*percent_gaz
            vp_karst = (percent_solid*vp_sol*rho_sol+percent_liquid*1500*1000+percent_gaz*340*1.292)/rho_karst
            vs_karst = percent_solid*vs_sol*rho_sol/rho_karst
            for i in range(layerbeg_water,layerend_water+1) :
                vp[i] = vp_karst
                vs[i] = vs_karst
                
        x = range(len(vp),0,-1)
        plt.clf()
        plt.plot(profondeur_interface[1:],vp,'x',label='vp')
        plt.plot(profondeur_interface[1:],vs,'x',label='vs')
        plt.xlabel("Profondeur (m)")
        plt.ylabel("Vitesses (m/s)")
        plt.legend()
        plt.savefig(stockage_data_name+"courbes_evolution_vp_vs.jpg")
        #plt.show()
       

if __name__ == '__main__':
    args =  read_args()
    main(args)
