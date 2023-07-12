
import shutil
import os
import argparse
### . ~/.bashrc
##################################### INITIALISATION ####################################################################


def read_args():
    
    parser = argparse.ArgumentParser()
                        
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
    parser.add_argument("--vp_middle",
                        default=3100.,
                        type = float,
                        help="vp of the break (3100->3800)")
    parser.add_argument("--layer_middle",
                        default=11,
                        type = int,
                        help="layer of the break (11->15)")
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
    parser.add_argument("--fichier_parfile_initial",
                        default=r"C:\Users\feufe\OneDrive\Bureau\Stage\Code\STAGE-GET\Test_src\Par_file",
                        help="Place where the initial pare_file (the one in Test_src) is stocked")
    parser.add_argument("--fichier_interface_initial",
                        default=r"C:\Users\feufe\OneDrive\Bureau\Stage\Code\STAGE-GET\Test_src\interfaces_simple_topo_flat.dat",
                        help="Place where the initial interfaces (the one in Test_src) is stocked")
    parser.add_argument("--fichier_source_initial",
                        default=r"C:\Users\feufe\OneDrive\Bureau\Stage\Code\STAGE-GET\Test_src\SOURCE",
                        help="Place where the initial source (the one in Test_src) is stocked")
    parser.add_argument("--chemin_stockage",
                        default=r"C:\Users\feufe\OneDrive\Bureau\Test_dest",
                        help="Folder where you want to stock the files generated")
    args = parser.parse_args()
    return args

args = read_args()

nb_layers = min(args.nb_layers,28)
thickness = args.thickness
nb_mesh_per_layer = args.nb_mesh_per_layer
xdeb = args.xdeb
xfin = args.xfin
vp_middle=args.vp_middle
layer_middle=args.layer_middle
fichier_parfile_initial = args.fichier_parfile_initial
fichier_interface_initial = args.fichier_interface_initial
fichier_source_initial = args.fichier_source_initial
chemin_stockage =  args.chemin_stockage
f0 = args.f0
burst_band_width = args.burst_band_width
nxmax= args.nxmax

material_number=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28]
rho= [2449.0, 2456., 2473., 2488., 2487., 2475., 2464., 2457., 2451., 2445., 2435., 2424., 2411., 2393., 2377., 2361., 2343., 2322, 2295., 2253., 2191., 2085., 1958., 1865., 1764., 1659., 1580., 1558.]

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
    if ((i+1)<layer_middle) :
        vp.append(4000-((4000-vp_middle)/(layer_middle-1))*i)
    elif ((i+1) == layer_middle) :
        vp.append(vp_middle)
    else :
        vp.append(vp_middle-((vp_middle-600)/(nb_layers-layer_middle))*(i-layer_middle+1))
    vs.append(vp[-1]/(3**(1/2)))
profondeur_interface = [-epaisseur]
for i in range(nb_layers) :
    profondeur_interface.append(profondeur_interface[-1]+thickness**(nb_layers-i))

#####################################  FIN INITIALISATION ###############################################################



with open("XDmod_28-layers.dep", "r+") as f:
    # Lire tout le contenu du fichier
    content = f.read()

     # Trouver l'emplacement du texte à modifier
    position = content.find("HR      VP      VS     RHO QP  QS  ETAP ETAS FREFP FREFS")
    if position == -1:
        print("Le texte à modifier n'a pas été trouvé.")
    else:
        # Déplacer le curseur de fichier à l'emplacement juste après le texte à modifier
        f.seek(position + len("HR      VP      VS     RHO QP  QS  ETAP ETAS FREFP FREFS"))
        for i in range (nb_layers-1,-1,-1) :
            f.write("\n"+str(thickness**(nb_layers-1-i)/1000)+" "+str(vp[i]/1000)+" "+str(vs[i]/1000)+" "+str(rho[i]/1000)+" "+" 0 0 0 0 1 1")
        f.write("\n0.00000000      1.20000      0.690000      2.00000       0.000000       0.0000000       0.00000000       1.00000000       1.00000000")
