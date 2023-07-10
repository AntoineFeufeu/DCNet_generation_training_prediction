import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from matplotlib import colors
from matplotlib.pyplot import figure
from scipy import interpolate
import os
from numpy import loadtxt,savetxt
import argparse


def read_args():
    parser = argparse.ArgumentParser()
                        
    parser.add_argument("--nb_modes",
                        default=6,
                        type=int,
                        help="number of modes or cuves")
    parser.add_argument("--fmin",
                        default=20,
                        type=float,
                        help="minimal frequency")
    parser.add_argument("--fmax",
                        default=270,
                        type=float,
                        help="maximal frequency")
    parser.add_argument("--vpmin",
                        default=10,
                        type=float,
                        help="minimal velocity")                  
    parser.add_argument("--vpmax",
                         default=1800,
                        type=float,
                        help="maximal velocity")
    parser.add_argument("--image_dim_x",
                        default=256,
                        type=int,
                        help="dimension of the image on X axis")
    parser.add_argument("--image_dim_y",
                        default=512,
                        type=int,
                        help="dimension of the image on Y axis")
    parser.add_argument("--stockage_data_name",
                        type=str,
                        help="Path to the directory of outputs")

    args = parser.parse_args()
    return args

# fonction de la variation de l'intervalle des fréquences [f-var(f),f+var(f)]
def variation(f,image_dim_x=256) :
  
  if f > 1e-10 :
    var = round((1.5 * math.exp(5))/(f+20)**1.5 + math.log(f)/3)
  else :
      var = 0

  if f > image_dim_x -3 :
    var =1
  
  return var

#fonction responsable de la générations des fichiers segmentation et embedding à partir des cps outputs
def generate_seg_emb(stockage_data_name,nb_modes,fmin,fmax,vpmin,vpmax,image_dim_x,image_dim_y) :
    
# num_cas : nombre des cas (modèles specfem)
#fmax: frequence maximale,fmin : fréquence mininmale,vpmin : vitesse de phase minimale,vpmax : vitesse de phase max
#nb_modes: nombres des courbes à utiliser
# image_dim_x : dimension de l'image selon l'axe X
# image_dim_y : dimension de l'image selon l'axe Y


    # chemin de fichier SDISPR.TXT dans le dossier des outputs destinés 
    sdispr = stockage_data_name + 'SDISPR.TXT'

    f = open(sdispr,'r+') 
    liste = f.readlines()
    lines=[]
    for x in liste :
        if x!= '\n' and 'FREQ' not in x :   # nettoyages des lignes vides et titres
            lines.append(x)
    freq =[]
    vp =[]
    i =-1
    for x, line in enumerate( lines) :

        if '#' in line :           #séparer les différents modes existants dans le fichier
            i = i +1
            freq.append([])
            vp.append([])
            continue
        
        if not '*' in line :
            f = float(line.split()[1])
            v = 1000*float(line.split()[2])
            if f > fmin and f < fmax and v > vpmin and v < vpmax  :
                freq[i].append(f)        #liste des fréquences(2ème colonne)
                vp[i].append(v)     #liste des vitesses de phase(1ème colonne)

    # scaling frequences et vitesses entre 0 et 255 et entre 0 et 511 pour dessiner les courbes dans une image512x256
    frequencies_rescaled =[]
    vitesses_rescaled = []
    for i in range(len(freq)) :

        freq_scaled = [ int ((image_dim_x-1) * ((f-fmin )/(fmax-fmin))) for f in  freq[i]]

        frequencies_rescaled.append(freq_scaled)

        vp_scaled = [ int ((image_dim_y-1) * ((v-vpmin )/(vpmax-vpmin))) for v in  vp[i]]

        vitesses_rescaled.append(vp_scaled)
    #linking tout les points pour avoir des courbes lisses
    new = []
    for s in range(nb_modes):
        f= interpolate.interp1d(list(map(float,frequencies_rescaled[s])), vitesses_rescaled[s])
        f_new = np.arange(min(frequencies_rescaled[s]),max(frequencies_rescaled[s]),step=0.001)
        new.append([f_new,f(f_new)])
        

    # ajout de l'envéloppe autour des courbes
    l1 = [] #nouvelle liste des fréquences 
    l2 = [] #nouvelle liste des vitesses
    
    #1ère courbe spéciale car elle doit etre fine au début
    e = new[0]
    l = int(len(e[0])/3)
    for f1,v1 in zip(e[0][1000:l],e[1][1000:l]) :
        
        v =np.arange(f1 -variation(f1,image_dim_x), f1 +variation(f1,image_dim_x) + 1)
        l1.append([])
        l2.append([])
        l1[0] += list(v)
        l2[0] += [v1] * len(v)
    for f1,v1 in zip(e[0][l:],e[1][l:]) :
        
        v =np.arange(f1 -variation(f1,image_dim_x), f1 +variation(f1,image_dim_x) + 1)
        l1.append([])
        l2.append([])
        l1[0] += list(v)
        if (variation(f1,image_dim_x) == 3) :
            l2[0] += [v1-1,v1-1,v1,v1,v1,v1+1,v1+1]
        elif (variation(f1,image_dim_x) == 2) :
            l2[0] += [v1-1,v1,v1,v1,v1+1]
        elif (variation(f1,image_dim_x) == 1) :
            l2[0] += [v1-1,v1,v1+1]
        else :
            l2[0] += [v1] * len(v) 
    #ajout les vitesses les plus proches dans l'intervalle des fréquences 
    for i in range(1,len(new)) :
        e = new[i]
        for f1,v1 in zip(e[0][:],e[1][:]) :

            v =np.arange(f1 -variation(f1,image_dim_x), f1 +variation(f1,image_dim_x) + 1)
            l1.append([])
            l2.append([])
            l1[i] += list(v)
            if (variation(f1,image_dim_x) == 3) :
                l2[i] += [v1-1,v1-1,v1,v1,v1,v1+1,v1+1]
            elif (variation(f1,image_dim_x) == 2) :
                l2[i] += [v1-1,v1,v1,v1,v1+1]
            elif (variation(f1,image_dim_x) == 1) :
                l2[i] += [v1,v1,v1+1]
            else :
                l2[i] += [v1] * len(v)
    #image de segmentation
    seg = np.zeros((image_dim_y,image_dim_x))
    for i in range(len(l1)):
        for x,y in zip(l1[i],l2[i]) :
            seg[image_dim_y-1-int(y)][int(x)] = 1

    #image de l'embedding
    embedding = np.zeros((image_dim_y,image_dim_x))
    j =1
    for i in range(len(l1)):
        for x,y in zip(l1[i],l2[i]) :
            embedding[image_dim_y-1 -int(y)][int(x)] = j
        j+=1
    plt.imshow(seg)
    plt.imsave(stockage_data_name+'segmentation.jpg',seg,cmap='gray')
    plt.clf()
    plt.imshow(embedding)
    plt.imsave(stockage_data_name+'embedding.jpg',embedding)
    plt.clf()
    savetxt(stockage_data_name+'segmentation.csv',seg,delimiter=',')
    savetxt(stockage_data_name+'embedding.csv',embedding,delimiter=',')




if __name__ == '__main__':
    
    args =  read_args()
    generate_seg_emb(stockage_data_name=args.stockage_data_name,nb_modes=args.nb_modes,
                    fmin=args.fmin,fmax=args.fmax,vpmin=args.vpmin,vpmax=args.vpmax,
                    image_dim_x = args.image_dim_x,image_dim_y=args.image_dim_y)
