import numpy as np
import scipy.fftpack as sf
import shutil
import matplotlib.pyplot as plt
import pandas as pd
import cmath

import scipy.ndimage
from numpy import loadtxt,savetxt
import os
from PIL import Image, ImageEnhance
import argparse

def read_args():
    parser = argparse.ArgumentParser()
                        
    parser.add_argument("--nrec",
                        default=96,
                        type=int,
                        help="number of receivers")
    parser.add_argument("--dt",
                        default=3e-5,
                        type=float,
                        help="time step")
    parser.add_argument("--nstep_output_receiver",
                        default=10,
                        type=int,
                        help="Number of steps between two measures for the receivers")
    parser.add_argument("--dx",
                        default= 1,
                        type=float,
                        help=" distance between receivers")
    parser.add_argument("--fmin",
                        default=20,
                        type=float,
                        help="minimal frequency")
    parser.add_argument("--fmax",
                        default=240,
                        type=float,
                        help="maximal frequency")
    parser.add_argument("--vpmin",
                        default=100,
                        type=float,
                        help="minimal velocity")
    parser.add_argument("--vpmax",
                         default=2000,
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
    parser.add_argument("--first_rec_on_source",
                        default=False,
                        help="First receiver on the source or not")
    parser.add_argument("--stockage_data_name",
                        type=str,
                        default="./ML_DATA/500008/",
                        help="Path to the directory of outputs")
    args = parser.parse_args()
    return args




#fonction d'ajout d'une fonction gaussienne à la fin de la série temporeille 
def add_quee(y_matrix):

    new_list = np.array([y_matrix[-1]*np.exp(-0.0001  *(i - len(y_matrix))**2 ) for i in range(len(y_matrix) +1,len(y_matrix) +1001)])
    new = np.concatenate((y_matrix, new_list), axis = 0)
    return new

#la formule de calcul de dispersion 
def calcul_dispersion (f,vp,fft, x) :


    V = (fft/abs(fft)) * cmath.exp(2j* cmath.pi * f * x/vp)
    

    return V

#fonction qui lit tout les fichiers .semv de chaque récepteur et appique la fft ensuite
def FFT(nrec,specfem_outputs):
    liste_fft =[]
    for i in range(1,nrec+1) :
        path = specfem_outputs + 'AA.S'+str(i).zfill(4)+'.BXZ.semv' ######
        with open(path,'r') as f :
            liste = f.readlines()
            y_t = np.array([float(x.split()[1]) for x in liste])
            y_t = add_quee(y_t)
            fft = sf.fft(y_t,len(y_t))[0:int(len(y_t))]
            liste_fft.append(fft)
    return liste_fft

def dispersion_matrix (nrec,specfem_outputs,dt,dx,fmin,fmax,vpmin,vpmax,image_dim_y,first_rec_on_source) :
    # vitesse de phase
    vitesses_phase=np.linspace(vpmin,vpmax,image_dim_y)
    liste_fft = FFT(nrec,specfem_outputs)
    frequences = sf.fftfreq(len(liste_fft[0]),dt)[0:int(len(liste_fft[0])/2)]
    f_min_indice = (np.abs(frequences - fmin)).argmin()
    f_max_indice = (np.abs(frequences - fmax)).argmin()
    frequences_borned = frequences[f_min_indice:f_max_indice+1]
    vp_len , f_len = (len(vitesses_phase),len(frequences_borned)) 
    disperion_matrix = np.zeros((vp_len,f_len), dtype=np.float64)
    
    if (first_rec_on_source):
        x = dx * np.arange(nrec).astype(float) 
    else :
        x = dx * np.arange(1,nrec+1).astype(float) 

    for v_id, vp in enumerate(vitesses_phase) :
        for id_f,f in enumerate(frequences_borned) :
            sum = 0
            for rec_id in range(len(liste_fft)) :
                fft = liste_fft[rec_id][id_f+f_min_indice]
                sum = sum + calcul_dispersion(f,vp,fft,x[rec_id])*dx
            disperion_matrix[vp_len -v_id -1 ][id_f] = abs(sum)
    disperion_matrix=(disperion_matrix/np.max(disperion_matrix))
    return disperion_matrix,frequences_borned

def visualize_dispersion (dispersion_matrix,fmin,fmax,vpmin,vpmax) : #fonction de visualisation des images de dispersion
    plt.figure(figsize=(10,10))
    plt.imshow(dispersion_matrix ,extent = [fmin,fmax,vpmin,vpmax] ,cmap ='jet') 
    plt.xlabel('fréquence(Hz)')
    plt.ylabel('vitesse_de_phase(m/s)')
    plt.locator_params(axis='y', nbins=40)
    plt.axis("auto")
    plt.show()

def generate_dispersion_image(matrix,stockage_data_name) : # générer une image RGB avec 3 channels 

  
  plt.imsave(stockage_data_name+'upsampled_image.jpg',matrix,cmap='jet')
  upsampled_image = Image.open(stockage_data_name+'upsampled_image.jpg')

  sharpness_enhancer = ImageEnhance.Sharpness(upsampled_image)
  sharpner_image = sharpness_enhancer.enhance(factor=1.1)
  color_enhancer = ImageEnhance.Color(sharpner_image)
  colored_image = color_enhancer.enhance(factor=1.1)
  contrast_enhancer = ImageEnhance.Contrast(colored_image)
  contrasted_image = contrast_enhancer.enhance(factor=1.1)
  contrasted_image.save(stockage_data_name+'dispersion.jpg')
  
 

  print('dispersion.jpg generated !')

def process_dispersion(matrix,image_dim_x): #transforme la taille de l'image en 512x256
  rate = image_dim_x/matrix.shape[1] # le coefficient de l'up_sampling sur l'axe X
  matrix = scipy.ndimage.zoom(matrix, (1,rate), order=1)
  return matrix

def generate_dispersion_matrix (matrix,stockage_data_name) : #générer le fichier csv de dispersion
    
    savetxt(stockage_data_name+'dispersion_matrix.csv',matrix,delimiter=',')
    print('dispersion_csv generated !')

def save_boundary_csv(fmin,fmax,vmin,vmax,stockage_data_name) : #généer le fichier boundary.csv
    fv=np.array([fmin,fmax,vmin,vmax])
    savetxt(stockage_data_name+'fv_boundary.csv',fv,delimiter=',')
    print('fv_boundary generated !')

def main(args):

    
    #### ICI MODIF 2DT
    output_dispersion= dispersion_matrix(args.nrec,specfem_outputs=args.stockage_data_name+"OUTPUT_FILES/",dt=(args.nstep_output_receiver*args.dt),dx=args.dx,fmin=args.fmin,fmax=args.fmax,vpmin=args.vpmin,vpmax = args.vpmax,image_dim_y= args.image_dim_y,first_rec_on_source=args.first_rec_on_source)

    frequencies = output_dispersion[1]

    matrix = output_dispersion[0]

    dispersion_matrix_256 = process_dispersion(matrix,args.image_dim_x)

    generate_dispersion_matrix(dispersion_matrix_256,stockage_data_name=args.stockage_data_name)

    generate_dispersion_image(dispersion_matrix_256,stockage_data_name=args.stockage_data_name)

    
    save_boundary_csv(frequencies[0],frequencies[-1],args.vpmin,args.vpmax,stockage_data_name=args.stockage_data_name)

if __name__ == '__main__':
   
    args =  read_args()
    main(args)









