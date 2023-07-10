
import matplotlib.pyplot as plt
from PIL import Image, ImageEnhance
import argparse
from numpy import loadtxt
## This code allows to recreate in a sample the good dispersion image with the dispersion matrix
stockage_data = "./ML_DATA/"


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
  


for i in range(9998,9999) :
    print(i)
    stockage_data_name = stockage_data + str(i).zfill(4) + "/"
    matrix = loadtxt(stockage_data_name+"dispersion_matrix.csv",delimiter=',')
    generate_dispersion_image(matrix,stockage_data_name)
