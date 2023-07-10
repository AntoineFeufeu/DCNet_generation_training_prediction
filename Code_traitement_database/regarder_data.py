import matplotlib.pyplot as plt
from numpy import loadtxt
import random as rd
# A little code to look at some samples you want to check(number of modes, shape,...)
nb_sample=10
ideb_rand = 1
ifin_rand = 20
stockage_data_name = "./ML_DATA/"
for j in range (nb_sample) :
    i = rd.randint(ideb_rand,ifin_rand)
    #i=j
    print(i)
    vf=loadtxt(stockage_data_name+"%04d/dispersion_matrix.csv"%i,delimiter=',')
    plt.imshow(vf,cmap='jet')
    plt.show()
    vf=loadtxt(stockage_data_name+"%04d/embedding.csv"%i,delimiter=',')
    plt.imshow(vf,cmap='inferno')
    plt.show()

# for j in range (1,10) :
#     i = rd.randint(1,2842)
#     # i=j
#     print(i)

#     # Ouvrir l'image
#     image = Image.open("./ML_DATA/%04d/superposition.jpg"%i)

#     # Afficher l'image
#     image.show()