import matplotlib.pyplot as plt
from numpy import loadtxt
import random as rd
from PIL import Image
import matplotlib.image as mpimg
import shutil
import os
# This code aims to sort good and bad predictions made by the model, in order to make semi-supervised learning. [y] : you keep it, [n] : you errase it
counter = 0
ideb = 1
ifin = 5
stockage_predictions = "./DCNet_output/"
stockage_ML_data = "./ML_DATA/"
for j in range (ideb,ifin) :
    
    if os.path.isdir(stockage_predictions+"%04d/"%j) :
        print(j)
        fig, ax = plt.subplots()
        if not os.path.exists(stockage_predictions+"%04d/embedding.csv"%j) :
            shutil.rmtree(stockage_predictions+"%04d/"%j)
            continue
        #vf=loadtxt("./DCnet_output/%04d/embedding.csv"%j,delimiter=',')
        img = mpimg.imread(stockage_predictions+"%04d/Fitted_Curves_Filtered.jpg"%j)
        ax.imshow(img)

        # Récupérer la taille de l'écran et calculer les coordonnées pour la position de la fenêtre
        screen_width, screen_height = plt.gcf().canvas.get_width_height()
        window_width, window_height = 640, 480  # Remplacer ces valeurs par les dimensions de la fenêtre que vous souhaitez
        window_x, window_y = 800, 100  # Modifier ces valeurs pour ajuster la position de la fenêtre

        # Définir la position de la fenêtre
        manager = plt.get_current_fig_manager()
        manager.window.setGeometry(window_x, window_y, window_width, window_height)

        while True:
            plt.show(block=False)
            user_input = input("Voulez-vous conserver le dossier correspondant (y/n) ? ")
            plt.close()
            if user_input.lower() in ("y", "n"):
                break

        if user_input.lower() == "y":
            counter+=1
            shutil.copyfile(stockage_ML_data+"%04d/fv_boundary.csv"%j,stockage_predictions+"%04d/fv_boundary.csv"%j)
            shutil.copyfile(stockage_ML_data+"%04d/log"%j,stockage_predictions+"%04d/log"%j)
            print("Le dossier " + str(j).zfill(4) + " sera conservé.")
        else:
            print("Le dossier " + str(j).zfill(4) + " sera supprimé.")
            shutil.rmtree(stockage_predictions+"%04d/"%j)
            # Supprimer le dossier
print(counter)
print(counter/(ifin-ideb))
