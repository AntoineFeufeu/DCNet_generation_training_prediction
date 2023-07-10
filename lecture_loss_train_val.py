import csv
import matplotlib.pyplot as plt
import numpy as np

# Nom du fichier CSV
nom_fichier = "./loss_train_val.csv"

# Listes pour stocker les données
loss_train = []
loss_validation = []

# Lecture du fichier CSV
with open(nom_fichier, 'r') as fichier:
    lecteur_csv = csv.reader(fichier)
    lignes = list(lecteur_csv)
    
    # Stockage des données dans les listes
    loss_train = list(map(float, lignes[0]))
    loss_validation = list(map(float, lignes[1]))

# Création de l'axe des abscisses
x = np.arange(len(loss_train))

# Tracé des courbes
plt.plot(x, loss_train, label='loss_train')
plt.plot(x, loss_validation, label='loss_validation')

# Ajout des légendes et d'un titre
plt.xlabel('Epoques')
plt.ylabel('Loss')
plt.title('Courbes de loss_train et loss_validation')

# Ajout de la légende
plt.legend()

# Affichage du graphe
plt.show()








# import csv
# import copy
# import matplotlib.pyplot as plt
# import numpy as np
# from numpy import loadtxt,savetxt
# def apply_blur_filter(image):
#     # Créer une copie de l'image
#     blurred_image = copy.deepcopy(image)

#     # Obtenir les dimensions de l'image
#     num_rows = len(image)
#     num_cols = len(image[0])

#     # Parcourir chaque pixel de l'image (à l'exception des bords)
#     for i in range(1, num_rows - 1):
#         for j in range(1, num_cols - 1):
#             # Appliquer le filtre de convolution
#             blurred_value = (
#                 image[i - 1][j - 1] + image[i - 1][j] + image[i - 1][j + 1] +
#                 image[i][j - 1] + image[i][j] + image[i][j + 1] +
#                 image[i + 1][j - 1] + image[i + 1][j] + image[i + 1][j + 1]
#             ) / 9

#             # Mettre à jour la valeur du pixel dans l'image floue
#             blurred_image[i][j] = blurred_value

#     return blurred_image

# def load_csv_image(file_path):
#     image = []
#     with open(file_path, 'r') as file:
#         reader = csv.reader(file)
#         for row in reader:
#             # Convertir les valeurs de la ligne en entiers
#             row = [float(value) for value in row]
#             image.append(row)
#     return image

# # Chemin vers l'image CSV
# csv_image_path = './ML_DATA/9997/dispersion_matrix.csv'

# # Charger l'image CSV
# image = load_csv_image(csv_image_path)

# # Appliquer le filtre flou
# blurred_image = apply_blur_filter(image)

# # Convertir l'image en un tableau NumPy
# blurred_array = np.array(blurred_image)

# # Afficher l'image originale
# plt.subplot(1, 2, 1)
# plt.imshow(image, cmap='jet')
# plt.title('Image originale')
# plt.axis('off')

# # Afficher l'image floue
# plt.subplot(1, 2, 2)
# plt.imshow(blurred_array, cmap='jet')
# plt.title('Image floue')
# plt.axis('off')

# # Afficher les deux images côte à côte
# plt.tight_layout()
# plt.show()


# savetxt('./ML_DATA/9999/dispersion_matrix.csv',blurred_array,delimiter=',')