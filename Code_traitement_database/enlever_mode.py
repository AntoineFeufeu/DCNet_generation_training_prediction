import os
import shutil
## This code aims to make space and errase some mode folders in a data base
stockage_data_name = "./ML_DATA/"
nb_mode_a_enlever = [6]

for n in nb_mode_a_enlever :
    print("numero_mode : ", n)
    for i in range (1,6) :
        
        if os.path.isdir(stockage_data_name+str(i).zfill(4)+"/mode_"+str(n)+"/") :
            print(i)
            shutil.rmtree(stockage_data_name+str(i).zfill(4)+"/mode_"+str(n)+"/")
