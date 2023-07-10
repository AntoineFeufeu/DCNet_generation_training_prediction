import os
## This code aims to change numbers in a data base, with a simple translation
stockage_data_name = "D:/STAGE/ML-DATA-WATER_repo/"
decalage = +3380
for i in range (2601,5481) :
    print (i)
    if os.path.isdir(stockage_data_name+str(i).zfill(4)) and not os.path.isdir(stockage_data_name+str(i+decalage).zfill(4)) :
        os.rename(stockage_data_name+str(i).zfill(4),stockage_data_name+str(i+decalage).zfill(4))
    elif not os.path.isdir(stockage_data_name+str(i).zfill(4)) :
        print (str(i) + " n'existe pas")
    else :
        print (str(i+decalage) + " existe déjà pour le remplacement de " + str(i)) 

