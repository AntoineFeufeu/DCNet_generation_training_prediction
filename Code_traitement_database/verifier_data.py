import os
import shutil
# A little code to see if sample folders are complete or not
stockage_data_name = "./ML_DATA/"
for i in range (1,5) :
    if not os.path.isdir(stockage_data_name+str(i).zfill(4)):
        print("Pas de dossier "+str(i).zfill(4))
        continue
    if not os.path.exists(stockage_data_name+str(i).zfill(4)+"/embedding.csv") :
        print("Pas de embedding "+str(i).zfill(4))
        #shutil.rmtree("nuwa/ML-DATA-NEW/"+str(i).zfill(4)+"/")
        #continue
    if not os.path.exists(stockage_data_name+str(i).zfill(4)+"/dispersion_matrix.csv") :
        print("Pas de dispersion "+str(i).zfill(4))
        #shutil.rmtree("nuwa/ML-DATA-NEW/"+str(i).zfill(4)+"/")
        #continue
    if not os.path.exists(stockage_data_name+str(i).zfill(4)+"/segmentation.csv") :
        print("Pas de segmentation "+str(i).zfill(4))
        #shutil.rmtree("nuwa/ML-DATA-NEW/"+str(i).zfill(4)+"/")
        #continue
