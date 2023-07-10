import subprocess
import logging
logging.getLogger('tensorflow').disabled = True
#l = [500,1500,2500,3000,3500,3700,4000,4500,1000010,1000011]
for i in range (1,2) :
    
    print(i)
    subprocess.run(["python", "predict_with_seg_clusters.py "]+ ["--sample="+str(i).zfill(4), "--nb_modes_voulu="+str(1),"--stockage_data="+"./ML_DATA/"])

