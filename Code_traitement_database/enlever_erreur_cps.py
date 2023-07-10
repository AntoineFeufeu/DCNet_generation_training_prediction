import matplotlib.pyplot as plt
from numpy import loadtxt
import shutil

## Sometimes, itappears that there is a bug in the CPS resolution, with some set of parameters, and the modes are incoherent (it espacially makes some nodes with horizontal lines, as if it was blocked at a velocity)
## This codes finds them because it is always loked at the same values, and errase samples that have this bug
ideb = 1
ifin = 3
stockage_ML_data= "./ML_DATA/"
for j in range (ideb,ifin) :
    print(j)
    f = open(stockage_ML_data+str(j).zfill(4)+"/SDISPR.TXT",'r+') 
    liste = f.readlines()
    lines=[]
    for x in liste :
        if x!= '\n' and 'FREQ' not in x :   # nettoyages des lignes vides et titres
            lines.append(x)
    vp =[]
    i =-1
    
    for x, line in enumerate( lines) :

        if '#' in line :           #séparer les différents modes existants dans le fichier
            i = i +1
            vp.append([])
            continue
        vp[i].append(1000*float(line.split()[2]))     #liste des vitesses de phase(1ème colonne)
        if i >= 3 and len(vp[i])>4 and  vp[i][-1] == vp[i][-3] == vp[i][-5] and abs(vp[i][-1]-690)< 1e-6 :
            f.close()
            print(j,i)
            vf=loadtxt(stockage_ML_data+"%04d/embedding.csv"%j,delimiter=',')
            plt.imshow(vf)
            plt.show()
            #shutil.rmtree(stockage_ML_data+str(j).zfill(4)+"/")
            print("Enlever" + str(j))
            break
