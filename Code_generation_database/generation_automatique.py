import subprocess
import os
import random as rd

## Here is a simple code to run many generation of samples in order to create the dataset

# LAUNCH AUTOMATISATION_GENERATION.PY
t = [(1+i*0.01) for i in range(1,4)]
v = [(3100+140*i) for i in range(1,5)]
l = [i for i in range (10,12)]
f= [[30,240],[10,230]]
vp = [[10,1600],[100,2200]]
i = 10000

lay = [(15,20)]
for thickness in t :
    for vp_l in v :
        for layer in l :
            for freq in f :
                for velo in vp :
                    for l_w in lay :
                        for j in range(1) :
                            i+=1
                            # beg = rd.randint(1,18)
                            # end = beg + rd.randint(1,7)
                            # # nxbeg = rd.randint(1,50)
                            # # nxend = nxbeg + rd.randint(5,35)
                            # nxbeg = 1
                            # nxend = 90
                            print(i)
                        
                            arguments = ["--vpmin="+str(velo[0]),"--vpmax="+str(velo[1]),"--fmin="+str(freq[0]),"--fmax="+str(freq[1]),"--thickness="+str(thickness), "--vp_middle="+str(vp_l), "--layer_middle="+str(layer),"--name="+str(i).zfill(4)]
                            subprocess.run(["python", "Code_generation_database/automatisation_generation.py"] + arguments)
