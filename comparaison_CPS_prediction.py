import matplotlib.pyplot as plt
import numpy as np
from numpy import loadtxt
from matplotlib.pyplot import figure
import matplotlib.pyplot as plt
## This code shows predicted curves
stockage_ML_data= "./ML_DATA/"
for i in range (1,3) :
    print(i)
    with open(stockage_ML_data+str(i).zfill(4)+"/log", 'r') as file:
            content = file.read()
            lines = content.split('\n')
            variable_dict = {}
            for line in lines:
                parts = line.split('=')
                if len(parts) != 2:
                    continue
                name = parts[0].strip()
                value = parts[1].strip()
                variable_dict[name] = value

            stockage_data_name = stockage_ML_data+str(i).zfill(4)+"/"

            
            fmin = float(variable_dict['fmin'])
            fmax = float(variable_dict['fmax'])
            vpmin = float(variable_dict['vpmin'])
            vpmax = float(variable_dict['vpmax'])
            image_dim_x = int(variable_dict['image_dim_x'])
            image_dim_y = int(variable_dict['image_dim_y'])

    with open(stockage_data_name+'SDISPR.TXT', 'r') as file:
        mode_freqs = []
        mode_phase_vels = []
        current_mode = -1
        for line in file:
            line = line.strip()
            if not line:
                # Ignore les lignes vides
                continue
            if 'RAYLEIGH' in line:
                
                current_mode = int(line.split()[-1])
                mode_freqs.append([])
                mode_phase_vels.append([])
            elif not '*' in line and not 'FREQ' in line and current_mode >= 0:
                freq, phase_vel = line.split()[1:]
                f = float(freq)
                v= 1000*float(phase_vel)
                if f > fmin and f < fmax and v > vpmin and v < vpmax  :
                    mode_freqs[current_mode].append(f)
                    mode_phase_vels[current_mode].append(v)



    
        for j in range (2) :
            x = np.array(np.array(mode_freqs[j]))
            y = np.array(mode_phase_vels[j])
            dy_dx = np.gradient(y, x)
            d2y_dx2 = np.gradient(dy_dx, x)
            d3y_dx3 = np.gradient(d2y_dx2, x)
            if j == 0 :
                plt.plot(x,y,color='black',label='CPS resolution')
            else :
                plt.plot(x,y,color='black')
            # if i == 70001 :
            #     plt.plot(x,y,color='black',label='Model without water')
            # if i == 70002 :
            #     plt.plot(x,y,color='black',label='Model with the layer full of water')
                # plt.plot(x,dy_dx,color='b',label='Profondeur max = '+str(i-3051+5)+'m',alpha=(i-3050)/50+0.5)
                # plt.plot(x,d2y_dx2,color='r',label='Profondeur max = '+str(i-3051+5)+'m',alpha=(i-3050)/50+0.5)
                # plt.plot(x,d3y_dx3,color='g')

    for number in [i]:
        curves=loadtxt('./DCNet_output/%04d/predicted_fv_curves.csv'%number, delimiter=',')
        print(curves.shape)
        u,_ = np.unique(curves[:,0], return_inverse=True)
        print(u)
        ii=0
        for k in u:
            c=curves[np.where(curves[:,0]==k)][:,1:]
            if ii ==0 :
                plt.plot(c[:,1],np.array(c[:,0]),color='r',label='Prediction')
            else :
                plt.plot(c[:,1],np.array(c[:,0]),color='r')
            ii+=1
    plt.title("Dispersion curves")
    plt.xlabel("Frequency")
    plt.ylabel("Phase velocity")
    plt.legend()
    plt.ylim(200,1800)
    
    plt.grid()
    #plt.savefig("C:/Users/feufe/Pictures/image à garder/passive/semi_plane/compare_"+str(int((i-60000)*10)))
    plt.show()
    plt.clf()
    