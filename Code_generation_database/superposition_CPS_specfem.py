import matplotlib.pyplot as plt
import numpy as np
import argparse
from matplotlib.pyplot import figure

def read_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--fmin",
                        default=20,
                        type=float,
                        help="minimal frequency")
    parser.add_argument("--fmax",
                        default=240,
                        type=float,
                        help="maximal frequency")
    parser.add_argument("--vpmin",
                        default=100,
                        type=float,
                        help="minimal velocity")
    parser.add_argument("--vpmax",
                         default=4000,
                        type=float,
                        help="maximal velocity")
    parser.add_argument("--image_dim_x",
                        default=256,
                        type=int,
                        help="dimension of the image on X axis")
    parser.add_argument("--image_dim_y",
                        default=512,
                        type=int,
                        help="dimension of the image on Y axis")
    parser.add_argument("--stockage_data_name",
                        type=str,
                        default="ML_DATA/3074/",
                        help="Path to the directory of outputs")

    args = parser.parse_args()
    return args


args = read_args()

image_dim_y = args.image_dim_y
image_dim_x = args.image_dim_x
fmin = args.fmin
fmax =  args.fmax
vpmin = args.vpmin
vpmax = args.vpmax
stockage_data_name = args.stockage_data_name
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
                mode_freqs[current_mode].append(int((f-fmin)*(image_dim_x)/(fmax-fmin)))
                mode_phase_vels[current_mode].append(int((v-vpmin)*image_dim_y/(vpmax-vpmin)))


data = np.loadtxt(stockage_data_name+"dispersion_matrix.csv", delimiter=',')
#print(len(data), len(data[0]))
for i in range (len(mode_freqs)) :
    #plt.plot(mode_freqs[i],mode_phase_vels[i])
    for j in range(len(mode_freqs[i])) :
        #print(mode_freqs[i][j],mode_phase_vels[i][j])
        data[image_dim_y-1-mode_phase_vels[i][j]][mode_freqs[i][j]] = 0.6
        #print(data[mode_freqs[i][j]][mode_phase_vels[i][j]])


x = (np.linspace(0, image_dim_x-1, 8)) #8
y = (np.linspace(0, image_dim_y-1, 12)) #12
str_x = [str(round((num/(image_dim_x-1))*(fmax-fmin)+fmin)) for num in x]
str_y = ([str(round((num/(image_dim_y-1))*(vpmax-vpmin)+vpmin)) for num in y])
str_y.reverse()
plt.xticks(x, str_x)
plt.yticks(y, str_y)
# Ajouter les étiquettes des axes
plt.xlabel("Fréquences")
plt.ylabel("Vitesses de phase")
plt.imsave(stockage_data_name+'superposition.jpg',data,cmap='jet')
#plt.show() 

plt.xticks(x, str_x)
plt.yticks(y, str_y)
plt.imshow(data, cmap='jet')
plt.colorbar()
plt.savefig(stockage_data_name+'superposition_labélisée.png')
# plt.show()
figure(2,figsize=(8, 8), dpi=80)
data2 = np.loadtxt(stockage_data_name+"dispersion_matrix.csv", delimiter=',')

plt.xticks(x, str_x)
plt.yticks(y, str_y)
plt.imshow(data2, cmap='jet')
plt.colorbar()
plt.savefig(stockage_data_name+'dispersion_labélisée.png')


