In all the python programs that I implemented, I import:
```
from numpy import loadtxt,savetxt
```
to save and load the csv files.

### Saving and reading dispersion in compatible .csv format - How to save dispersion matrix as 3-channel images with cmap='jet' in Python 

Let's say you have a dispersion matrix in numpy array format. In our case, the shape of the matrix is 512x256. To save this in the csv format compatible with the programs the I implemented, you can save it as follows:
```
savetxt('path/dispersion_matrix.csv',matrix,delimiter=',')
```
Where: 'path/dispersion_matrix.csv' is the path and name of file that we choose. matrix is the dispersion matrix in numpy array format.

Then, if you want to load this to use in your programs:
```
matrix=loadtxt('path/dispersion_matrix.csv',delimiter=',')
```
with matrix being the name of the variable that you can choose at your convenience.
So, we have matrix being the dispersion matrix read from 'path/name_of_file.csv' and you can work normally with it in python.

Next, to save dispersion matrix as 3-channel images with cmap='jet' in Python:
```
import matplotlib.pyplot as plt
plt.imsave('path/ dispersion.jpg',matrix,cmap='jet')
```
### Saving and reading embedding.csv and segmentation.csv in compatible format
In fact, they are just 2D-matrices. So, to save and read these files in compatible format, we do the same things as for dispersion matrix.
### Saving and reading fv_boundary.csv in compatible format
To save file:
```
import numpy as np
fmin=0.
fmax=270.
vmin=10.
vmax=1500.
fv=np.array([fmin,fmax,vmin,vmax])
savetxt('path/fv_boundary.csv',fv,delimiter=',')
```
To load file:
```
fv=loadtxt('path/fv_boundary.csv',delimiter=',')
```
### Reading and presenting the evolution of loss function from file loss_train_val.csv created by train.py

Let's say we have run the program train.py. At the end of training, we have a file called loss_train_val.csv which contains the evolution of loss function on training and validating set. To read this file from python: 

```
import matplotlib.pyplot as plt
loss= loadtxt('some_path/loss_train_val.csv',delimiter=',')
train_loss=loss[0,:]
valid_loss=loss[1,:]
```
By doing this, we have 2 variables train_loss and valid_loss containing evolution of loss function on training and validating set, respectively. So, you can work normally with these variables in python. For example, if we want to plot train_loss and valid_loss:
```
import matplotlib.pyplot as plt
plt.plot(train_loss, label='training')
plt.plot(valid_loss, label='validating')
plt.legend()
plt.show()
```
### Reading and presenting the predicted dispersion curves created by predict.py

Let's say we have predicted for sample=500 by using the program predict.py. To read and present them in Python:
```
import matplotlib.pyplot as plt
for number in [500]:
    curves=loadtxt('./DCNet_output/%04d/predicted_fv_curves.csv'%number, delimiter=',')
    print(curves.shape)
    u,_ = np.unique(curves[:,0], return_inverse=True)
    print(u)
    for i in u:
        c=curves[np.where(curves[:,0]==i)][:,1:]
        plt.plot(c[:,1],c[:,0])
    plt.grid()
    plt.show()
```


















