colonne1 = []
colonne2 = []
## A little code to sees the receivers output
import matplotlib.pyplot as plt

with open('./OUTPUT_DATA/0001/OUTPUT_FILES/AA.S0004.BXZ.semv', 'r') as f:
    for line in f:
        cols = line.split()
        colonne1.append(float(cols[0]))
        colonne2.append(float(cols[1]))
plt.plot(colonne1,colonne2)
plt.show()


