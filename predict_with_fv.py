import tensorflow as tf
from tensorflow.keras.models import Sequential,Model
from tensorflow.keras.layers import Dense,Flatten,Conv2D,MaxPooling2D,Dropout,BatchNormalization,\
Lambda,Input,Input,Conv2DTranspose,Add
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from tensorflow.keras.optimizers import Adam
import numpy as np 
import os
import matplotlib.pyplot as plt
from numpy import loadtxt,savetxt
from matplotlib import colors
from matplotlib.pyplot import figure
from structure import DCNet
import argparse
import logging
from sklearn.cluster import DBSCAN
import collections
from scipy.ndimage import gaussian_filter1d
#from Read_Dispersion import read_dispersion

def preprocess_image(img_path):
    im = plt.imread(img_path)
    im = np.array(im) / 255.0 
    return im
def load_csv(path):
    data = loadtxt(path, delimiter=',')
    data=np.array(data)
    return data
def output(model,idx,dataset_path): 
    image_file = dataset_path+"/%04d/dispersion.jpg"%idx
    im=preprocess_image(image_file)
    ims=im[np.newaxis,...]
    seg,emb=model.predict(ims)
    seg=seg[-1,...,-1]
    emb=emb[-1,...]
    seg=(seg>0.5).astype(int)
    im=im/np.max(im)*255
    im=im.astype(int)
    dispersion_matrix=load_csv(dataset_path+"/%04d/dispersion_matrix.csv"%idx)
    fv=loadtxt(dataset_path+"/%04d/fv_boundary.csv"%idx, delimiter=',')
    return im,seg,emb,dispersion_matrix,fv
def curve(final_label,dis_mat,k):
    curve=[]
    for i in range(final_label.shape[1]):
        velo=np.where(final_label[:,i]==k)
        if len(velo[0])==0:
            continue
        j=velo[0][np.argmax(dis_mat[velo,i])]
        curve+=[[511-j,i]]
    return np.array(curve)
def super_pose(list_fv,number):  ### this function plot predicted curves superimposed on initial input dispersion image.
    vf=loadtxt("D:\poly\internship\ML_DATA\%04d\\fv_boundary.csv"%number,delimiter=',')
    f=vf[:2]
    v=vf[2:]
    #fig, ax = plt.subplots(figsize=(6,6))
    figure(7,figsize=(6, 6), dpi=80)
    dispersion=loadtxt('D:\poly\internship\ML_DATA\%04d\dispersion_matrix.csv'%number,delimiter=',')
    #art=ax.imshow(dispersion, extent=[f[0],f[-1],v[0],v[-1]],cmap='jet')
    plt.imshow(dispersion, extent=[f[0],f[-1],v[0],v[-1]],cmap='jet')
    plt.axis("auto")
    for i in range(len(list_fv)):
        mode=list_fv[i]
        plt.plot(mode[:,-1],mode[:,-2],color='white')
    plt.title("Predicted curves superimposed on input dispersion image")
    plt.savefig('./DCNet_output/%04d/superimposed.jpg'%number)
    plt.show()
def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample",
                        default=1,
                        type=int,
                        help="Choosing the sample number and get from default file")
    parser.add_argument("--weight_f",
                        default=0.5,
                        type=float,
                        help="Choosing the weight for f")
    parser.add_argument("--weight_v",
                        default=0.5,
                        type=float,
                        help="Choosing the weight for v")
    parser.add_argument("--weight_dir",
                        default='./checkpoints_best_only/checkpoint',
                        #type=str,
                        help="Choosing the weight directory for model; by default: ./checkpoints_best_only/checkpoint")
    args = parser.parse_args()
    return args
def main(args):
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
    ###load model
    model=DCNet().build(512,256)
    model.load_weights(args.weight_dir)
    number=args.sample
    weight_f=args.weight_f
    weight_v=args.weight_v
    number=args.sample
    im,seg,emb,dis_mat,fv=output(model,number,"./ML_DATA/")
    if not os.path.exists('./DCNet_output/%04d'%number):
        os.makedirs('./DCNet_output/%04d'%number)
    cmap1=colors.ListedColormap(['black','white'])
    norm1=colors.Normalize(vmin=0,vmax=1)
    ### show input dispersion image
    figure(0,figsize=(8, 8), dpi=80)
    plt.imshow(im,cmap='jet')
    plt.axis("auto")
    plt.title("Input dispersion image")
    plt.show()
    ### show output of segmentation branch
    figure(1,figsize=(8, 8), dpi=80)
    plt.imshow(seg,cmap=cmap1,norm=norm1)
    plt.axis("auto")
    plt.title("Segmentation")
    plt.savefig('./DCNet_output/%04d/Segmentation.jpg'%number)
    plt.show()
    ### show output of embedding branch
    figure(2,figsize=(8, 8), dpi=80)
    emb1=emb/np.max(emb)*255
    emb1=emb1.astype(int)
    plt.imshow(emb1,cmap='jet')
    plt.axis("auto")
    plt.title("Embedding")
    plt.savefig('./DCNet_output/%04d/Embedding.jpg'%number)
    plt.show()
    ### DBSCAN
    print("Executing DBSCAN...")
    seg_plat=np.reshape(seg,(-1,))
    energy=np.where(seg_plat==1)
    emb_plat=[]
    for i in range(emb.shape[0]):
        for j in range(emb.shape[1]):
            temp=np.concatenate((emb[i,j,:],np.array([np.max(emb)*weight_v*i/emb.shape[0],np.max(emb)*weight_f*j/emb.shape[1]])))
            emb_plat+=[temp]
    emb_plat=np.array(emb_plat)
    emb_concerned=emb_plat[energy]
    clustering = DBSCAN(eps=0.5, min_samples=50).fit(emb_concerned)
    clust=clustering.labels_
    clust=clust-min(clust)+1
    final_label=np.repeat(0,int(emb.shape[0]*emb.shape[1]))
    final_label[energy]=clust
    count=collections.Counter(final_label.reshape(-1,))
    threshold=1000#800
    for i in range(1,len(count)):
        if count[i]<threshold:
            final_label[np.where(final_label==i)]=0
    unique = np.unique(final_label, return_counts=False)
    for i in range(len(unique)):
        final_label[np.where(final_label==unique[i])]=i
    final_label=final_label.reshape(seg.shape)
    #cmap2 = colors.ListedColormap(["black",'grey','blue','green','violet','yellow','red','brown'])
    #norm2=colors.Normalize(vmin=0,vmax=7)
    ### show different modes by DBSCAN
    figure(3,figsize=(6, 6), dpi=80)
    plt.imshow(final_label,cmap='inferno')#cmap=cmap2,norm=norm2)
    plt.axis("auto")
    plt.title("Different modes by DBSCAN")
    plt.savefig('./DCNet_output/%04d/groups.jpg'%number)
    plt.show()
    #### Show fitted curved found from maximum peak.
    figure(4,figsize=(6, 6), dpi=80)
    plt.imshow(final_label,extent=[0,255,0,511],cmap='inferno')#cmap=cmap2,norm=norm2),cmap=cmap2,norm=norm2)
    plt.axis("auto")
    for i in range(1,len(unique)):
        mode=curve(final_label,dis_mat,i)
        plt.plot(mode[:,1],mode[:,0],color='white')
    plt.title('Fitted Curves')
    plt.savefig('./DCNet_output/%04d/Fitted_Curves.jpg'%number)
    plt.show()
    ##### Show rectified curves after using gaussian filter.
    figure(5,figsize=(6, 6), dpi=80)
    plt.imshow(final_label,extent=[0,255,0,511],cmap='inferno')#cmap=cmap2,norm=norm2)cmap=cmap2,norm=norm2)
    plt.axis("auto")
    for i in range(1,len(unique)):
        mode=curve(final_label,dis_mat,i)
        plt.plot(mode[:,1],gaussian_filter1d(mode[:,0],sigma=5),color='white')
    plt.title('Fitted Curves After Filtering')
    plt.savefig('./DCNet_output/%04d/Fitted_Curves_Filtered.jpg'%number)
    plt.show()
    list_fv=[]
    #### show the curves on real scale after rescaling.
    figure(6,figsize=(6, 6), dpi=80)
    for i in range(1,len(unique)):
        mode=curve(final_label,dis_mat,i)
        mode[:,0]=gaussian_filter1d(mode[:,0],sigma=5)
        mode[:,0]=mode[:,0]/511*(fv[-1]-fv[-2])+fv[-2]
        mode[:,1]=mode[:,1]/255*(fv[1]-fv[0])+fv[0]
        fv_temp=np.zeros((int(len(mode[:,0])),3))
        fv_temp[:,0]=i
        fv_temp[:,1:]=mode
        list_fv+=[fv_temp]
        plt.plot(mode[:,1],mode[:,0])
    plt.grid()
    plt.xlabel('Frequency(Hz)')
    plt.ylabel("Phase velocity(m/s)")
    plt.axis("auto")
    plt.title("Dispersion Curves")
    plt.savefig('./DCNet_output/%04d/Predicted_Fitted_Curves.jpg'%number)
    plt.show()
    #### show predicted curves superimposed on initial input dispersion image with function super_pose()
    super_pose(list_fv,number)
    for i in range(1,len(list_fv)):
        list_fv[0]=np.concatenate((list_fv[0], list_fv[i]))
    if not os.path.exists('./DCNet_output/%04d'%number):
        os.makedirs('./DCNet_output/%04d'%number)
    #### save curves under form of .csv first column: mode number, 2nd column: phase velocity, 3rd column: frequency.
    savetxt('./DCNet_output/%04d/predicted_fv_curves.csv'%number, list_fv[0], delimiter=',')
    return
if __name__ == '__main__':
  args = read_args()
  main(args)


