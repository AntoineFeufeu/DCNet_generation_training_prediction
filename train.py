import tensorflow as tf
from tensorflow.keras.models import Sequential,Model
from tensorflow.keras.layers import Dense,Flatten,Conv2D,MaxPooling2D,Dropout,BatchNormalization,\
Lambda,Input,Input,Conv2DTranspose,Add
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from tensorflow.keras.optimizers import Adam
import numpy as np
from numpy import savetxt 
#import pandas as pd
import os
#import matplotlib.pyplot as plt
from numpy import loadtxt
from loss import discriminative_loss,segmentation_loss
from Feed_Data import DataGenerator
from structure import DCNet
import logging
import argparse
#from data_reader import Config, DataReader
from keras import backend as K

def dice_coef(y_true, y_pred,smooth=1):
    y_true_f = K.flatten(y_true)
    y_pred_f = K.flatten(y_pred)
    intersection = K.sum(y_true_f * y_pred_f)
    return (2. * intersection + smooth) / (K.sum(y_true_f) + K.sum(y_pred_f) + smooth)

def compile_model(model,init_lr = 0.0005, epochs = 15):
    opt = Adam(learning_rate=init_lr)#, decay=init_lr / epochs)
    model.compile(optimizer=opt, 
            loss={
                  'segmentation_output': segmentation_loss, 
                  'embedding_output': discriminative_loss},
            loss_weights={
                  'segmentation_output': 0.5, 
                  'embedding_output': 0.5},
            metrics={'segmentation_output': ['accuracy',dice_coef,tf.keras.metrics.MeanIoU(num_classes=2),tf.keras.metrics.Recall(),tf.keras.metrics.Precision()],
                      'embedding_output': 'accuracy'
                  }
                 )

def get_callbacks(path='./checkpoints_best_only/checkpoint'):
    checkpoint_path=path
    checkpoint=ModelCheckpoint(filepath=checkpoint_path,frequency="epoch",save_weights_only=True,
                               verbose=1,save_best_only=True,monitor='val_loss',mode='min')
    early=tf.keras.callbacks.EarlyStopping(monitor="loss",patience=5,mode='min',verbose=1)
    plateau= tf.keras.callbacks.ReduceLROnPlateau(monitor="loss",patience=3,factor=0.2)
    return [checkpoint,early,plateau]


def read_args():
    parser = argparse.ArgumentParser()
                        
    parser.add_argument("--epochs",
                        default=5,
                        type=int,
                        help="number of epochs")
    
    parser.add_argument("--batch_size",
                        default=8,
                        type=int,
                        help="batch size")
    parser.add_argument("--valid_batch_size",
                        default=4,
                        type=int,
                        help="valid batch size")
    
    parser.add_argument("--learning_rate",
                        default=0.0005,
                        type=float,
                        help="learning rate")
    parser.add_argument("--data_dir",
                        default="./ML_DATA",
                        help="Data set directory with format precised in the doc")

    # parameters for transfer learnings                    
    parser.add_argument("--weight_dir",
                        default= None,
                        help="model directory used for loading pretrained weights. If None: training from scratch\
                            (default: ./checkpoints_best_only/checkpoint)")  
    args = parser.parse_args()
    return args

def train_model(args):
    model=DCNet().build(512,256)
    compile_model(model,args.learning_rate,args.epochs)
    model.summary()
    if args.weight_dir is not None: 
        logging.info("loading pretrained to model for training ...")
        try: # restore the model
            model.load_weights(args.weight_dir)
        except:
            logging.info("please set a right path for weight_dir!")
            exit()
    callbacks=get_callbacks()
    batch_size = args.batch_size
    valid_batch_size = args.valid_batch_size
    data_generator = DataGenerator()
    # train_idx, valid_idx, test_idx = data_generator.generate_split_indexes(100,1-args.valid)
    train_idx=loadtxt('./train_idx.csv', delimiter=',')
    #test_idx=loadtxt('./test_idx.csv', delimiter=',')
    valid_idx=loadtxt('./valid_idx.csv', delimiter=',')
    train_gen = data_generator.generate_images(train_idx, is_training=True, batch_size=batch_size,dataset_path=args.data_dir)
    valid_gen = data_generator.generate_images(valid_idx, is_training=True, batch_size=valid_batch_size,dataset_path=args.data_dir)
    epochs = args.epochs
    history = model.fit(train_gen,
                    steps_per_epoch=len(train_idx)//batch_size,
                    epochs=epochs,
                    callbacks=callbacks,
                    validation_data=valid_gen,
                    validation_steps=len(valid_idx)//valid_batch_size)
    train_loss = history.history['loss']
    val_loss = history.history['val_loss']
    loss_train_val=np.zeros((2,int(epochs)))
    loss_train_val[0,:]=train_loss
    loss_train_val[1,:]=val_loss
    savetxt('./loss_train_val.csv', loss_train_val, delimiter=',')
    val_embedding_output_accuracy = history.history['val_embedding_output_accuracy']
    val_segmentation_output_mean_io_u = history.history['val_segmentation_output_mean_io_u']
    val_segmentation_output_dice_coef = history.history['val_segmentation_output_dice_coef']
    val_segmentation_output_accuracy = history.history['val_segmentation_output_accuracy']
    embedding_output_accuracy = history.history['embedding_output_accuracy']
    segmentation_output_mean_io_u = history.history['segmentation_output_mean_io_u']
    segmentation_output_dice_coef = history.history['segmentation_output_dice_coef']
    segmentation_output_accuracy = history.history['segmentation_output_accuracy']
    segmentation_output_recall = history.history['segmentation_output_recall']
    val_segmentation_output_recall = history.history['val_segmentation_output_recall']
    val_segmentation_output_precision = history.history['val_segmentation_output_precision']
    segmentation_output_precision = history.history['segmentation_output_precision']
    metrics =np.zeros((12,int(epochs)))
    metrics[0,:]= segmentation_output_accuracy
    metrics[1,:]= segmentation_output_dice_coef
    metrics[2,:]= segmentation_output_mean_io_u
    metrics[3,:]= embedding_output_accuracy
    metrics[4,:]= val_segmentation_output_accuracy
    metrics[5,:]= val_segmentation_output_dice_coef
    metrics[6,:]= val_segmentation_output_mean_io_u
    metrics[7,:]= val_embedding_output_accuracy
    metrics[8,:]= segmentation_output_recall
    metrics[9,:]= val_segmentation_output_recall
    metrics[10,:]= segmentation_output_precision
    metrics[11,:]= val_segmentation_output_precision
    savetxt('./metrics.csv', metrics, delimiter=',')
    return 0

def main(args):
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
    train_model(args)
    return

if __name__ == '__main__':
  args = read_args()
  main(args)




