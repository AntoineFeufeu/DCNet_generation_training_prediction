import numpy as np 
#import pandas as pd
import os
#import matplotlib.pyplot as plt
from PIL import Image
from numpy import loadtxt
class DataGenerator():
    """
    Data generator for the DCNet dataset (ML_DATA). This class should be used when training our Keras multi-output model.
    """
        
    def generate_split_indexes(self,deb,length,TRAIN_TEST_SPLIT = 0.8):
        p = np.random.permutation(length)
        train_up_to = int(length * TRAIN_TEST_SPLIT)
        train_idx = p[:train_up_to] 
        test_idx = p[train_up_to:]
        train_up_to = int(train_up_to * TRAIN_TEST_SPLIT)
        train_idx, valid_idx = train_idx[:train_up_to], train_idx[train_up_to:] 
        train_idx=list(np.array(train_idx)+1+deb)
        valid_idx=list(np.array(valid_idx)+1+deb)
        test_idx=list(np.array(test_idx)+1+deb)
        print(train_idx, valid_idx, test_idx)
        return train_idx, valid_idx, test_idx
    
    def preprocess_image(self, img_path):
        """
        Used to perform some minor preprocessing on the image before inputting into the network.
        """
        #im = plt.imread(img_path)
        im = np.array(Image.open(img_path))
        im = np.array(im) / 255.0 
        return im
    def load_csv(self, path):
        """
        Used to load csv file.
        """
        data = loadtxt(path, delimiter=',')
        data=np.array(data)
        return data

        
    def generate_images(self, image_idx, is_training, batch_size,dataset_path):
        """
        Used to generate a batch with images when training/testing/validating our Keras model.
        """
        
        # arrays to store our batched data
        images, segmentations, embeddings = [], [], []
        while True:
            for idx in image_idx:
                segmentation_file = dataset_path+"/%04d/segmentation.csv"%idx
                embedding_file = dataset_path+"/%04d/embedding.csv"%idx
                image_file = dataset_path+"/%04d/dispersion.jpg"%idx
                im = self.preprocess_image(image_file)
                segmentation=self.load_csv(segmentation_file)
                embedding=self.load_csv(embedding_file)
                segmentations.append(segmentation)
                embeddings.append(embedding)
                images.append(im)
                
                # yielding condition
                if len(images) >= batch_size:
                    yield np.array(images), [np.array(segmentations), np.array(embeddings)]
                    images, segmentations, embeddings = [], [], []
                    
            if not is_training:
                break
