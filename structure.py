import tensorflow as tf
from tensorflow.keras.models import Sequential,Model
from tensorflow.keras.layers import Dense,Flatten,Conv2D,MaxPooling2D,Dropout,BatchNormalization,\
Lambda,Input,Conv2DTranspose,Add
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
class DCNet():
    """
    Used to generate our multi-output model. This CNN contains two branches, one for segmentation, other for 
    embedding. 
    """
    def make_encoder(self, inputs):

        x=Conv2D(16, (3, 3), strides=(1, 1), padding="same",activation="relu", kernel_initializer = 'he_normal')(inputs)
        x=BatchNormalization()(x)
        x=Conv2D(16, (3, 3), strides=(1, 1), padding="same",activation="relu", kernel_initializer = 'he_normal')(x)
        x1=BatchNormalization()(x)
        x=MaxPooling2D((2,2))(x1)
        x=Conv2D(32, (3, 3), strides=(1, 1), padding="same",activation="relu", kernel_initializer = 'he_normal')(x)
        x=BatchNormalization()(x)
        x=Conv2D(32, (3, 3), strides=(1, 1), padding="same",activation="relu", kernel_initializer = 'he_normal')(x)
        x2=BatchNormalization()(x)
        x=MaxPooling2D((2,2))(x2)
        x=Conv2D(64, (3, 3), strides=(1, 1), padding="same",activation="relu", kernel_initializer = 'he_normal')(x)
        x=BatchNormalization()(x)
        x=Conv2D(64, (3, 3), strides=(1, 1), padding="same",activation="relu", kernel_initializer = 'he_normal')(x)
        x=BatchNormalization()(x)
        x=Conv2D(64, (3, 3), strides=(1, 1), padding="same",activation="relu", kernel_initializer = 'he_normal')(x)
        x3=BatchNormalization()(x)
        x=MaxPooling2D((2,2))(x3)
        x=Conv2D(128, (3, 3), strides=(1, 1), padding="same",activation="relu", kernel_initializer = 'he_normal')(x)
        x=BatchNormalization()(x)
        x=Conv2D(128, (3, 3), strides=(1, 1), padding="same",activation="relu", kernel_initializer = 'he_normal')(x)
        x=BatchNormalization()(x)
        x=Conv2D(128, (3, 3), strides=(1, 1), padding="same",activation="relu", kernel_initializer = 'he_normal')(x)
        x4=BatchNormalization()(x)
        x=MaxPooling2D((2,2))(x4)
        return x,x1,x2,x3,x4
    def make_default_hidden_layers(self, outputs_4,x1,x2,x3,x4):
        """
        Used to generate a default set of hidden layers.
        """
        x = Conv2D(128, (3, 3), strides=(1, 1), padding="same",activation="relu", kernel_initializer = 'he_normal')(outputs_4)
        x=BatchNormalization()(x)
        x = Conv2D(128, (3, 3), strides=(1, 1), padding="same",activation="relu", kernel_initializer = 'he_normal')(x)
        x=BatchNormalization()(x)
        x = Conv2D(128, (3, 3), strides=(1, 1), padding="same")(x)
        x=BatchNormalization()(x)
        x=Conv2DTranspose(128, (4,4), strides=(2, 2),activation="relu", padding="same", kernel_initializer = 'he_normal')(x)
        x=BatchNormalization()(x)
        x=Add()([x, x4])
        x=Conv2DTranspose(64, (4,4), strides=(2, 2),activation="relu", padding="same", kernel_initializer = 'he_normal')(x)
        x=BatchNormalization()(x)
        x=Add()([x, x3])
        x=Conv2DTranspose(32, (4,4), strides=(2, 2),activation="relu", padding="same", kernel_initializer = 'he_normal')(x)
        x=BatchNormalization()(x)
        x=Add()([x, x2])
        x=Conv2DTranspose(16, (4,4), strides=(2, 2),activation="relu", padding="same", kernel_initializer = 'he_normal')(x)
        x=BatchNormalization()(x)
        x=Add()([x, x1])
        return x
    def build_segmentation_branch(self, outputs_4,x1,x2,x3,x4):
        x= self.make_default_hidden_layers(outputs_4,x1,x2,x3,x4)
        x=Conv2DTranspose(1, (1,1), strides=(1, 1), padding="same",\
                          activation="sigmoid", kernel_initializer = 'he_normal',name="segmentation_output")(x)
        return x
    def build_embedding_branch(self, outputs_4,x1,x2,x3,x4):
        x= self.make_default_hidden_layers(outputs_4,x1,x2,x3,x4)
        x=Conv2DTranspose(3, (1,1), strides=(1, 1), padding="same",\
                          activation="relu", kernel_initializer = 'he_normal',name="embedding_output")(x)
        return x
    def build(self,height,width):
        """
        Used to assemble our multi-output model CNN.
        """
        input_shape = (height, width, 3)
        inputs = Input(shape=input_shape)
        outputs_4,x1,x2,x3,x4=self.make_encoder(inputs)
        segmentation_branch = self.build_segmentation_branch(outputs_4,x1,x2,x3,x4)
        embedding_branch = self.build_embedding_branch(outputs_4,x1,x2,x3,x4)
        model = Model(inputs=inputs,
                     outputs = [segmentation_branch, embedding_branch],
                     name="DCNet")
        return model