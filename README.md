# Dispersion_Curves
This project is within my internship at OMP, Toulouse. The purpose of this subject is to find a machine learning/deep learning method to extract dispersion curves from dispersion images, which is especially important for inversion processus to study underground structure. The achitecture of our neural net is composed of 2 U-nets in parallel, that share a commun part in the encoder. In this project, you have codes to create the database and train the model.

├── __Code_generation_database__: Folder which contains scripts to create a new sample\
│   ├── __Test_src__: virgin parameter files for SPECFEM2D and CPS. Do not change them except if you know what you are doing \
│   ├── __alleger_XX.py__: Erases receiver's output files in the horizontal direction \
│   ├── __automatisation_generation.py__: Main code to generate a sample \
│   ├── __automatisation_generation_nrec_inferieur.py__: Used to create a new sample with another one, by only changing the number of receivers used to make the                                                            dispersion image \
│   ├── __dispersion.py__: Makes the dispersion matrix, which is the input of the model. \
│   ├── __ecriture_totale.py__: Codes that automatically change parameter files to create the simulation you want \
│   ├── __generation_automatique.py__: Allows processing several samples \
│   ├── __nombre_mode.py__: Creates an embedding and segmentation image with the number of modes you want, to train the model with the great number of modes\
│   ├── __nombre_mode_fast.py__: Creates an embedding and segmentation image with the number of modes you want. Faster but needs to have the 6 modes version                                        already created to run. \
│   ├── __run_specfem.sh__: The bash used to launch specfem2d. Do not modify it except if you know what you are doing \
│   ├── __segmentation-embedding.py__: Creates segmentation and embedding image, initially with 6 modes, which are the outputs of the model \
│   ├── __sortir_bon_mode.py__: Change in the Data folder, for each sample, the number of modes that the model will use. The segmentation and embedding matrix                                     with the number of modes you wants need to be created before with nombre_mode.py or nombre_mode_fast.py codes. \
│   └── __superposition_CPS_specfem.py__: Creates an image of the dispersion image superposed with CPS predictions. \
├── __Code_semi_supervised_learning__: Folder which contains shortcodes to simplify semi-supervised learning\
│   ├── __DC_net_to_ML_data.py__: Put predictions on the database\
│   ├── __afficher_curve.py__: Displays dispersion curves that have been predicted \
│   └── __trier_predictions.py__: Manual assisted sorting to only choose predictions that seem to be correct  \
├── __Code_traitement_database__: Folder which contains shortcodes to modify or check the database \
│   ├── __Output_to_output_data.py__: Put all the output files in a different folder than the database, to lighten it \
│   ├── __alleger_XX_int.py__: Erases receiver's output files in the horizontal direction, for several samples. Possible to use when their names are int. \
│   ├── __bon_nombre.py__: Changes sample's folder names in a database to have samples folder names from 1 to number_of_samples, to train greatly.  \
│   ├── __decalage_numero.py__: Shifts the sample's folder names \
│   ├── __enlever_erreur_cps.py__: Find samples where CPS had a bug, and erase these sample folders \
│   ├── __enlever_mode.py__: Erase segmentation and embedding images with several modes on them which is no longer useful \
│   ├── __extraire_vp_vs.py__: Creates an image to visualize vp and vs apportionments in the soil you generated. \
│   ├── __generer_disp_jpeg.py__: Creates a dispersion image \
│   ├── __regarder_data.py__: Allows looking at some dispersion, embedding, and segmentation images of the database \
│   ├── __verifier_data.py__: Checks if each sample in the database has a dispersion matrix, a segmentation matrix, and an embedding one. \
│   └── __visionner.py__: Allows looking at some receivers vertical outputs of a sample. \
├── __specfem2d__: SPECFEM program folder \
├── __PROGRAMS.330__: CPS program folder \
├── __DCNet_output__: The folder where all the predictions are stocked\
├── __ML_DATA__: The folder to stock all the database \
├── __Feed_Data.py__: help to feed data into a model for training \
├── __generate_predict.py__: Used to generate predictions (the one with the segmentation clustering) for several samples \
├── __loss.py__: defining loss function for the model \
├── __predict.py__: using model to predict \
├── __predict_with_fv.py__: using the model to predict (adding f and v in DBSCAN) \
├── __predict_with_seg_clusters.py__: using the model to predict (adding segmentation clustering in DBSCAN) \
├── __requirements.txt__: packages/libraries necessary for our project \
├── __structure.py__: defining the architecture of the neural net \
├── __train.py__: training model \
└── __train_validate_test_split.py__: spliting train-test-validate, creating 3 file: __test_idx.csv__, __train_idx.csv__, __valid_idx.csv__ 

## Some instructions for using files:
For these files to work correctly, please put on the files in the same directory. To use data, please contact me or Mr. Roland Martin (OMP) for the directory called __ML_DATA__. Once you have this __ML_DATA__, put the folder __ML_DATA__ in the same directory as all the files above. 

Before working with these files, please unzip __checkpoints_best_only.zip__ and put it in the same directory as all the files above. 

### Preliminary step: set up a virtual environment to work with python
Create a virtual environment :\
- with conda :
```
conda create --name DCNet python=3.10
```
- with Python:
```
python3 -m venv DCNet
```
Notice that the environment here is named DCNet, you can choose a name at your convenience.

To work in this environment, we need to activate it:
- with conda :
```
conda activate DCNet
```
- with Python:
```
source DCNet/bin/activate
```
Now this environment is empty. We load packages/libraries necessary for our project:
- with conda :
```
conda install --file requirements.txt
```
- with Python:
```
pip install -r requirements.txt
```
When you are in this virtual environment (in other words, the DCNet environment is activated), you can work with .py files above.
After performing all the tasks and you want to quit this environment, type:
- with conda :
```
conda deactivate
```
- with Python:
```
deactivate
```
Next, follow the instructions below:
## Build DATA
### Importing CPS and specfem2d programs
If you have the ML_DATA folder, you can train the DCNet model or make predictions on this dataset. However, this project has integrated an automatic way to generate samples, which can be added to the dataset or made to be predicted by the model. This part runs only in Linux, so if you are in a Windows environment, you can use WSL and if you are in an OS environment, you can use a virtual Unix environment.

First, you must download the specfem2d and CPS algorithm in this folder to use this feature.
For specfem2d, type in the terminal command : 
```
git clone --recursive --branch devel https://github.com/SPECFEM/specfem2d.git
```
Go to the specfem2d folder and type :
```
./configure FC=gfortran CC=gcc
```
If you want to run in parallel, i.e., using more than one processor core, then you would type
```
./configure FC=gfortran CC=gcc MPIFC=mpif90 --with-mpi
```
Finally, after choosing your configuration :
```
make
```

For CPS, install it in the main folder by following the installation process:  http://www.eas.slu.edu/People/RBHerrmann/ComputerPrograms.html

You can find in the doc subfolder the complete documentation to install these two programs (manual_SPECFEM2D and manual_CPS).

Finally, you must have two subfolders specfem2d and PROGRAMS.330 in your main folder, which each runs correctly.

### Generating data
The main code for generating data is automatisation-generation.py, in the Code_generation_database subfolder. You can find in the doc subfolder an explanation of how this code works, and which codes it uses to generate the database.

Warning: if you use the topographic features to make non-flat soils, the CPS result will not be accurate because CPS only process flat layers model: so you can only use these samples to make predictions on them or to make semi-supervised learning, but don't directly train your model on these samples!

## DCNet

### First, before training, we need to split the train-validate-test.
We use the file train_validate_test_split.py. In the terminal command, type:
```
python train_validate_test_split.py --total=5746 --ratio=0.8
```
Where: total is the total number of samples, in our case, 5746; ratio is the split ratio, by default, it gets 0.8. Notice that:
<img src="https://render.githubusercontent.com/render/math?math=ratio \in ]0,1["> \
By doing this, it will create 3 files: __test_idx.csv, train_idx.csv, valid_idx.csv__, that contain indexes of samples for testing, training, and validating.
In this repository, you do not need to do this step as I already upload these 3 files.
Another important notice is that you just do it only one time before training because it will give different results randomly each time.
### For training:
In the terminal command, type if you want to train your model from scratch:
```
python train.py --epochs=75 --batch_size=32 --learning_rate=0.0005 --data_dir=./ML_DATA
```
In the terminal command, type if you have a pre-trained model:
```
python train.py --epochs=5 --batch_size=32 --learning_rate=0.0005 --data_dir=./ML_DATA --weight_dir=./checkpoints_best_only/checkpoint
```
In the terminal command, type if you want to train your model in the background:
```
nohup python train.py --epochs=75 --batch_size=32 --learning_rate=0.0005 --data_dir=./ML_DATA &
```
Where: epochs: number of epochs; data_dir: directory containing data; weight directory: the directory containing pre-trained weight.
#### By default:
epochs=5, batch_size=32, learning_rate=0.0005, data_dir=./ML_DATA, weight_dir=./checkpoints_best_only/checkpoint

To avoid overfitting, we save only the best model based on its performance on validating set. The weight is saved in the directory ./checkpoints_best_only/checkpoint.

At the end of the training, it will save a file called __loss_train_val.csv__ containing the evolution of the loss function on the training and validating set.
Also, this will generate a file called metrics.csv containing the evolution of metrics on training and validating sets. \
metrics.csv : 
- 1er line : segmentation_output_accuracy
- 2nd line : segmentation_output_dice_coef
- 3rd line : segmentation_output_mean_io_u
- 4th line : embedding_output_accuracy
- 5th line : val_segmentation_output_accuracy
- 6th line : val_segmentation_output_dice_coef
- 7th line : val_segmentation_output_mean_io_u
- 8th line : val_embedding_output_accuracy
- 9th line : segmentation_output_recall
- 10th line : val_segmentation_output_recall
- 11th line : segmentation_output_precision
- 12th line : val_segmentation_output_precision

### For predicting:
```
python predict.py --sample=600
```
Where the sample is the sample number we want to predict.
By typing this, the file will perform the following tasks:
1. Showing input dispersion image.
2. Showing output of the segmentation branch and at the same time saving this in __Segmentation.jpg__.
3. Showing output of embedding branch and at the same time saving this in __Embedding.jpg__.
4. Showing different modes obtained from DBSCAN (normally, we have to wait for seconds for results) and at the same time saving this in __groups.jpg__.
5. Showing fitted curved found from the maximum peak and at the same time saving this in __Fitted_Curves.jpg__.
6. Showing rectified curves after using a Gaussian filter and at the same time saving this in __Fitted_Curves_Filtered.jpg__.
7. Showing the curves on the real scale after rescaling and at the same time saving this in __Predicted_Fitted_Curves.jpg__.
8. Showing predicted curves superimposed on the initial input dispersion image and at the same time saving this in __superimposed.jpg__.
9. Saving curves under the form of .csv first column: mode number, 2nd column: phase velocity, 3rd column: frequency. The file is named __predicted_fv_curves.csv__.

Noting that all the output files above are saved in folder ./DCNet_output/number, with the number being the sample number.

* Noting also that numbers smaller than 2600 are solid, between 2601 and 5480 are a mixture of solid, water, and air, and numbers from 5481 are karsts.

Example of prediction with karst :
```
python predict.py --sample=5117
```
I also implement two other versions of prediction.

One that uses f and v in addition to the output of the embedding branch in DBSCAN. To use this program:
```
python predict_with_fv.py --sample=600 --weight_f=0.5 --weight_v=0.5
```
By default, weight_f=weight_v=0.5

One that uses segmentation clustering in addition to the output of embedding branches in DBSCAN. To use this program:
```
python predict_with_seg_clusters.py --sample=600 --nb_modes_voulu=1 --stockage_data=./ML_DATA/
```
Put at nb_modes_voulu the same number of modes as the number used to train the neural network.


## Good links to have an idea about the metrics that I used in my project :

- https://towardsdatascience.com/metrics-to-evaluate-your-semantic-segmentation-model-6bcb99639aa2 
- https://ilmonteux.github.io/2019/05/10/segmentation-metrics.html
- https://www.kaggle.com/code/yassinealouini/all-the-segmentation-metrics


