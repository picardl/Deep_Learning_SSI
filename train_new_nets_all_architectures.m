function train_new_net_clst(param_path)
%This script simulates a set of training images and uses them to train an
%autoencoder network to reduce the dimensionality of the data and a
%feedforward network to classifiy central sites as occupied or unoccupied

%The script generates the following .mat data files:
%     neural_net_training_data - containing training images and
%       corresponding 3 by 3 lattice occupation patterns
%     weights_3by3 - containing fine-tuned layer weights for autoencoder,
%       where w1-w4 are encoding layers and w5-w8 are decoding.
%     classifier_net - containing a trained neural network object for encoded 
%       image classification
%     error_3by3 - containing test and training errors of the autoencoder 
%       for each training batch

%clear all
close all 

load(param_path);
cd(params.savefolder)

disp(params)

%Generate training data
sim_3by3_clst(params)

%Determine threshold for naive reconstruction and display naive fidelity
naive_matrix_experimental

%Train and save autoencoder
%deepauto

clear pic_data

%Train and save classifier
%classifier_integrated