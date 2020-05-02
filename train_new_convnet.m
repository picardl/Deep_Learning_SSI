function train_new_convnet(param_path)
%This script simulates a set of training images and uses them to train a
%convolutional neural network to classify 3x3 lattice images

%The script generates the following .mat data files:
%     training_data - training images and corresponding 3 by 3 lattice 
%       occupation patterns
%     pointspreadfunction - empirical isolated atom PSF
%     convolutional_net - trained convolutional network classifier

close all 

load param_path params
cd(params.savefolder)

disp(params)

%Generate training data
sim_3by3_clst(params)

%Train and save convolutional neural network
convolutional

clear pic_data