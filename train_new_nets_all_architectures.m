function train_new_nets_all_architectures(param_path)
%This script simulates a set of training images and uses them to train a
%series of models to classify 3x3 lattice 

%The script generates the following .mat data files:
%     training_data - training images and corresponding 3 by 3 lattice 
%       occupation patterns
%     pointspreadfunction - empirical isolated atom PSF
%     gaussian_net - trained neural network model based fitting
%       gaussians to all 512 possible lattice configurations
%     three_layer_net - trained three-layer feedforward classifier
%     convolutional_net - trained convolutional network classifier

close all 

load param_path params
cd(params.savefolder)

disp(params)

%Generate training data
sim_3by3_clst(params)

%Train naive, gauss fit and gauss net models. Save gauss net
gauss_naive_and_net

%Train and save feedforward three-layer network
trilayer_feedforward

%Train and save convolutional neural network

clear pic_data