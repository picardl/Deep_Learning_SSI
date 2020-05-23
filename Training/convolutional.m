%This script defines and trains a convolutional network to classify the 
%occupancy of 3x3 segments of optical lattice images
%
%Requires data files:
%   params: struct containing imaging simulation parameters
%   training_data.mat: .mat file containing arrays of simulated 3x3 lattice
%       images and corresponding occupation patterns

load_data %Script to load simulation data and parameters

%Integer number of images making up approx 10% of dataset
%Used to select subset of data for calculating generalisation error
N_GE = round(N/10);

%Generate N_GE random indices to remove from training data
training_ind = linspace(1,N,N);
GE_ind = randperm(N, N_GE);
training_ind(GE_ind) = [];

%Generate N_val random indices to remove from training data for validation
%and testing
N_val = round(N/5);
val_perm = randperm(size(training_ind,2), N_val);
val_ind = training_ind(val_perm);
training_ind(val_perm) = [];

%Set aside another 5th of those images for post-training generalization
%error test
GE_ind = val_ind(1:round(0.2*N_val));
val_ind((N_val-round(0.2*N_val)):end) = [];

%Generate training and validation datasets
training = square_pics(:,:,1, training_ind);
target = categorical(nn_patts(:, training_ind));
val_pics = square_pics(:,:,1, val_ind);
val_patt = categorical(nn_patts(:, val_ind));
validation_set = {val_pics,val_patt};

%Generate generalization dataset
GE_pics = square_pics(:,:,1, GE_ind);
GE_patt = nn_patts(:, GE_ind);

%Define architecture of convolutional neural net
layers = [
    imageInputLayer([dim1 dim1 1])
    
    %First convolutional layer with batchNorm and reLu
    convolution2dLayer(filtersize,numfilters)
    batchNormalizationLayer
    reluLayer
    
    %Average pooling
    averagePooling2dLayer(poolsize,'Stride',stride)
    
    %Second convolutional layer with batchNorm and reLu
    convolution2dLayer(filtersize,numfilters*2)
    batchNormalizationLayer
    reluLayer
    
    %Average pooling
    averagePooling2dLayer(poolsize,'Stride',stride)
    
    %Third convolutional layer with batchNorm and reLu
    convolution2dLayer(filtersize,numfilters*5)
    batchNormalizationLayer
    reluLayer
    
    %Fully connected output layer with softmax
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];

%Training parameters
options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',4, ...
    'Shuffle','every-epoch', ...
    'ValidationData',validation_set, ...
    'ValidationFrequency',30, ...
    'Verbose',true, 'VerboseFrequency', 500);

%Train network
net = trainNetwork(training,target,layers,options);

%Check generalization error
GE_out = classify(net, GE_pics);
binary_out = grp2idx(GE_out) - 1;
errs = abs(binary_out' - GE_patt);
fidelity = 1 - mean(errs);

fprintf('\nConvolutional neural network fidelity: %.3f\n\n', fidelity*100)

save convolutional_net net