%This script implements four versions of optical lattice image
%reconstruction, beginning with manual naive reconstruction, and building
%on this method using three possible neural network architectures for image
%classification.
%Requires data files:
%   params: struct containing imaging simulation parameters
%   training_data.mat: .mat file containing arrays of simulated 3x3 lattice
%       images and corresponding occupation patterns





%-------------------------------------------------------------------------
%FIND BEST CLASSIFICATION FIDELITY USING JUST KNOWN PSF
%-------------------------------------------------------------------------

%Integer number of images making up approx 10% of dataset
%Used to select subset of data for calculating generalisation error
N_GE = round(N/10);

%Generate N_GE random indices to remove from training data
training_ind = linspace(1,N,N);
GE_ind = randperm(N, N_GE);
training_ind(GE_ind) = [];

%Generate training and testing datasets
training = pic_data(:, training_ind);
target = nn_patts(:, training_ind);
GE_pics = pic_data(:, GE_ind);
GE_patt = nn_patts(:, GE_ind);

%Output of layer 1 of network applied to training data
L1 = PSF_row*training;
train_err = [];

%Estimate maximum necessary range of hidden neuron weights to test
s_max = 20/( (sum(sum( GE_pics(:,1:5) ))/5)*mean(PSF_row));

%Generate array of 4000 evenly spaced values to test
s_arr = linspace(0, s_max, 4000);

%Determine average training error over all simulated images
for i = 1:4000
    L2 = 1./(1 + exp(-L1*s_arr(i) + 5));
    train_err = [train_err, mean(abs(round(L2) - target))];
end

%Find value of multiplier s which minimizes error
[min_err, indmin] = min(train_err);
best_s = s_arr(indmin);

%Determine generalization error
L1_GE = PSF_row*GE_pics;
L2_GE = 1./(1 + exp(-L1_GE*best_s + 5));
GE = mean(abs(round(L2_GE) - GE_patt));

fprintf('\nMax naive fidelity (training): %.3f\n', 100*(1-min_err))
fprintf('Max naive fidelity (generalization): %.3f\n', 100*(1-GE))
fprintf('\nWith layer weight of: %.3f\n', best_s)

%-------------------------------------------------------------------------
%CREATE AND TRAIN BILAYER NEURAL NETWORK
%-------------------------------------------------------------------------

%Create new network, with 1 visible and 1 hidden layer
net  = network(1, 2, [0; 1], [1; 0], [0 0;1 0], [0 1]);

%Define transfer functions for hidden and output layers
net.layers{1}.transferFcn = 'purelin';
net.layers{2}.transferFcn = 'logsig';

net.inputs{1}.processFcns = {}; %No preprocessing function

%Number of neurons in each layer
net.input.size = size(training,1); %Number of elements in input vector
net.layers{1}.size = 1;
net.layers{2}.size = 1;

%Initialise each layer individually
net.initFcn = 'initlay';

%Initialise weight matrices and bias vectors to naive values
net.IW{1} = PSF_row;
net.LW{2,1} = best_s;
net.b{2} = -5;

%Training algorithm
net.trainFcn = 'traincgp';

net.performFcn = 'crossentropy';

%Choose random training, test and validations sets
net.divideFcn = 'dividerand';

%Plot performance during  training
%net.plotFcns = {'plotperform'};

net.trainParam.showWindow = false;

%Set training parameters
net.performParam.normalization = 'standard';
net.performParam.regularization = 0;
net.trainParam.epochs = 1000;
net.trainParam.max_fail = 30;

%Integer number of images making up approx 10% of dataset
%Used to select subset of data for calculating generalisation error
N_GE = round(N/10);

%Generate N_GE random indices to remove from training data
training_ind = linspace(1,N,N);
GE_ind = randperm(N, N_GE);
training_ind(GE_ind) = [];

%Generate training and testing datasets
training = pic_data(:, training_ind);
target = nn_patts(:, training_ind);
GE_pics = pic_data(:, GE_ind);
GE_patt = nn_patts(:, GE_ind);

%Configure neural net using first element of dataset as example
%net = configure(net, training(:, 1), target(:, 1));

breakout = 1; %counter to break out of while loop
GE = 1; %Reset generalization error

%Train neural network, retraining up to 5 times if performance is worse than naive
while (GE > min_err) && (breakout < 6)

	%Train neural network
	[net, tr] = train(net, training, target);

	%Determine generalisation error, with binary cost function rather than MSE
	GE_test = net(GE_pics);
	GE = mean(abs(round(GE_test) - GE_patt));

	breakout = breakout + 1;
end

mean_fidel = (1 - GE)*100;

fprintf('\Three-layer network fidelity (manual initialisation): %.3f\n', mean_fidel)
fprintf('Best training performance: %.3f\n', tr.best_perf)
fprintf('Best validation performance: %.3f\n', tr.best_vperf)
fprintf('Best testing performance: %.3f\n\n', tr.best_tperf)

save three_layer_net net


%-------------------------------------------------------------------------
%FIND MAXIMUM FIDELITY USING GAUSSIAN FIT
%-------------------------------------------------------------------------
gauss_center_fits = zeros(1,N);
gauss_center_confs = zeros(2,N);
gauss_center_sum = zeros(1,N);
gauss_all_fits = zeros(9,N);

for j = 1:N
    [gausamps, gaussconf, fitZ] = GaussiansAmplitudeFit(square_pics(:,:,1,j),pixelspersite, PSFwidth);
    gauss_center_fits(j) = gausamps(5);
    gauss_all_fits(:,j) = gausamps;
    gauss_center_confs(:,j) = gaussconf(:,5);
end

training_fits = gauss_center_fits(:, training_ind);
GE_fits = gauss_center_fits(:, GE_ind);

t_max = max(max(training_fits));
t_arr = linspace(0,t_max, 4000);

train_err_gauss = [];

for i = 1:4000
    t = t_arr(i);
    fitout = training_fits;
    fitout(fitout < t) = 0;
    fitout(fitout >= t) = 1;
    train_err_gauss = [train_err_gauss, mean(mean(abs(fitout - target)))];
end

[min_err_gauss, indmin] = min(train_err_gauss);
best_t = t_arr(indmin);

GE_fits(GE_fits < best_t)  = 0;
GE_fits(GE_fits >= best_t) = 1;
GE_gauss = mean(abs(round(GE_fits) - GE_patt));

fprintf('\nMax Gaussian fit fidelity (training): %.3f\n', 100*(1-min_err_gauss))
fprintf('Max Gaussian fit fidelity (generalization): %.3f\n', 100*(1-GE_gauss))
fprintf('\nWith threshold of: %.3f\n', best_t)

%-------------------------------------------------------------------------
%CREATE NETWORK INCLUDING EACH POSSIBLE COMBINATION OF GAUSSIAN FITS
%-------------------------------------------------------------------------

%Create new average PSF for isolated atom in top left corner
PSF_data = simulate_setpattern([1,0,0; 0,0,0;0,0,0], pixelspersite, 10000, NA, latticespacing, imagingpulse, lossrate, recoilvel, scatteringrate, atomicmass, lambda, 0, 0.99, 0);
PSF_corner = sum(PSF_data, 3);

PSF_int = sum(sum(PSF_corner));
PSF_corner = PSF_corner/PSF_int;

clear PSF_data;

%Create new average PSF for isolated atom on top edge
PSF_data = simulate_setpattern([0,1,0; 0,0,0;0,0,0], pixelspersite, 10000, NA, latticespacing, imagingpulse, lossrate, recoilvel, scatteringrate, atomicmass, lambda, 0, 0.99, 0);
PSF_edge = sum(PSF_data, 3);

PSF_int = sum(sum(PSF_edge));
PSF_edge = PSF_edge/PSF_int;

site_specific_PSF = zeros(round(3*pixelspersite),round(3*pixelspersite),9);

%Generate subsequent rows by rotating PSFs centered on corner and edge
for orientation = 1:2
	site_specific_PSF(:,:,2*orientation - 1) = PSF_corner;
    site_specific_PSF(:,:,2*orientation) = PSF_edge;
    
	PSF_corner = rot90(PSF_corner);
	PSF_edge = rot90(PSF_edge);
end
site_specific_PSF(:,:,5) = PSF;
for orientation = 3:4
	site_specific_PSF(:,:,2*orientation - 1) = PSF_corner;
    site_specific_PSF(:,:,2*orientation) = PSF_edge;
    
	PSF_corner = rot90(PSF_corner);
	PSF_edge = rot90(PSF_edge);
end

%Create new network, with 1 visible and 1 hidden layer
net  = network(1, 2, [1; 1], [1; 0], [0 0;1 0], [0 1]);

%Define transfer functions for hidden and output layers
net.layers{1}.transferFcn = 'purelin';
net.layers{2}.transferFcn = 'logsig';

net.inputs{1}.processFcns = {}; %No preprocessing function

%Number of neurons in each layer
net.input.size = size(pic_data,1); %Number of elements in input vector
net.layers{1}.size = 512;
net.layers{2}.size = 1;

%Initialise each layer individually
net.initFcn = 'initlay';

%Initial value for layer 1, weighting filled central sites positively and
%empty negatively
init_L1 = zeros(1,512);
init_L1(1,1:256) = -1;
init_L1(1,257:512) = 1;

%Initialize layer 1 bias to zero
net.b{1} = zeros(512,1);

%Initialise weight matrices and layer 2 bias vector to naive values

combs = dec2base(0:power(2,9) - 1,2) - '0';
empty = combs(combs(:,5) == 0,:);
filled = combs(combs(:,5) == 1,:);
all_input_weights = zeros(512,round(3*pixelspersite)^2);
indices = repmat(1:9,256,1);
empty = empty.*indices;
filled = filled.*indices;

for i = 1:256
    mat1_all = site_specific_PSF(:,:,nonzeros(empty(i,:)));
    mat1_sum = sum(mat1_all, 3);
    all_input_weights(i,:) = reshape(mat1_sum,1,[]);
    
    mat2_all = site_specific_PSF(:,:,nonzeros(filled(i,:)));
    mat2_sum = sum(mat2_all, 3);
    all_input_weights(i+256,:) = reshape(mat2_sum,1,[]);
end
net.IW{1} = all_input_weights;
net.LW{2,1} = init_L1;
net.b{2} = 1;

%Training algorithm
net.trainFcn = 'traincgp';

net.performFcn = 'crossentropy';

%Choose random training, test and validations sets
net.divideFcn = 'dividerand';

%Plot performance during  training
%net.plotFcns = {'plotperform'};

net.trainParam.showWindow = false;

%Set training parameters
net.performParam.normalization = 'standard';
net.performParam.regularization = 0;
net.trainParam.epochs = 1000;
net.trainParam.max_fail = 30;

%Integer number of images making up approx 10% of dataset
%Used to select subset of data for calculating generalisation error
N_GE = round(N/10);

%Generate N_GE random indices to remove from training data
training_ind = linspace(1,N,N);
GE_ind = randperm(N, N_GE);
training_ind(GE_ind) = [];

%Generate training and testing datasets
training = pic_data(:, training_ind);
target = nn_patts(:, training_ind);
GE_pics = pic_data(:, GE_ind);
GE_patt = nn_patts(:, GE_ind);

%Configure neural net using first element of dataset as example
%net = configure(net, training(:, 1), target(:, 1));

breakout = 1; %counter to break out of while loop
GE = 1; %Reset generalization error

%Train neural network
[net, tr] = train(net, training, target);

%Determine generalisation error, with binary cost function rather than MSE
GE_test = net(GE_pics);
GE = round(GE_test) - GE_patt;

mean_fidel = (1 - mean(abs(GE)))*100;

fprintf('\nNetwork with 512-neuron visible layer fidelity: %.3f\n', mean_fidel)
fprintf('Best training performance: %.3f\n', tr.best_perf)
fprintf('Best validation performance: %.3f\n', tr.best_vperf)
fprintf('Best testing performance: %.3f\n\n', tr.best_tperf)

save gaussian_net net

clear('pic_data', 'training', 'target')

%-------------------------------------------------------------------------
%CREATE CONVOLUTIONAL NEURAL NETWORK
%-------------------------------------------------------------------------
%Generate N_GE random indices to remove from training data
N_val = round(N/5);
val_perm = randperm(size(training_ind,2), N_val);
val_ind = training_ind(val_perm);
training_ind(val_perm) = [];

GE_ind = val_ind(1:round(0.2*N_val));
val_ind((N_val-round(0.2*N_val)):end) = [];

%Generate training and validation datasets
training = square_pics(:,:,1, training_ind);
target = categorical(nn_patts(:, training_ind));
val_pics = square_pics(:,:,1, val_ind);
val_patt = categorical(nn_patts(:, val_ind));
validation_set = {val_pics,val_patt};

GE_pics = square_pics(:,:,1, GE_ind);
GE_patt = nn_patts(:, GE_ind);

%Define architecture of convolutional neural net
layers = [
    imageInputLayer([dim1 dim1 1])
    
    convolution2dLayer(filtersize,numfilters)
    batchNormalizationLayer
    reluLayer
    
    averagePooling2dLayer(poolsize,'Stride',stride)
    
    convolution2dLayer(filtersize,numfilters*2)
    batchNormalizationLayer
    reluLayer
    
    averagePooling2dLayer(poolsize,'Stride',stride)
    
    convolution2dLayer(filtersize,numfilters*5)
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];

options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',4, ...
    'Shuffle','every-epoch', ...
    'ValidationData',validation_set, ...
    'ValidationFrequency',30, ...
    'Verbose',true, 'VerboseFrequency', 500);

net = trainNetwork(training,target,layers,options);

GE_out = classify(net, GE_pics);
binary_out = grp2idx(GE_out) - 1;
errs = abs(binary_out' - GE_patt);
fidelity = 1 - mean(errs);

fprintf('\nConvolutional neural network fidelity: %.3f\n\n', fidelity*100)

save convolutional_net net