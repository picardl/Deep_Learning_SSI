%This script defines and trains a feedforward network with 1 visible layer,
%1 hidden layer and 1 output to classify the occupancy of 3x3 segments of
%optical lattice images
%
%Requires data files:
%   params: struct containing imaging simulation parameters
%   training_data.mat: .mat file containing arrays of simulated 3x3 lattice
%       images and corresponding occupation patterns

load_data %Script to load simulation data and parameters
gen_PSF %Generate point spread function, or load it if already saved

%Create new network, with 1 visible and 1 hidden layer
net  = network(1, 2, [0; 1], [1; 0], [0 0;1 0], [0 1]);

%Define transfer functions for hidden and output layers
net.layers{1}.transferFcn = 'purelin';
net.layers{2}.transferFcn = 'logsig';

net.inputs{1}.processFcns = {}; %No preprocessing function

%Number of neurons in each layer
net.input.size = size(pic_data,1); %Number of elements in input vector
net.layers{1}.size = 1;
net.layers{2}.size = 1;

%Initialise each layer individually
net.initFcn = 'initlay';

%Initialise weight matrices and bias vectors to naive values
net.IW{1} = PSF_row;
net.LW{2,1} = best_s;
net.b{2} = -5;

%Training algorithm and loss function
net.trainFcn = 'traincgp';
net.performFcn = 'crossentropy';

%Choose random training, test and validations sets
net.divideFcn = 'dividerand';

%Plot performance during  training
%net.plotFcns = {'plotperform'}; %REMOVE FOR CLUSTER

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

%Train neural network
[net, tr] = train(net, training, target);

%Determine generalisation error, with binary cost function rather than MSE
GE_test = net(GE_pics);
GE = mean(abs(round(GE_test) - GE_patt));

mean_fidel = (1 - GE)*100;

fprintf('\n Three-layer network fidelity (manual initialisation): %.3f\n', mean_fidel)
fprintf('Best training performance: %.3f\n', tr.best_perf)
fprintf('Best validation performance: %.3f\n', tr.best_vperf)
fprintf('Best testing performance: %.3f\n\n', tr.best_tperf)

save three_layer_net net