%This script is designed to test lattice reconstruction using simulated
%images

%Path to saved neural network
net_path_1 = 'C:\Users\lewis\Documents\Current Work\Innsbruck HR Imaging\Reconstruction\convolutional_net.mat';
net_path_2 = 'C:\Users\lewis\Documents\Current Work\Innsbruck HR Imaging\Trained Neural Nets\erbium_266_1-5us\classifier_net.mat';
net_path_3 = 'C:\Users\lewis\Documents\Current Work\Innsbruck HR Imaging\Trained Neural Nets\erbium_266_1-5us\single_visible_net.mat';

min_confidence = 0.5;
load 'C:\Users\lewis\Documents\Current Work\Innsbruck HR Imaging\Trained Neural Nets\erbium_266_1-5us\params_nn_integrated_id3300.mat'
%SIMULATION PARAMETERS 
sites = 10; %simulated pattern has dimensions sites x sites
fraction_filled = [0.1 0.5 0.9]; %Fraction of lattice sites occupied
num_pics = 5; %Number of pictures to simulate in each batch
NA = params.NA; %Numerical aperture
latticespacing = params.latticespacing; %in nm
imagingpulse = params.imagingpulse; %in s
lossrate = params.lossrate; %loss probability per scattering event
recoilvel = params.recoilvel; %scattering recoil velocity in m / s
scatteringrate = params.scatteringrate; %per s
pixelspersite = params.pixelspersite; %Distance between lattice sites in px
atomicmass = params.atomicmass; %Mass in amu of lattice atom
lambda = params.lambda; %Imaging wavelength
addnoise = params.addnoise; %Noise intensity (0 to 1)
collectedphotonsratio = params.collectedphotonsratio;

%Number of different conditions in which to test network
num_conditions = 3;

%SIMULATE IMAGES AND RECONSTRUCT OCCUPATIONS

%Array to hold fidelities
fidelities_1 = zeros(num_conditions, 1);
fidelities_2 = zeros(num_conditions, 1);
fidelities_3 = zeros(num_conditions, 1);

for i=1:num_conditions
    
    %Randomly generated lattice occupation
     patterns = round( rand(sites,sites, num_pics) + (fraction_filled(i) - 0.5));
    
     %Array to hold testing pictures
    test_pics = zeros(sites*pixelspersite, sites*pixelspersite, num_pics);
     
    %Simulate pictures
    for j = 1:num_pics
        test_pics(:,:,j) = simulate_setpattern( patterns(:,:,j), pixelspersite, 1,...
            NA, latticespacing, imagingpulse, lossrate,recoilvel, scatteringrate,...
            atomicmass, lambda, addnoise, collectedphotonsratio,0);
    end
    
    %Binarize simulated images
    test_pics(test_pics > 0) = 1;
    
    %Pad images with a border of zeros
    test_pics = padarray(test_pics, [pixelspersite pixelspersite 0], 0,'both');
    
    %Generate lattice coordinates
    mask = ones(size(test_pics(:,:,1)));
    [ lattice_coords, lattice_indices ] = gen_lattice( test_pics(:,:,1),...
        pixelspersite, pixelspersite, pixelspersite, mask);
    
    %Reconstruct lattice
    [binary_1, ~, ~] = reconstruct_lattice(test_pics, lattice_coords, lattice_indices,...
        pixelspersite, 0.5, net_path_1, true, -1);
    fidelities_1(i) = 1 - mean(mean(mean(abs(binary_1 - patterns))));
    
    [binary_2, ~, ~] = reconstruct_lattice(test_pics, lattice_coords, lattice_indices,...
        pixelspersite, 0.5, net_path_2, true, -1);
    fidelities_2(i) = 1 - mean(mean(mean(abs(binary_2 - patterns))));
    
    [binary_3, ~, ~] = reconstruct_lattice(test_pics, lattice_coords, lattice_indices,...
        pixelspersite, 0.5, net_path_3, true, -1);
    fidelities_3(i) = 1 - mean(mean(mean(abs(binary_3 - patterns))));
end
    
all_fidels = [fidelities_1(1:num_conditions),fidelities_2(1:num_conditions),fidelities_3(1:num_conditions)];
b = bar(fraction_filled, all_fidels*100, 0.9);
b(1).FaceColor = [0.12 0.56 1];
b(2).FaceColor = [1 0.6 0];
b(3).FaceColor = [1 0.9 0.3];
ylim([85, 100]);
xlim([0,1]);
ax = gca;
ax.YTick = [90 95 100];
ax.XTick = [0.1 0.3 0.5 0.7 0.9];
xlabel('Lattice filling fraction');
ylabel('Reconstruction fidelity (%)');
legend('Convolutional', 'Autoencoder', 'Three-layer');

auto_falsepos = 1 - mean(mean(mean(abs(binary_1 - patterns))));
        