%This script is designed to test lattice reconstruction for pretrained networks using simulated
%images, and compare trained neural networks to reconstruction using a fit
%of 9 Gaussians (state of the art method without deep learning)
%Paths of network, saved point spread function and simulation parameters files must be manually specified

threshold = 0.107; %Threshold fit amplitude for Gaussian reconstruction
bPlot = 0; %Bool determining whether to plot Gaussian fit for each segment

%Path to saved neural network
net_paths = {};

%Add as many networks as we want to test
net_paths{end+1} = '{user_path}\test_net.mat';
net_paths{end+1} = '{user_path}\test_net.mat';
net_paths{end+1} = '{user_path}\test_net.mat';
num_nets = length(net_paths);

%Load saved simulation params. Alternatively they can be manually defined
%below
param_path = "{user_path}\params.mat";
load(param_path, 'params')

%Load PSF
try
    load("{user_path}\pointspreadfunction.mat", 'PSF');
    PSFcoefs = gauss_fit2D(PSF, size(PSF,1)/2, size(PSF,2)/2, pixelspersite);
    PSFwidth = PSFcoefs(2);
catch
    PSFwidth = params.pixelspersite/2;
    warning('Failed to load a PSF, using PSFwidth of 1/2 lattice spacing')
end

%Number of different filling fractions in which to test network
fraction_filled = [0.5, 0.9]; %Fraction of lattice sites occupied
num_conditions = length(fraction_filled);

min_confidence = 0.5; %Minimum confidence for positive reconstruction classification, in range 0.5 - 1

%SIMULATION PARAMETERS 
sites = 10; %simulated pattern has dimensions sites x sites
num_pics = 2; %Number of pictures to simulate in each batch
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

%SIMULATE IMAGES AND RECONSTRUCT OCCUPATIONS

%Cell array to hold fidelities and reconstruction outputs
fidelities = {};
binary_recons = {};
sigmoid_recons = {};
confidence_recons = {};
recons_lists = {};


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
    
    %Reconstruct lattice using neural networks
    for j = 1:num_nets
        [binary_recons{j,i}, sigmoid_recons{j,i}, confidence_recons{j,i},recons_lists{j,i}] =...
            reconstruct_lattice(test_pics, lattice_coords, lattice_indices,...
            pixelspersite, 0.5, net_paths{j}, true, -1);
        fidelities{j,i} = 1 - mean(mean(mean(abs(binary_recons{j,i} - patterns))));
    end
    
    %Reconstruct using gaussian fit
    [binary_recons{j + 1,i}, recons_lists{j + 1,i}] = reconstruct_lattice_gaussfit(test_pics, lattice_coords, lattice_indices,...
        PSFwidth, pixelspersite, threshold, true);
     fidelities{j + 1,i} = 1 - mean(mean(mean(abs(binary_recons{j + 1,i} - patterns))));

end

disp('Fidelities of all networks tested (rows) in each condition (columns):')
disp(fidelities)

%Plot fidelities for different filling fractions
if length(fraction_filled) > 1
    cat_fidels = zeros(num_conditions, num_nets);
    for i = 1:(num_nets + 1)
        for j = 1:num_conditions
            cat_fidels(j,i) = fidelities{i,j};
        end
    end
    bar(fraction_filled, cat_fidels*100, 0.9);
    ylim([85, 100]);
    xlim([0,1]);
    ax = gca;
    ax.YTick = [90 95 100];
    ax.XTick = [0.1 0.3 0.5 0.7 0.9];
    xlabel('Lattice filling fraction');
    ylabel('Reconstruction fidelity (%)');
    %legend();
end

        