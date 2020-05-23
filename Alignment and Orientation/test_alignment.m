%This script is used to test the various lattice alignment processes
%required to process real lattice images.
%Simulate a set of lattice pictures with given paramters, find angle and
%spacing, rotate and align and generate lattice coordinates.

%Path to neural network to be used for reconstruction
net_path = 'C:\Users\lewis\Documents\Current Work\Innsbruck HR Imaging\Trained Neural Nets\accurate_kozuma_binarized\classifier_net.mat';

sites = 15; %simulated pattern has dimensions sites x sites
fraction_filled = 0.4;
init_angle = 0; %in degrees
n_pics = 3;
NA = 0.81; %Numerical aperture
latticespacing = 543.5000; %in nm
imagingpulse = 4.0000e-05; %in s
lossrate = 4.87e-5; %loss probability per scattering event
recoilvel = 0.0057; %scattering recoil velocity in m / s
scatteringrate = 13000000; %per s
pixelspersite = 5; %Distance between lattice sites in px
atomicmass = 174; %Mass in amu of lattice atom
lambda = 399; %Imaging wavelength
addnoise = 0.0520; %Noise intensity (0 to 1)
collectedphotonsratio = 0.0660;
latticedepth = 435;
angle = 0; %Angle of rotation of lattice with respect to image axes
SNR_guess = 5; %Estimate of signal to noise ratio, used for determining image bounds

%3 dimensional array containing simulated pictures, associated lattice patterns and PSF
[sim_pics, sim_patterns, PSFwidth] = simulate_npictures(sites, fraction_filled, angle, pixelspersite, n_pics, NA, latticespacing, imagingpulse, lossrate, recoilvel, scatteringrate, atomicmass, lambda, addnoise,collectedphotonsratio, 435);

%[angle, spacing] = find_angle_and_spacing(sim_pics, sites, pixelspersite, latticespacing, fraction_filled, init_angle);
angle = 0; spacing = 5;

[ crop_pics, x_offset, y_offset, mask ] = rotate_and_align(sim_pics, angle, spacing, 5);

[ lattice_coords, lattice_indices ] = gen_lattice( crop_pics(:,:,1), spacing, x_offset, y_offset, mask); %Generate lattice site coordinates

[binary, sigmoid, confidence, reconstruction_list] = reconstruct_lattice(crop_pics, lattice_coords, lattice_indices, spacing, 0.5, net_path, true, -1);

lin_bin = reshape(binary(:,:,1),[],1);
lin_pat = reshape(sim_patterns(:,:,1),[],1);
errs = find(abs(lin_bin-lin_pat)>0);

occ_firstpic = find(reconstruction_list(:,1) > 0.5);
occ_coords = lattice_coords(occ_firstpic,:);
err_coords = lattice_coords(errs,:);

mean_bright = mean(mean(mean(crop_pics(crop_pics>0))));
crop_pics = crop_pics/mean_bright; %tanh(atanh(0.5)*crop_pics/mean_bright);

figure(9) 
image(30*crop_pics(:,:,1))
hold on
set(gca, 'Units', 'Pixels');
plot(occ_coords(:,1),occ_coords(:,2),'wo','MarkerSize',15, 'LineWidth', 2);
plot(err_coords(:,1),err_coords(:,2),'rx','MarkerSize',15, 'LineWidth', 2);
hold off 

figure(10)
image(63*sigmoid(:,:,1));

figure(11)
image(63*sim_patterns(:,:,1));