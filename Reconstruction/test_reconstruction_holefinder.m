% %This script is designed to test lattice reconstruction using simulated
% %images
% 
%Path to saved neural network
net_path_1 = 'C:\Users\lewis\Documents\Current Work\Innsbruck HR Imaging\Trained Neural Nets\Greiner\classifier_net_filling_id9076.mat';
% 
% min_confidence = 0.5;
% 
% %SIMULATION PARAMETERS
% sites = 5; %simulated pattern has dimensions sites x sites
% num_pics = 300; %Number of pictures to simulate in each batch
% NA = 0.89; %Numerical aperture
% latticespacing = 266; %in nm
% imagingpulse = 3e-6; %in s
% lossrate = 1e-4; %loss probability per scattering event
% recoilvel = 5.9e-3; %scattering recoil velocity in m / s
% scatteringrate = 2*pi*30e6; %per s
% pixelspersite = 8; %Distance between lattice sites in px
% atomicmass = 168; %Mass in amu of lattice atom
% lambda = 401; %Imaging wavelength
% addnoise = 0.005; %Noise intensity (0 to 1)
% collectedphotonsratio =  0.09; %QE of Andor and 20%loss
% 
% %Number of different conditions in which to test network
% num_conditions = 6;
% 
% %SIMULATE IMAGES AND RECONSTRUCT OCCUPATIONS
%     
% %Randomly generated lattice occupation
%  patterns = zeros(sites,sites,num_pics);
%  patterns(3,3,1:num_pics/2) = 1;
% 
%  %Array to hold testing pictures
% test_pics = zeros(sites*pixelspersite, sites*pixelspersite, num_pics);
% 
% %Simulate pictures
% for j = 1:num_pics
%     test_pics(:,:,j) = simulate_setpattern( patterns(:,:,j), pixelspersite, 1,...
%         NA, latticespacing, imagingpulse, lossrate,recoilvel, scatteringrate,...
%         atomicmass, lambda, addnoise, collectedphotonsratio,0);
% end
% 
% %Binarize simulated images
% test_pics(test_pics > 0) = 1;
% 
% %Pad images with a border of zeros
% test_pics = padarray(test_pics, [pixelspersite pixelspersite 0], 0,'both');
% 
% %Generate lattice coordinates
% mask = ones(size(test_pics(:,:,1)));
% mask(1:pixelspersite-1,:) = 0;
% mask(:,1:pixelspersite-1) = 0;
% mask(size(mask,1)-pixelspersite+1:size(mask,1),:) = 0;
% mask(:,size(mask,2)-pixelspersite+1:size(mask,2)) = 0;
% 
% [ lattice_coords, lattice_indices ] = gen_lattice( test_pics(:,:,1),...
%     pixelspersite, pixelspersite, pixelspersite, mask);

%Reconstruct lattice
[binary_1, ~, ~] = reconstruct_lattice(test_pics, lattice_coords, lattice_indices,...
    pixelspersite, 0.5, net_path_1, true, -1);
errs  = binary_1(3,3,:) - patterns(3,3,:);
global_fidelity = 1 - mean(abs(errs));
false_atom = subplus(errs);
false_hole = subplus(-errs);

false_atom_rate = mean(false_atom)*2;
false_hole_rate = mean(false_hole)*2;

intens_bin = sum(test_pics(24:32,24:32,:),[1,2]);
histogram(intens_bin)
test_fids = zeros(40,1);
test_false_atom = zeros(40,1);
test_false_hole = zeros(40,1);
for x = 1:40
    test_class = intens_bin;
    thresh = x;
    test_class(test_class <= thresh) = 0;
    test_class(test_class >= thresh) = 1;
    test_errs  = test_class - patterns(3,3,:);
    test_fids(x) = 1 - mean(abs(test_errs));
    test_false_atom(x) = mean(subplus(test_errs))*2;
    test_false_hole(x) = mean(subplus(-test_errs))*2;
end

plot(test_fids)
hold on
plot(test_false_atom)
plot(test_false_hole)
xlabel('Photon count threshold')
ylabel('Fidelity')
legend('Hole identification fidelity','False atom rate','False hole rate','Global fidelity')
hold off

figure(2)
plot(1 - (0.89*test_false_hole + 0.11*test_false_atom))
xlabel('Photon count threshold')
ylabel('Lattice-wide fidelity')

all_site_errs  = binary_1(2:4,2:4,:) - patterns(2:4,2:4,:);
all_site_fidelity = 1 - mean(mean(mean(abs(all_site_errs))));