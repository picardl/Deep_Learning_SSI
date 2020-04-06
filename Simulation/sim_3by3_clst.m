function sim_3by3_clst(params)
%Simulate a set of 3 by 3 lattices containing all possible combinations of
%neighbours for central site, for neural net training or testing
%External definitions:
%   reps - number of copies of each pattern to be simulated
% 	pixelspersite – width of site in pixels.
%   NA - numerical aperture (between 0 and 1)
%   latticespacing - distance between neighbouring lattice sites, in nm
%   in nm
%   imagingpulse - length of imaging laser pulse, in s
%   scatteringrate - rate at which atoms scatter photons, in s^-1
%   lossrate - probability of atom loss in every scattering event
%   recoilvel - atom recoil velocity from scattering in m/s

%Number of copies of each occupation pattern to be simulated for training
%data
reps = params.reps;

%Set simulation parameters
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
binarize = params.binarize;
latticedepth = params.latticedepth;

%3D matrix to hold training pictures
nn_training_pics = zeros(round(3*pixelspersite), round(3*pixelspersite), reps*512);

%All ossible lattice occupation patterns, reshaped to 3x3 array
combs = dec2base(0:power(2,9) - 1,2) - '0';
nn_patts = reshape(combs', 3, 3, 512);

%Matrix of occupation patterns to be saved, with a pattern matrix for every
%training image
nn_training_patts = zeros(3, 3, reps*512);

mod_filling = ones(5,5,5*reps);
fillings = linspace(0,1,5*reps);
for i = 1:5*reps
    mod_filling(:,:,i) = mod_filling(:,:,i)*fillings(i);
end

%Randomly generated lattice occupation
patterns = round( rand(5,5, 5*reps) + mod_filling - 0.5);
patterns(2:4,2:4,:) = 0;

 %Array to hold testing pictures
test_pics = zeros(round(5*pixelspersite),round(5*pixelspersite), 5*reps);

%Simulate pictures
for j = 1:5*reps
    test_pics(:,:,j) = simulate_setpattern( patterns(:,:,j), pixelspersite, 1,...
        NA, latticespacing, imagingpulse, lossrate,recoilvel, scatteringrate,...
        atomicmass, lambda, 0, collectedphotonsratio,latticedepth);
end

%Binarize simulated images
test_pics(test_pics > 0) = 1;
center_pics = test_pics(round(pixelspersite) +1:round(4*pixelspersite),round(pixelspersite)+1:round(4*pixelspersite),:);

for i = 1:512
    %Generate training pictures
    nn_training_pics(:,:, (i-1)*reps + 1 : i*reps) = simulate_setpattern(nn_patts(:,:,i), pixelspersite, reps, NA, latticespacing, imagingpulse, lossrate, recoilvel, scatteringrate, atomicmass, lambda, addnoise, collectedphotonsratio, latticedepth);
    
    nn_training_patts(:, :, (i-1)*reps + 1 : i*reps) = repmat(nn_patts(:,:,i), 1, 1, reps);
    
    if rem(i,50) == 0
        c = clock;
        fprintf('%i/512 patterns simulated at time %i:%i:%i\n', i, c(4), c(5),round(c(6)));
    end
    
    randBorderInds = randi(5*reps, reps,1);
    randBorders = center_pics(:,:,randBorderInds);
    nn_training_pics(:,:, (i-1)*reps + 1 : i*reps) = nn_training_pics(:,:, (i-1)*reps + 1 : i*reps) + randBorders;
end

%Normalise images so all pixels lie in range [0 1]
nn_training_pics = nn_training_pics/max(max(max(nn_training_pics)));

%Binarize images if switch active, otherwise normalise them
if binarize
    nn_training_pics(nn_training_pics > 0) = 1;
else
    for i=1:reps*512
	im = nn_training_pics(:,:,i);
        mean_bright = mean(mean(im(im>0)));
        nn_training_pics(:,:,i) = tanh(atanh(0.5)*im/mean_bright);
    end
end

%Save data
save(sprintf('%straining_data',params.savefolder), 'nn_training_pics', 'nn_training_patts', '-v7.3')

%Save 5 random simulated pictures
for i = 1:5
    im = nn_training_pics(:,:,randi(reps*512));
    imwrite(im, sprintf('%ssample_image_%i.png', params.savefolder, i))
end

clear nn_training_pics
