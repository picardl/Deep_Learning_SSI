function train_new_convnet(reps, binarize, NA, latticespacing, imagingpulse, lossrate, recoilvel, scatteringrate, pixelspersite, atomicmass, lambda, addnoise, collectedphotonsratio, latticedepth, poolsize, stride, filtersize, numfilters)

rng('shuffle') %Seed random number generator
jobid = datestr(clock, 'ddmmmyy_HHMMSS'); %Job ID as datetime
jobname = sprintf('convnet_%s', jobid);

workdir = pwd; %working directory

%Folder to save output data
params.savefolder = strcat(workdir,'/',jobname);

system(sprintf('mkdir "%s"', params.savefolder));

params.savefolder = strcat(params.savefolder,'/');

%Number of copies of each occupation pattern to be simulated for training
%data
params.reps = reps;

%Set simulation parameters
params.NA = NA; %Numerical aperture
params.latticespacing = latticespacing; %in nm
params.imagingpulse = imagingpulse; %in s
params.lossrate = lossrate; %loss probability per scattering event
params.recoilvel = recoilvel; %scattering recoil velocity in m / s
params.scatteringrate = scatteringrate; %per s
params.pixelspersite = pixelspersite; %Distance between lattice sites in px
params.atomicmass = atomicmass; %Mass in amu of lattice atom
params.lambda = lambda; %Imaging wavelength
params.addnoise = addnoise; %Noise intensity (0 to 1)
params.collectedphotonsratio = collectedphotonsratio; %Fraction of photons which are collected by CCD
params.binarize = binarize; %True if images are to be binarized, false otherwise
params.latticedepth = latticedepth; %Depth of pinning lattice, in units of E_r
params.poolsize = poolsize; %Convolutional pooling size
params.stride = stride; %Convolutional pooling stride
params.filtersize = filtersize; %Convolutional filter size
params.numfilters = numfilters; %Convolutional numfilters

% the parameters structure is done, save it to file
param_filename = sprintf('params_%s.mat', jobname);
param_path = strcat(params.savefolder,param_filename);
save(param_path, 'params');

%Make sure all folders in working directory are in path
addpath(genpath(pwd))

%Change to new working directory
cd(params.savefolder)

%Create log file and start logging output
logfile = sprintf('convnet_log_%s', jobid);
diary(logfile)

disp(params)

%Generate training data
sim_3by3_clst(params)

%Train and save convolutional neural network
convolutional

clear pic_data
diary off