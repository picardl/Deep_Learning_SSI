function create_nn_job_integrated(reps, binarize, maxepoch, numopen, rt, mem, NA, latticespacing, imagingpulse, lossrate, recoilvel, scatteringrate, pixelspersite, atomicmass, lambda, addnoise, collectedphotonsratio, latticedepth, poolsize, stride, filtersize, numfilters)

rng('shuffle') %Seed random number generator
jobid = randi([1,10000],1,1); %random job id to prevent output overwrite

jobname = sprintf('nn_integrated_id%i', jobid);

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

%If first layer output size default bigger than input image, shrink all
%layers by a common factor
shrink = 1;
if (3*pixelspersite)^2 < 1000
    shrink = ceil(1000/(3*pixelspersite)^2);
end

%Set autoencoder parameters
params.maxepoch = maxepoch; %Max number of epochs for pretraining
params.numhid = round(1000/shrink); %Number of output neurons of each layer
params.numpen = round(500/shrink); 
params.numpen2 = round(250/shrink); 
params.numopen = numopen; %Number of visible units
params.fine_maxepoch = 100; %Max epochs for backprop finetuning
params.batchsize = 100; %Number of pictures per batch for finetuning

%Maximum expected cross entropy of classifier
%If cross entropy exceeds this, classifier assumes local minimum, restarts
%training (up to 10 times)
params.min_perf = 0.5; 

% the parameters structure is done, save it to file
param_filename = sprintf('params_%s.mat', jobname);
param_path = strcat(params.savefolder,param_filename);
save(param_path, 'params');

% it remains to create the queuescript
filename = sprintf('qs_%s', jobname);
fid = fopen(sprintf('%s', filename), 'w');

% write commands for queue and resources etc.
fprintf(fid, '#!/bin/sh\n');
fprintf(fid, '#$ -q std.q\n');
fprintf(fid, '#$ -N %s\n', jobname);
fprintf(fid, '#$ -j yes\n');
fprintf(fid, '#$ -cwd\n');
fprintf(fid, sprintf('#$ -l h_rt=%i:00:00\n', rt));
fprintf(fid, sprintf('#$ -l h_vmem=%.1fg\n', mem));
fprintf(fid, '#$ -M lewispicard@gmail.com\n');
fprintf(fid, '#$ -m e\n');
fprintf(fid, sprintf('#$ -o %s%s.o\n', params.savefolder, jobname));

% avoid MCR cache collision (see http://undocumentedmatlab.com/blog/speeding-up-compiled-apps-startup)
fprintf(fid, 'export MCR_CACHE_ROOT=/tmp/mcr$RANDOM.$USER\n');

fprintf(fid, 'module load matlab/R2018b\n');

% use multiple CPUs if we can
fprintf(fid, '#$ -pe openmp 8\n');
fprintf(fid, 'export OMP_NUM_THREADS=$NSLOTS\n');

% call the executable!
fprintf(fid, 'cd $HOME/scratch/lewis_SSI\n');
fprintf(fid, sprintf('%s/train_new_net_integrated "%s" \n', workdir, param_path));

% clean up the MCR cache collision fix
fprintf(fid, 'rm -rf $MCR_CACHE_ROOT\n');

% finish the script file (close it)
fclose(fid);
    
% submit the job!
system(sprintf('qsub %s', filename));