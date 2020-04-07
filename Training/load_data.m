%Loads training data from saved files
%Unpacks simulation paremeters struct 'params' which must be loaded into
%instance

if exist('params','var') == 0
    error('Must load params struct before calling load_data.m')
end

%Number of copies of each occupation pattern in training data
reps = params.reps;
N = reps*512;

%Load simulation parameters
NA = params.NA; %Numerical aperture
latticespacing = params.latticespacing; %in nm
imagingpulse = params.imagingpulse; %in s
lossrate = params.lossrate; %loss probability per scattering event
recoilvel = params.recoilvel; %scattering recoil velocity in m / s
scatteringrate = params.scatteringrate; %per s
pixelspersite = params.pixelspersite; %Distance between lattice sites in px
atomicmass = params.atomicmass; %Mass in amu of lattice atom
lambda = params.lambda; %Imaging wavelength
poolsize = params.poolsize; %Convolutional pooling size
stride = params.stride; %Convolutional pooling stride
filtersize = params.filtersize; %Convolutional filter size
numfilters = params.numfilters; %Convolutional numfilters

%Open training saved training data file
savedata = matfile(sprintf('training_data.mat'));

%If training data not already loaded, load and vectorise all images
if exist('pic_data', 'var') == 0
    
    N_pics = int32(size(savedata.nn_training_pics, 3));
    dim1 = size(savedata.nn_training_pics, 1);
    dim2 = size(savedata.nn_training_pics, 2);
    pic_data = zeros(dim1*dim2, N_pics);
    square_pics = zeros(dim1, dim2,1, N_pics);
    
    %Number of images to load from training data in each batch
    %Needed to minimise memory usage
    loadsize = 1000;
    N_batches = int32(ceil(N_pics/loadsize));

    %Load training data in batches
    for i = 1:N_batches
        endind = int32(min(i*loadsize, N_pics));
        pic_data(:,(i - 1)*loadsize + 1 : endind) = reshape(savedata.nn_training_pics(:,:,(i - 1)*loadsize + 1 : endind), [dim1*dim2, int32(min(loadsize, N_pics - (i-1)*loadsize))]);
        square_pics(:,:,1, (i - 1)*loadsize + 1 : endind) = savedata.nn_training_pics(:,:,(i - 1)*loadsize + 1 : endind);
    end
end

nn_patts = reshape(savedata.nn_training_patts(2,2,:), 1, []);