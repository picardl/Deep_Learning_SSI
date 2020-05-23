function [binary_reconstruction, sigmoid_reconstruction, confidence_reconstruction, reconstruction_list] = reconstruct_lattice(images, lattice_coords, lattice_indices, spacing, min_confidence, neural_net_path, binarize, filling)
%Reconstructs a set of optical lattice images using a pretrained neural
%network
%
%Arguments
%   images: Array containing images to be reconstructed
%           Dimensions 1 and 2 are image axes, dimension 3 corresponds to
%           distinct images stored in same array
%           For good reconstruction, lattice should not drift between
%           images.
%           Lattice axes should be aligned with image axes. If this is not
%           already the case, images should be rotated before
%           reconstruction.
%   lattice_coords: N-by-2 array containing coordinates of lattice site
%           centres in images, in [x y] (i.e. [col row]) format
%   lattice_indices: N-by-2 array containing unique index of lattice sites
%           corresponding to coordinates in the same row of lattice_coords
%           in [row col] format.
%   spacing: Lattice spacing, in pixels
%   neural_net_path: string, path to a saved matlab neural network object
%   min_confidence: minimum confidence required for assignment of each site
%           in range [0.5 1].
%   binarize: boolean, if true preprocess image by setting all pixels with
%           values greater than 0 to 1. Choice of this input depends on
%           whether network was trained with binary or continuous images.
%   filling: float in range {0 1} corresponding to estimate of filling
%           fraction around site. If not used, set to -1.
%
%Returns:
%   binary_reconstruction: Array containing reconstructed images as a stack
%           of 2D binary arrays.
%   sigmoid_reconstruction: Array containing reconstructed images as raw
%           neural network outputs in the range [0 1]
%   confidence_reconstruction: Array containing reconstructed images with
%           site values of 0, 0.25, 0.75 or 1, corresponding to confidently unoccupied,
%           uncertain unoccupied, uncertain occupied and confidently occupied,
%           respectively
%   reconstruction_list: 2D array containing reconstruction outputs, where
%       first dimension is index in lattice coordinate vector and second
%       dimension is index of input image.

%#function network 

try
    load(neural_net_path, 'net');
catch
    warning('Failed to load net object, please check path contains a MATLAB neural network object named net');
end

%Get input size from either convolutional or network object
if class(net) == "network"
    inputsize = net.input.size;
elseif class(net) == "SeriesNetwork"
    inputsize = prod(net.Layers(1).InputSize);
end

%Check that provided lattice coordinates match training parameters of network
if (inputsize > (round(3*spacing))^2 + 1)||(inputsize < (round(3*spacing))^2)
    disp('Network input is not same size as supplied images. Make sure network training parameters match images');
    return
end

%Keep confidence within range
if min_confidence < 0.5
    min_confidence = 0.5;
    disp('Confidence should be in range [0.5 1]. Rounding to 0.5')
elseif min_confidence > 1
    min_confidence = 1;
    disp('Confidence should be in range [0.5 1]. Rounding to 0.99')
end

%Binarize or normalize images, depending on input switch (Normalization
%comes later)
if binarize
    images(images > 0) = 1;
end

%Number of images to reconstruct
N = size(images, 3);

%Create empty array to hold network outputs
sigmoid_reconstruction = zeros(max(lattice_indices(:,1)), max(lattice_indices(:,2)), N);
reconstruction_list = zeros(size(lattice_coords, 1),N);%Column vectors to hold occupation of each coordinate pair

for i = 1:N %Loop over all images
    for j = 1:size(lattice_coords, 1) %Loop over all sites in each image
        
        %Coordinates of central site
        xcord = lattice_coords(j,1);
        ycord = lattice_coords(j,2);
        
        %Create 3-by-3 site image segment
        try
            im_segment = images(round(ycord - 3*spacing/2 + 1):round(ycord + 3*spacing/2),...
                round(xcord -  3*spacing/2 + 1):round(xcord + 3*spacing/2), i);
        catch
            error('Unable to slice image segment. Ensure that for all lattice coordinates there is a border of at least 1.5 lattice spacings before the edge of the image')
        end
        
        %Reshape segment to vector for input to network
        im_vec = reshape(im_segment, [], 1);
        if filling >= 0
            im_vec = [im_vec; filling];
        end
        
        if size(im_vec) ~= inputsize
            error('Wrong input size');
        end
        
        %Normalise image segment
        if ~binarize
            mean_bright = mean(mean(im_vec(im_vec>0)));
            im_vec = tanh(atanh(0.5)*im_vec/mean_bright);
        end
        
        %Classify central site and store output in array
        if class(net) == "network"
            site_out = net(im_vec);
        elseif class(net) == "SeriesNetwork"
            site_out = grp2idx(classify(net, im_segment)) - 1;
        else
            disp('Network type not recognized')
            site_out = net(im_vec);
        end
        sigmoid_reconstruction(lattice_indices(j, 1), lattice_indices(j, 2), i) = site_out;
        reconstruction_list(j, i) = site_out;
    end
end

%Binarize sigmoid outputs
binary_reconstruction = round(sigmoid_reconstruction);

%Create reconstruction with four categories corresponding to occupied,
%unoccupied, uncertainly occupied and uncertainly unoccupied
confidence_reconstruction = sigmoid_reconstruction;
confidence_reconstruction(confidence_reconstruction < (1 - min_confidence)) = 0;
confidence_reconstruction(confidence_reconstruction < (1 - min_confidence) & ...
    confidence_reconstruction < 0.5) = 0.25; 
confidence_reconstruction(confidence_reconstruction > 0.5 & ...
    confidence_reconstruction < min_confidence) = 0.75;
confidence_reconstruction(confidence_reconstruction > min_confidence) = 1;