function [ resampled_images ] = square_resample( images, compression_factor, lattice_separation )
%RESAMPLE Resample a set of images by averaging nearby pixels
%   
%   Inputs:     images - i*j*k matrix, where each k is a separate image and
%               i and j are coordinates in a 2D image
%               compression_factor - factor by which to compress image (in
%               one dimension)
%               lattice_separation - integer number of pixels separating
%               neighbouring lattice sites
% 
%   Outputs:    resampled_images - 3D matrix of resampled images

num_pics = size(images, 3);

%Calculate new pixel separation, which must be an integer factor of input
resamp_separation = ceil(lattice_separation/compression_factor);
round_factor = lattice_separation/resamp_separation;

%Side length, in pixels, of resampled square image
resamp_length = size(images, 1)*resamp_separation/lattice_separation;

resampled_images = zeros(resamp_length, resamp_length, num_pics);

for k = 1:num_pics
    for i = 1:resamp_length
        for j = 1:resamp_length
            resampled_images(i, j, k) = ...
                mean(mean(images((i - 1)*round_factor + 1 : i*round_factor, ...
                (j - 1)*round_factor + 1 : j*round_factor, k)));
        end
    end
end

end
