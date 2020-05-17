function [binary_reconstruction, reconstruction_list] = reconstruct_lattice_gaussfit(images, lattice_coords, lattice_indices, PSFwidth, spacing, threshold, binarize)
%Reconstructs a set of optical lattice images using fit to 9 Gaussian PSFs
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
%   binarize: boolean, if true preprocess image by setting all pixels with
%           values greater than 0 to 1. Choice of this input depends on
%           whether network was trained with binary or continuous images.
%   filling: float in range {0 1} corresponding to estimate of filling
%           fraction around site. If not used, set to -1.
%
%Returns:
%   binary_reconstruction: Array containing reconstructed images as a stack
%           of binary arrays.


%Binarize or normalize images, depending on input switch
if binarize
    images(images > 0) = 1;
end

%Number of images to reconstruct
N = size(images, 3);

%Create empty array to hold network outputs
binary_reconstruction = zeros(max(lattice_indices(:,1)), max(lattice_indices(:,2)), N);
reconstruction_list = zeros(size(lattice_coords, 1),N);%Column vectors to hold occupation of each coordinate pair

for i = 1:N %Loop over all images
    for j = 1:size(lattice_coords, 1) %Loop over all sites in each image
        
        %Coordinates of central site
        xcord = lattice_coords(j,1);
        ycord = lattice_coords(j,2);
        
        %Create 3-by-3 site image segment
        im_segment = images(round(ycord - round(3*spacing/2) + 1):round(ycord + round(3*spacing/2)),...
            round(xcord - round(3*spacing)/2 + 1):round(xcord + round(3*spacing)/2), i);

        %Classify central site and store output in array
        [all_amps_out, ~] = GaussiansAmplitudeFit(im_segment,spacing, PSFwidth);
        center_amp = all_amps_out(5);
        if center_amp >= threshold
            site_out = 1;
        elseif center_amp < threshold
            site_out = 0;
        end
        binary_reconstruction(lattice_indices(j, 1), lattice_indices(j, 2), i) = site_out;
        reconstruction_list(j, i) = site_out;
    end
end

