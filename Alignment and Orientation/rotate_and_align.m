function [ crop_pics, x_offset, y_offset, mask ] = rotate_and_align( input_pics, angle, spacing)
%ROTATE_AND_ALIGN Takes a set of input images with known lattice angle and
%spacing, rotates them to align the lattice and image axes, and crops out
%any regions of the image outside the bounds of the lattice.
% 
% 	Inputs:	input_pics – set of pictures in 3D matrix
% 		angle – lattice angle relative to image axes
% 		spacing - lattice spacing in pixels
% 
% 	Outputs:crop_pics – input images rotated and cropped
% 		x_offset – offset, in pixels, between left edge of image and left edge of lattice.
% 		y_offset – offset, in pixels, between left edge of image and left edge of lattice.
% 		mask – mask, used in generation of lattice coordinates, has the same size as the 
%         input pictures and a value of 0 at any position representing a potential lattice site 
%         which is not fully in the bounds of the image.


%Create a mask of ones with a border of zeros which will later be used to
%filter out lattice sites only partially captured by the camera
mask = ones(size(input_pics(:,:,1)));
border = round(spacing/4);
mask([1:border, end - border : end], :) = 0;
mask(:, [ 1 : border, end - border : end] ) = 0;

%Rotate image to align lattice and image axes
rot_pics = imrotate(input_pics, -angle);
mask = imrotate(mask, -angle);

%Add all pictures together and project along x- and y-axes
xproj_sum = sum( sum(  rot_pics , 3 ), 1 );
yproj_sum = sum( sum(  rot_pics , 3 ), 2 )';

%Determine region of image actually occupied by lattice, using an estimate
%of the signal-to-noise ratio.
SNR = 5;
[left, right, top, bottom] = lattice_bounds(xproj_sum, yproj_sum, SNR);

%Crop images to only retain area occupied by lattice
crop_pics = rot_pics(top:bottom, left:right, :);
mask = mask(top:bottom, left:right);

%Reevaluate projections for cropped images
xproj_sum_crop = sum( sum(  crop_pics , 3 ), 1 );
yproj_sum_crop = sum( sum(  crop_pics , 3 ), 2 )';

%Determine horizontal and vertical alignment of lattice relative to image
%by fitting series of equally spaced gaussians
x_offset = alignment_fit(xproj_sum_crop, spacing);
y_offset = alignment_fit(yproj_sum_crop, spacing);

end

