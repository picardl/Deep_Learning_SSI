% Note whole Imaging Simulation folder must be included in path to run,
% except for the original version of CCDDetection_singlepicture

sites = 10; %simulated pattern has dimensions sites x sites
pixels_per_site = 20;
fraction_filled = 0.2;
init_angle = 30; %in degrees
n_pics = 1000;
NA = 0.85;

[sim_pics, sim_patterns, PSFwidth] = simulate_npictures(sites, pixels_per_site, fraction_filled, init_angle, NA, n_pics)  ; %3 dimensional array containing simulated pictures, associated lattice patterns and PSF width to be used later

[angle, spacing] = find_angle_and_spacing( sim_pics, sites, pixels_per_site, fraction_filled, init_angle);


%Create a mask of ones with a border of zeros which will later be used to
%filter out lattice sites only partially captured by the camera
mask = ones(size(sim_pics(:,:,1)));
border = round(spacing/4);
mask([1:border, end - border : end], :) = 0;
mask(:, [ 1 : border, end - border : end] ) = 0;

%Rotate image to align lattice and image axes
rot_pics = imrotate(sim_pics, -angle);
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

%Generate lattice site coordinates
[ lattice_coords, lattice_indices ] = gen_lattice( crop_pics(:,:,1), spacing, x_offset, y_offset, mask);