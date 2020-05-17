function [ int_avg ] = mean_int( images, lattice_coords, lattice_spacing )
%MEAN_INT Summary of this function goes here
%   Detailed explanation goes here
intensity = zeros(size(images, 3)*length(lattice_coords), 1);

%Integer lattice coordinates and spacing
int_coords = round(lattice_coords);
int_sp = round(lattice_spacing);

%Bounds of each lattice site, with bounds outside range of image rounded to
%boundary
x_bounds = [int_coords(:,1) - int_sp/2, int_coords(:,1) + int_sp/2];
x_bounds(x_bounds < 1) = 1;
x_bounds(x_bounds > size(images,2)) = size(images,2);

y_bounds = [int_coords(:,2) - int_sp/2, int_coords(:,2) + int_sp/2];
y_bounds(y_bounds < 1) = 1;
y_bounds(y_bounds > size(images,1)) = size(images,1);

%Finds intensity at each lattice site, weighted by the Gaussian PSF
for j = 1:size(images, 3)
    for i = 1:size(lattice_coords,1)
            site_isol = images(y_bounds(i,1) : y_bounds(i,2), x_bounds(i,1) : x_bounds(i,2));
            site_isol = reshape(site_isol, size(site_isol,1)*size(site_isol, 2), 1);
            site_isol = site_isol(site_isol > 0);
            intensity((j-1)*length(int_coords) + i) = sum(site_isol) / length(site_isol);
    end
    
end
nbins = 40;
histogram(intensity, nbins)

int_avg = mean(intensity);


