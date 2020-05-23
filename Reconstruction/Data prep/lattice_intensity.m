function [ intensity ] = lattice_intensity( image, lattice_coords, lattice_spacing, PSFwidth, deconv, PSF)
%   Calculates the total intensity in each lattice site of an input image

intensity = zeros([size(lattice_coords,1), 1]);

%Integer lattice coordinates and spacing
int_coords = round(lattice_coords);
int_sp = round(lattice_spacing);

%Bounds of each lattice site, with bounds outside range of image rounded to
%boundary
x_bounds = [int_coords(:,1) - int_sp/2, int_coords(:,1) + int_sp/2];
x_bounds(x_bounds < 1) = 1;
x_bounds(x_bounds > size(image,2)) = size(image,2);

y_bounds = [int_coords(:,2) - int_sp/2, int_coords(:,2) + int_sp/2];
y_bounds(y_bounds < 1) = 1;
y_bounds(y_bounds > size(image,1)) = size(image,1);


[X,Y] = meshgrid(1:size(image,2), 1:size(image, 1)); %Coordinates of each pixel

if deconv
    d_im = deconvlucy(image, PSF);
end

%Finds intensity at each lattice site, weighted by the Gaussian PSF
for i = 1:size(lattice_coords,1)
    if deconv
        intensity(i) = sum(sum(d_im(y_bounds(i,1) : y_bounds(i,2), x_bounds(i,1) : x_bounds(i,2))));   
    else
        weight = exp(-(X - int_coords(i, 1)) .^2 / (2 * PSFwidth ^ 2) - (Y - int_coords(i, 2)) .^2 / (2 * PSFwidth ^2));
        w_im = image .* weight;
        intensity(i) = sum(sum(w_im(y_bounds(i,1) : y_bounds(i,2), x_bounds(i,1) : x_bounds(i,2))));
    end
end


end

