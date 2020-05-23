function [ deg_output, output_spacing ] = find_angle_and_spacing( calibration_pictures, sites_guess, pixelspersite_guess, spacing_guess, filling_guess, angle_guess )
%Given a set of CCD pictures of a lattice stored as slices of a
%3-dimensional array, this function returns the lattice angle and spacing
%(in pixels).
%   sites_guess : expected number of sites per side of lattice
%   pixelspersite_guess : expected lattice spacing in pixels
%   spacing_guess : expected lattice spacing in nm
%   filling_guess : expected fraction of sites filled
%   angle_guess : expected lattice angle, in degrees

[centers_store, midpoint] = calculate_atom_centers(calibration_pictures, sites_guess, pixelspersite_guess, filling_guess);

rad_angle = degtorad(angle_guess);

%Find lattice angle by maximising lattice frequency component of Fourier
%transform of mutual distance histogram
f = @(angle)rotate_centers(angle, centers_store, midpoint, pixelspersite_guess, spacing_guess); %anonymous function for optimisation
%options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf);
output_angle = -patternsearch(f, -rad_angle, [], [], [], [], -rad_angle - 0.3, -rad_angle + 0.3); %Optimise in 0.6 radian region around expected angle
[~, pixelspersite_guess, histogram] = rotate_centers(-output_angle, centers_store, midpoint, pixelspersite_guess, spacing_guess);

deg_output = radtodeg(output_angle);

%Determine lattice spacing to greater precision than FT output
output_spacing = find_spacing( histogram, pixelspersite_guess );

end
