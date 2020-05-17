%This script is used to test the various lattice alignment processes
%required to process real lattice images.
%Simulate a set of lattice pictures with given paramters, find angle and
%spacing, rotate and align and generate lattice coordinates.

sites = 10; %simulated pattern has dimensions sites x sites
fraction_filled = 0.3;
init_angle = 22; %in degrees
n_pics = 50;
NA = 0.85; %Numerical aperture
latticespacing = 256; %in nm
imagingpulse = 3e-6; %in s
lossrate = 1e-4; %loss probability per scattering event
recoilvel = 5.9e-3; %scattering recoil velocity in m / s
scatteringrate = 2*pi*30e6; %per s
pixelspersite = 20; %Distance between lattice sites in px
atomicmass = 168; %Mass in amu of lattice atom
lambda = 401; %Imaging wavelength
AddNoise = 0.005; %Noise intensity (0 to 1)

angle = 20; %Angle of rotation of lattice with respect to image

%3 dimensional array containing simulated pictures, associated lattice patterns and PSF
%[sim_pics, sim_patterns, PSFwidth] = simulate_npictures(sites, fraction_filled, angle, pixelspersite, n_pics, NA, latticespacing, imagingpulse, lossrate, recoilvel, scatteringrate, atomicmass, lambda, AddNoise);

%[angle, spacing] = find_angle_and_spacing(sim_pics, sites, pixelspersite, latticespacing, fraction_filled, init_angle);

[ crop_pics, x_offset, y_offset, mask ] = rotate_and_align(sim_pics, angle, spacing);

[ lattice_coords, lattice_indices ] = gen_lattice( crop_pics(:,:,1), spacing, x_offset, y_offset, mask); %Generate lattice site coordinates