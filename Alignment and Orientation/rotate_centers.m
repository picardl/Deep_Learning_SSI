function [FT_max, max_spacing, h1] = rotate_centers(angle, centers, midpoint, pixelspersite_guess, spacing_guess)
%ROTATE_centers Computes new positions of atom centers upon rotation
%   Rotates centers by given angle, calculates mutual distances along x
%   axis of original and rotated centers, returns negative of maximum value of Fourier
%   Transform of mutual distance histogram and its index

%Transform coordinate system to have origin at center of picture
centers(:,1) = centers(:, 1) - midpoint(2);
centers(:,2) = midpoint(1) - centers(:, 2);

%Rotate centers about midpoint
radii = sqrt(centers(:, 1) .^ 2 + centers(:, 2) .^ 2);
init_angles = atan2(centers(:, 2),centers(:, 1));
new_angles = init_angles + angle;
rot_centers = [radii .* cos(new_angles), radii .* sin(new_angles)];

%Define mutual distance array
cent_length = length( rot_centers );
mut_length = cent_length * ( cent_length + 1 ) / 2;
rot_mut_dist = zeros( [mut_length 1] );
k = 1;

%Calculate mutual distances in x-direction
for i = 1:( cent_length - 1 )
    for j = (i + 1):cent_length
        rot_mut_dist(k) = abs( rot_centers(j,1) - rot_centers(i,1) );
        k = k + 1;
    end
end

%Bin mutual distances
edges = 0:1:10*pixelspersite_guess;
[h1, ~] = histcounts(rot_mut_dist, edges);

figure(2)
histogram(rot_mut_dist, edges);
xlabel('Mutual distance between centers / px')
ylabel('Frequency of occurence')

FT1 = abs( fft( detrend(h1) ) ) .^ 2;
FT1_real = FT1( 1 : round(length(FT1) / 2));

FT_x = 1:1:length(FT1_real);
wavelengths = spacing_guess * size(edges, 2) ./ ((FT_x - 1) * pixelspersite_guess);

figure(3)
plot(wavelengths, FT1_real)
xlabel('Lattice spacing / nm')
ylabel('Fourier coefficient')

%Maximise Fourier coefficient in a region around expected periodicity of
%lattice
highlamba_ind = round(length(FT1)/ (0.5*pixelspersite_guess));
lowlamba_ind = round(length(FT1)/ (2.5*pixelspersite_guess));
[FT_max, max_index] = max( FT1_real(lowlamba_ind:highlamba_ind) );
max_spacing = length(FT1) / (max_index - 2 + lowlamba_ind); %lattice spacing in pixels
%NB: -2 is due to MATLAB indexing from 1 not 0, but starting from lower
%index inclusive when slicing
FT_max = -1 * FT_max;


end

