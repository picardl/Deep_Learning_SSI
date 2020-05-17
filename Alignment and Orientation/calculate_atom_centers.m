function [ all_centers, midpoint ] = calculate_atom_centers( all_pictures, sites, pixel_spacing, fraction_filled )
%Returns array containing center coordinates of identifiable atoms
%   Takes a 3 dimensional array containing simulated images as input, returns an
%   array containing x and y coordinates of all clearly identifiable
%   isolated atoms. The three scalar arguments are approximate values of
%   the number of sites per side, spacing in pixels and filling fraction of the
%   lattice.

if isempty(sites), sites = 10; end
if isempty(pixel_spacing), pixel_spacing = 20; end
if isempty(fraction_filled), fraction_filled = 0.5; end

all_centers = [;];

input_size = size( all_pictures );

input_length = input_size(3);

midpoint = ( double( input_size(1:2) ) / 2) + 0.5;

%Define square neighbourhood around maximum, with sides 3 times lattice
%spacing, to ignore any atoms in neighbouring sites
neighbourhood_boundary = 3 * pixel_spacing;
if mod(neighbourhood_boundary,2) == 0
    neighbourhood_boundary = neighbourhood_boundary + 1; %square centred on single point must have odd integer side length
end

find_atoms = vision.LocalMaximaFinder; %Local maxima finder object
find_atoms.MaximumNumLocalMaxima = fraction_filled*sites^2; %Find at most the expected number of atoms
find_atoms.NeighborhoodSize = [neighbourhood_boundary neighbourhood_boundary]; %Suppress other maxima within neighbourhood of detected maximum
find_atoms.IndexDataType = 'double';


for i = 1:input_length
    
    %Apply Gaussian blur to approximate broad PSF from isolated detection peaks
    %This c
    h = fspecial('gaussian',10,5);
    blur = imfilter( all_pictures(:, :, i), h );
    
    max_value = max( max( blur ) );
    
    %Filter out low intensity background in order to separate peaks
    minus_background = blur > 0.3 * max_value;
    blur = blur .* minus_background;
    
    %Filter out objects in binarized image with areas below estimated noise
    %threshold and above max area of isolated atom, assuming 20 pixels per site
    binary = im2bw(blur, 0.1);
    singles = bwareafilt(binary, [0.25*pixel_spacing^2, 1.25*pixel_spacing^2] );
    blur = blur.*singles;
    
    find_atoms.Threshold = round( 0.4 * max_value ); %Ignore any peaks below 40% of tallest peak height
    centers = step( find_atoms, blur );
    centers_size = size( centers );
    
    adjust = 0;
    
    %Test whether detected centers are actually maxima, since zeroing of
    %neigbourhood around detected maxima can lead to later false detections
    if ~isempty( centers )
        for j = 1:centers_size(1)
            
            if j + adjust > size(centers, 1)
                break
            end
            
            if not_max(blur, centers(j + adjust, :))
                centers(j + adjust, :) = [];
                adjust = adjust - 1;
            end
        end
    end
     
%     figure(1)
%     surf(blur)
% 
%     figure(2)
%     image(63*blur/max_value)
%     hold on
%     plot(centers(:,1),centers(:,2),'r.');
%     hold off
%     
%     figure(3)
%     image(binary)
    
    all_centers = [all_centers ; centers]; %#ok<AGROW>
end

end

