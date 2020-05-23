function [ center_state, neighbours ] = neighbour_combns( input_image, sites, spacing, PSF)
%NEIGHBOUR_PERMS Finds all possible permutations of neighbour occupancies
%for a single uncertain site, returns central site occupancy which
%minimises difference with real picture

center = 0.5 + size(input_image)/2;

%Reshape 3x3 matrix of site states to a vector and determine number of
%known occupancy sites and their positions
lin_sites = reshape(sites, [], 1);
certain = find(lin_sites == 0 | lin_sites == 1);
uncertain_sz = length(lin_sites) - length(certain);

%Generate all possible combinations of filled and empty states for
%uncertain sites
combs = dec2base(0:power(2,uncertain_sz) - 1,2) - '0';

y_width = size(PSF, 1)/2;
x_width = size(PSF, 2)/2;


if size(PSF, 1) < spacing || size(PSF, 2) < spacing
    disp('Error, neighbour_perms can only operate on images where the width of the PSF is greater than the lattice spacing') %If this condition is false, the function should be unecessary anyway
    return
end

if ~isempty(certain)
%Create vectors with all combinations of uncertain and certain sites in
%correct order

        %First fill any leading uncertain sites with all possible
        %combinations
        neighbourhoods(:, 1:certain(1) - 1) = combs(:,1:certain(1) - 1);
        neighbourhoods(:, certain(1)) = lin_sites(certain(1));

        i = 1;

    if length(certain) > 1
        %Insert certain sites in their original positions in every
        %uncertain site combination
        for i = 2:length(certain);
            neighbourhoods(:, certain(i-1) + 1 : certain(i) - 1) = combs(:, certain(i-1) + 2 - i: certain(i) - i);
            neighbourhoods(:, certain(i)) = lin_sites(certain(i));
        end
    end
    
    %Fill any remaining uncertain sites
    neighbourhoods(:, certain(i) + 1 : 9) = combs(:, certain(i) - i + 1: 9 - i);
else
    neighbourhoods = combs; %If all sites are uncertain, simply use all combinations for 9 sites
end

non_z_ind = find(input_image > 0); %Positions of photon detections in original image

N = sum(sum(input_image)); %Total recorded intensity counts

chi_sq = zeros([size(neighbourhoods, 1), 1]); %Holds chi squared values for each possible filling combination

for i = 1:size(neighbourhoods, 1)
        %Add the PSF centred on a particular site determined by vector
        %neighbourhoods
        
    filled_image = zeros(size(input_image)); %Will hold image with PSFs added at filled sites
        
    if neighbourhoods(i, 1) == 1
        cent_pos = floor(center - [spacing, spacing]);
        add_PSF = PSF((y_width - cent_pos(1)):end, (x_width - cent_pos(2)):end);
        add_PSF = padarray(add_PSF, [size(input_image,1)-size(add_PSF,1), 0], 'post');
        add_PSF = padarray(add_PSF, [0, size(input_image,2)-size(add_PSF,2)], 'post');
        filled_image = filled_image + add_PSF;
    end

    if neighbourhoods(i, 2) == 1
        cent_pos = floor(center - [0, spacing]);
        add_PSF = PSF(:, (x_width - cent_pos(2)):end);
        add_PSF = padarray(add_PSF, [(size(input_image,1)-size(add_PSF,1))/2, 0]);
        add_PSF = padarray(add_PSF, [0, size(input_image,2)-size(add_PSF,2)], 'post');
        filled_image = filled_image + add_PSF;
    end

    if neighbourhoods(i, 3) == 1
        cent_pos = floor(center - [-spacing, spacing]);
        add_PSF = PSF(1: end - (cent_pos(1) + y_width - size(input_image,1)), (x_width - cent_pos(2)):end);
        add_PSF = padarray(add_PSF, [size(input_image,1)-size(add_PSF,1), 0], 'pre');
        add_PSF = padarray(add_PSF, [0, size(input_image,2)-size(add_PSF,2)], 'post');
        filled_image = filled_image + add_PSF;
    end

    if neighbourhoods(i, 4) == 1
        cent_pos = floor(center - [spacing, 0]);
        add_PSF = PSF((y_width - cent_pos(1)):end, :);
        add_PSF = padarray(add_PSF, [size(input_image,1)-size(add_PSF,1), 0], 'post');
        add_PSF = padarray(add_PSF, [0, (size(input_image,2)-size(add_PSF,2))/2]);
        filled_image = filled_image + add_PSF;
    end

    if neighbourhoods(i, 5) == 1
        cent_pos = floor(center);
        add_PSF = PSF;
        add_PSF = padarray(add_PSF, [(size(input_image,1)-size(add_PSF,1))/2, (size(input_image,2)-size(add_PSF,2))/2]);
        filled_image = filled_image + add_PSF;
    end

    if neighbourhoods(i, 6) == 1
        cent_pos = floor(center - [-spacing, 0]);
        add_PSF = PSF(1: end - (cent_pos(1) + y_width - size(input_image,1)), :);
        add_PSF = padarray(add_PSF, [size(input_image,1)-size(add_PSF,1), 0], 'pre');
        add_PSF = padarray(add_PSF, [0, (size(input_image,2)-size(add_PSF,2))/2]);
        filled_image = filled_image + add_PSF;
    end

    if neighbourhoods(i, 7) == 1
        cent_pos = floor(center - [spacing, -spacing]);
        add_PSF = PSF((y_width - cent_pos(1)):end, 1: end - (cent_pos(2) + x_width - size(input_image,2)));
        add_PSF = padarray(add_PSF, [size(input_image,1)-size(add_PSF,1), 0], 'post');
        add_PSF = padarray(add_PSF, [0, size(input_image,2)-size(add_PSF,2)], 'pre');
        filled_image = filled_image + add_PSF;
    end

    if neighbourhoods(i, 8) == 1
        cent_pos = floor(center - [0, -spacing]);
        add_PSF = PSF(:, 1: end - (cent_pos(2) + x_width - size(input_image,2)));
        add_PSF = padarray(add_PSF, [(size(input_image,1)-size(add_PSF,1))/2, 0]);
        add_PSF = padarray(add_PSF, [0, size(input_image,2)-size(add_PSF,2)], 'pre');
        filled_image = filled_image + add_PSF;
    end

    if neighbourhoods(i, 9) == 1
        cent_pos = floor(center - [-spacing, -spacing]);
        add_PSF = PSF(1: end - (cent_pos(1) + y_width - size(input_image,1)), 1: end - (cent_pos(2) + x_width - size(input_image,2)));
        add_PSF = padarray(add_PSF, [size(input_image,1)-size(add_PSF,1), 0], 'pre');
        add_PSF = padarray(add_PSF, [0, size(input_image,2)-size(add_PSF,2)], 'pre');
        filled_image = filled_image + add_PSF;
    end
    
    %Normalise integral (approximated by sum of whole image) to 1
    
    filled_image = filled_image / max(max(filled_image));
    filled_image = filled_image + 0.005;
    
    norm_img = input_image / max(max(input_image));
    
    %Calculate chi squared statistic for each possible neighbourhood
    chi_sq(i) = sum(((norm_img(non_z_ind) - filled_image(non_z_ind)).^2)./(filled_image(non_z_ind)));
    
end

%Select site with minimum chi squared statistic
[~, min_ind] = min(chi_sq);
neighbours = reshape(neighbourhoods(min_ind, :), 3, 3);

center_state = neighbours(2, 2);

end

