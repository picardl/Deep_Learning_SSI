function [ lattice_sites, lattice_indices ] = gen_lattice( image, spacing, x_off, y_off, mask)
%Generates an array of lattice site center coordinates.
    %Based on the dimensions of a sample image, the function determines the
    %centers of lattice sites for a particular spacing and translational
    %alignment
%     
%   Inputs:	image – The image for which the lattice coordinates are generated. May be a single image sample from a set, as it is only used for dimensions.
%           spacing - lattice spacing in pixels
%           x_off –  x offset, output from rotate_and_align
%           y_off – y offset, output from rotate_and_align
%           mask – mask, output from rotate_and_align
% 
% 	Outputs:lattice_sites – coordinates of all lattice sites in image, in [x y] format.
%           lattice_indices – indices of pattern matrix corresponding to each pair of lattice coordinates.


%Expected number of lattice sites in image in x and y directions
x_site_num = round((size(image, 2)-2*x_off)/spacing);
y_site_num = round((size(image, 1)- 2*y_off)/spacing);

%Vectors to hold lattice coordinates in image and corresponding indices in pattern matrix
lattice_sites = zeros([x_site_num*y_site_num, 2]);
lattice_indices = zeros([x_site_num*y_site_num, 2]);

k=1;

for i = 1:x_site_num
    for j = 1:y_site_num

		%Do not generate coordinates for any lattice site that would have more than 25% of area outside image
        if mask(round((i - 0.5) * spacing + x_off), round((j - 0.5) * spacing + y_off)) == 0
            continue
        end
        
        if ((i + 0.5) * spacing + x_off > size(image, 2)) || ((j - 0.5) * spacing + y_off > size(image, 1))
            continue
        end
            
        
		%Generate lattice site coordinate
        lattice_sites(k, :) = [(i - 0.5) * spacing + x_off, (j - 0.5) * spacing + y_off];
        lattice_indices(k, :) = [j, i];
        
        k = k + 1;
    end
end

lattice_sites = lattice_sites(1 : k - 1, :);
lattice_indices = lattice_indices(1 : k - 1, :);

end

