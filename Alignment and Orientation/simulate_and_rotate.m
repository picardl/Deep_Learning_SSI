 function [ all_pictures_rot, all_patterns, all_pictures_nonrot ] = simulate_and_rotate(sitelength, n, angle, pixelspersite, NA, latticespacing, imagingpulse, lossrate, recoilvel, scatteringrate, atomicmass, lambda, addnoise, collectedphotonsratio, latticedepth)
%Simulates n square lattice images for a given pattern of occupied / unoccupied
%sites and then rotates them by a given angle, returns images in a 3-dimensional array of height n
% 
% 	Inputs: sitelength - edge length of square lattice to simulate, in
%               number of sites
%           n - number of images to simulate
%           angle - angle (in deg) through which to rotate images after
%               simulation
% 			pixelspersite – width of site in pixels.
% 			n – number of pictures to simulate.
%           NA - numerical aperture (between 0 and 1)
%           latticespacing - distance between neighbouring lattice sites,
%           in nm
%           imagingpulse - length of imaging laser pulse, in s
%           scatteringrate - rate at which atoms scatter photons, in s^-1
%           lossrate - probability of atom loss in every scattering event
%           recoilvel - atom recoil velocity from scattering in m/s
%           collectedphotonsratio - fraction of scattered photons which are
%               collected on CCD
%           latticedepth - Depth of pinning lattice, in Hz
% 
% 	Outputs:
%           all_pictures_rot – 3D matrix, containing all simulated pictures,
%               rotated by specified angle.
%           all_patterns – 3D matrix of binary patterns from which pictures
%               were simulated
%           all_pictures_nonrot – 3D matrix, containing all unrotated
%               pictures
all_pictures_nonrot = zeros(sitelength*pixelspersite,sitelength*pixelspersite,n);
all_patterns = zeros(sitelength,sitelength,n);

    for i = 1:n
        %Generate random pattern and simulate
        all_patterns(:,:,i) = randi(2,sitelength,sitelength) - 1;
        all_pictures_nonrot(:,:,i) = simulate_setpattern( all_patterns(:,:,i), pixelspersite, 1, NA, latticespacing, imagingpulse, lossrate, recoilvel, scatteringrate, atomicmass, lambda, addnoise, collectedphotonsratio, latticedepth);
        
        %Rotate through angle specified leaving dark border around edges
        rotpic = imrotate(all_pictures_nonrot(:,:,i),angle,'loose'); 
        
        %Stack rotated pictures in matrix
        if i == 1
            all_pictures_rot = rotpic;
        else
            all_pictures_rot = cat(3,all_pictures_rot,rotpic);
        end
    end
end