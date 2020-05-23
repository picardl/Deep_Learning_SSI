%This script generates a new empirical point spread function (PSF) for an
%isolated atom, or loads one if a corresponding .mat file exists

if exist('PSF', 'var') == 0 %If PSF already loaded in instance, end script
    
    %If not, check if saved PSF data exists
    if isfile('pointspreadfunction.mat') == 1
        load('pointspreadfunction.mat');
        
        %Reshape normalised PSF to row vector
        PSF_row = reshape(PSF, 1, []);
        
        %Find width of PSF fitted by a Gaussian
        PSFcoefs = gauss_fit2D(PSF, size(PSF,1)/2, size(PSF,2)/2, pixelspersite);
        PSFwidth = PSFcoefs(2);
        
    %If no saved PSF, generate one
    else
        %Generate 10000 sample images of central isolated atoms
        PSF_pattern = zeros(3,3);
        PSF_pattern(2,2) = 1;
        PSF_data = simulate_setpattern(PSF_pattern, pixelspersite, 10000, NA, latticespacing, imagingpulse, lossrate, recoilvel, scatteringrate, atomicmass, lambda, 0, 0.99, 0);

        %Sum images to get a single 2D image
        PSF = sum(PSF_data, 3);
        clear PSF_data 

        %Normalise PSF
        PSF_int = sum(sum(PSF));
        PSF = PSF/PSF_int;
        
        %Fit 2D Gaussian to empirical PSF
        PSFcoefs = gauss_fit2D(PSF, size(PSF,1)/2, size(PSF,2)/2, pixelspersite);
        PSFwidth = PSFcoefs(2);

        %Reshape normalised PSF to row vector
        PSF_row = reshape(PSF, 1, []);

        save pointspreadfunction PSF
    end
end