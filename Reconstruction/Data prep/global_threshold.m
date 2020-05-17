function [ threshold] = global_threshold( input_images, lattice_coords, spacing, PSFwidth, deconv, PSF )
%global_threshold determines a single threshold intensity for naive reconstruction
%from a set of input images.
% 
% 	Inputs: 	input_images - 3D matrix of lattice images.
%               lattice_coords - coordinates of all lattice sites in images, in [x y] format.
%               spacing - lattice spacing in pixels.
%               bound_confidence - desired confidence level (between 0 and 1) for which upper and lower bounds should be determined.
%               deconv - boolean, true = use deconvolution, false = do not use deconvolution
%               PSF - 2D matrix containing a normalised point spread function
% 
% 	Outputs:	threshold - intensity threshold for naive assignment.

intensity_store = zeros(size(lattice_coords, 1)*size(input_images,3), 1);

%Find intensty of every site in image set and store as column vector
for index = 1:size(input_images,3)
    intensity_store((index-1)*size(lattice_coords, 1) + 1 : index*size(lattice_coords, 1)) = lattice_intensity(input_images(:,:,index), lattice_coords, spacing, PSFwidth, deconv, PSF);
end

edges = linspace(0, max(intensity_store), 50);

intensity_hist = histcounts(intensity_store, edges);
x_data = edges(1:end-1) + (edges(2)-edges(1))/2;

% % Get user input for fitting bounds
% disp('Enter bounds for positions of unoccupied and occupied peaks.');
% disp('Enter 999 to use all default values.');
% unocc_lb = input('Enter unoccupied peak lower bound: ');
unocc_lb = 999;

if unocc_lb == 999
    unocc_lb = 0;
    occ_lb = 0.1*x_data(end);
    unocc_ub = 0.5*x_data(end);
    occ_ub = x_data(end);
else
    unocc_ub = input('Enter unoccupied peak upper bound: ');
    occ_lb = input('Enter occupied peak lower bound: ');
    occ_ub = input('Enter occupied peak upper bound: ');
end

%Fit 2 gaussians to intensity histogram
[init_fit, ~] = assignment_fit(intensity_hist, x_data, unocc_lb, unocc_ub, occ_lb, occ_ub );
coeffs = coeffvalues(init_fit);

%Output of fitting, used to plot functions
gauss1 = @(x) coeffs(1)*exp(-((x - coeffs(2))/coeffs(3))^2);
gauss2 = @(x) coeffs(4)*exp(-((x - coeffs(5))/coeffs(6))^2);
gauss_sum = @(x) gauss1(x) + gauss2(x);
xinterval = [0 x_data(end)];

% Threshold intensity distinguishing occupied and unnoccupied sites to
% first order
threshold = fminbnd(gauss_sum, coeffs(2), coeffs(5)); %Minimum of sum of two gaussians

% %Plot histogram
% figure(1)
% histogram(intensity_store, edges);
% hold on
% %Plot Gaussian fits
% fplot(gauss1, xinterval);
% fplot(gauss2, xinterval);
% fplot(gauss_sum, xinterval);
% legend('', 'unoccupied', 'occupied', 'sum')
% xlabel('site intensity')
% ylabel('frequency of occurence')
% hold off


end

