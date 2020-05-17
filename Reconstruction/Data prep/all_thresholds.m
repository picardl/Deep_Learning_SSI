function [ all_thresh, lower_bounds, upper_bounds ] = all_thresholds( images_3by3, patterns_3by3, lattice_spacing, PSFwidth)
%ALL_THRESHOLDS finds threshold intensities for different combinations of
%neighbouring site occupations
%   
% 	Inputs:       images_3by3 - 3D matrix of input 3x3 images containing all possible combinations of neighbours
%                 patterns_3by3 - 3D matrix of filling patterns of input images
%                 lattice_spacing - lattice spacing in pixels.
%                 PSFwidth – width of the point spread function (in pixels) for the given imaging parameters.
% 
% 	Outputs:      threshold - intensity threshold for naive assignment.
%                 lower_uncert_bound - lower intensity bound below which sites can be assigned as unoccupied with confidence greater than bound_confidence.
%                 upper_uncert_bound - upper intensity bound above which sites can be assigned as occupied with confidence greater than bound_confidence.

%Number of pictures in input data set
num_pics = size(images_3by3, 3);

%Initialise arrays to 0
diag_neighb = zeros(num_pics, 1);
near_neighb = zeros(num_pics, 1);
all_thresh = zeros(5, 5);
lower_bounds = zeros(5, 5);
upper_bounds = zeros(5, 5);

%Determine number of diagonal and nearest neighbours for each input picture
for i = 1:num_pics
    diag_neighb(i) = patterns_3by3(1, 1, i) + patterns_3by3(1, 3, i) + patterns_3by3(3, 1, i) + patterns_3by3(3, 3, i);
    near_neighb(i) = patterns_3by3(1, 2, i) + patterns_3by3(2, 1, i) + patterns_3by3(2, 3, i) + patterns_3by3(3, 2, i);
end

for j = 1:5
    for k = 1:5
        
        %Group all images with a given combination of nearest and diagonal
        %neighours and determine central site intensity for each image
        threshpics = images_3by3(:, :, near_neighb == j - 1 & diag_neighb == k - 1);
        intens = zeros(size(threshpics, 3), 1);
        for m = 1:size(threshpics, 3);
            intens(m) = lattice_intensity( threshpics(:, :, m), [1.5*lattice_spacing, 1.5*lattice_spacing], lattice_spacing, PSFwidth, 0, 0);
        end
        
        %Generate histogram of intensities
        edges = linspace(0, max(intens), 50);
        intensity_hist = histcounts(intens, edges);
        x_data = edges(1:end-1) + (edges(2)-edges(1))/2;
        
        %Bounds of means for occupied and unoccpied intensity distribution
        unocc_lb = 0;
        occ_lb = 0.1*x_data(end);
        unocc_ub = 0.7*x_data(end);
        occ_ub = x_data(end);
        
         %Fit 2 gaussians to intensity histogram
        [init_fit, ~] = assignment_fit(intensity_hist, x_data, unocc_lb, unocc_ub, occ_lb, occ_ub );
        coeffs = coeffvalues(init_fit);

        %Output of fitting, used to plot functions
        gauss1 = @(x) coeffs(1)*exp(-((x - coeffs(2))/coeffs(3))^2);
        gauss2 = @(x) coeffs(4)*exp(-((x - coeffs(5))/coeffs(6))^2);
        gauss_diff = @(x) abs(gauss1(x) - gauss2(x));
        gauss_sum = @(x) gauss1(x) + gauss2(x);
        xinterval = [0 x_data(end)];
        
        %If peaks are in wrong order (cannot always be avoided), switch
        %their indices
        if coeffs(2) > coeffs(5)
            temp = coeffs(4:6);
            coeffs(4:6) = coeffs(1:3);
            coeffs(1:3) = temp;
        end
            
        % Threshold intensity distinguishing occupied and unnoccupied sites
        all_thresh(j, k) = fminbnd(gauss_diff, coeffs(2), coeffs(5)); %Minimum of sum of two gaussians
        lower_bounds(j, k) = coeffs(2) - coeffs(3)/sqrt(2);
        upper_bounds(j, k) = coeffs(5) + coeffs(6)/sqrt(2);
        
%         figure((j - 1)*5 + k)
%         histogram(intens, edges);
%         hold on
%         %Plot Gaussian fits
%         fplot(gauss1, xinterval);
%         fplot(gauss2, xinterval);
%         fplot(gauss_sum, xinterval);
%         legend('', 'unoccupied', 'occupied', 'sum')
%         xlabel('site intensity')
%         ylabel('frequency of occurence')
%         hold off
    end
end
     
