function [fitresult, gof] = assignment_fit(intensity_hist, x_data, unocc_lb, unocc_ub, occ_lb, occ_ub)
%Fit a pair of gaussians to a histogram of site intensity.

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x_data, intensity_hist );

% Set up fittype and options.
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [max(yData)/10, unocc_lb, x_data(end)/20, max(yData)/20, occ_lb, x_data(end)/20]; %Both peaks must be at least 10% and 5% of max height of histogram, lower bound for occupied peak center given as function argument
opts.StartPoint = [max(yData), 0.1*x_data(end), x_data(end)/50, max(yData), 0.6*x_data(end), x_data(end)/50];
opts.Upper = [1.5*max(yData), unocc_ub, x_data(end), 1.5*max(yData), occ_ub, x_data(end)]; %Neither peak width can be more than half the full span of the x_data

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData, excludedPoints );
% legend( h, 'intensity_hist', 'Excluded intensity_hist', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% ylabel intensity_hist
% grid on


