function [offset] = alignment_fit(dist, spacing)
%Fit a series of up to 20 identical gaussians to an input array, output the
%horizontal offset of the fit.
%The input array corresponds to the projection of intensity data aliong the
%x or y axis of an image of a lattice.


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( [], dist );

expected_sites = round(length(dist)/spacing);
suppress = zeros( [1, 20] );
suppress(1:expected_sites) = 1;

% Set up fittype and options.
ft = fittype( @(a, b, c, x) suppress(1)*a.*exp(-((x-b)./c).^2) + suppress(2)*a.*exp(-((x-(b + spacing))./c).^2) + suppress(3)*a.*exp(-((x-(b + 2*spacing))./c).^2) + suppress(4)*a.*exp(-((x-(b + 3*spacing))./c).^2) + suppress(5)*a.*exp(-((x-(b + 4*spacing))./c).^2) + suppress(6)*a.*exp(-((x-(b + 5*spacing))./c).^2) + suppress(7)*a.*exp(-((x-(b + 6*spacing))./c).^2) + suppress(8)*a.*exp(-((x-(b + 7*spacing))./c).^2) + suppress(9)*a.*exp(-((x-(b + 8*spacing))./c).^2) + suppress(10)*a.*exp(-((x-(b + 9*spacing))./c).^2) + suppress(11)*a.*exp(-((x-(b + 10*spacing))./c).^2) + suppress(12)*a.*exp(-((x-(b + 11*spacing))./c).^2) + suppress(13)*a.*exp(-((x-(b + 12*spacing))./c).^2) + suppress(14)*a.*exp(-((x-(b + 13*spacing))./c).^2) + suppress(15)*a.*exp(-((x-(b + 14*spacing))./c).^2) + suppress(16)*a.*exp(-((x-(b + 15*spacing))./c).^2) + suppress(17)*a.*exp(-((x-(b + 16*spacing))./c).^2) + suppress(18)*a.*exp(-((x-(b + 17*spacing))./c).^2) + suppress(19)*a.*exp(-((x-(b + 18*spacing))./c).^2) + suppress(20)*a.*exp(-((x-(b + 19*spacing))./c).^2) , 'independent', 'x', 'dependent', 'y', 'coefficients', {'a', 'b', 'c'} );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [5e5, 10, 10];
opts.Lower = [0, -spacing, 0];
opts.Upper = [Inf, spacing, spacing];

% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft, opts );
coeffs = coeffvalues(fitresult);
offset = coeffs(2) - (spacing / 2);

% % Plot fit with data.
% figure(1)
% h = plot( fitresult, xData, yData );
% legend( h, 'hist_freq', 'untitled fit 1', 'Location', 'South' );
% % Label axes
% ylabel hist_freq
% grid on

