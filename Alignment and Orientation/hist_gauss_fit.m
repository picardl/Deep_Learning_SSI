function [fitresult, gof] = hist_gauss_fit(hist, expected_spacing)
%CREATEFIT(HIST_FREQ)
%  Create a fit.
%
%  Data for 'hist gaussian fit' fit:
%      Y Output: hist_freq
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 03-Aug-2016 17:19:21


%% Fit: 'hist gaussian fit'.
[xData, yData] = prepareCurveData( [], hist );

a_start = max(hist);

% Set up fittype and options.
ft = fittype( 'a*exp(-((x-b)/c)^2) + a*exp(-((x-(b + d))/c)^2) + a*exp(-((x-(b + 2*d))/c)^2) + a*exp(-((x-(b + 3*d))/c)^2) + a*exp(-((x-(b + 4*d))/c)^2) + a*exp(-((x-(b + 5*d))/c)^2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [a_start, 10, 10, expected_spacing];
opts.Lower = [0, -10, 0, expected_spacing - 5];
opts.Upper = [Inf, Inf, Inf, expected_spacing + 5];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure(1)
h = plot( fitresult, xData, yData );
% legend( h, 'hist_freq', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
ylabel('Histogram bin height')
xlabel('Mutual distance / px')
grid on


