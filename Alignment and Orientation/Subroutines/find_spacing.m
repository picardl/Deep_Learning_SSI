function [ spacing ] = find_spacing( h1, max_lambda )
%Fit 5 identical, equidistant Gaussians to histogram to give more precise spacing

hist = detrend(h1(0.5 * max_lambda : 6.5 * max_lambda));
hist = hist - min(hist);

[fitresult, ~] = hist_gauss_fit(hist, max_lambda);
coeffs = coeffvalues(fitresult);
spacing = coeffs(4);

end

