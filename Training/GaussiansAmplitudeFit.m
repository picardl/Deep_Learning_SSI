function [coefs, confs, FitZ] = GaussiansAmplitudeFit(inputpicture,pixelspersite, PSFwidth)
%Fits a sum of 9 Gaussians to an input image, which is assumed to be a 3x3
%lattice segment exactly centered on the central site.
%
%Arguments:
%   inputpicture: Array containing image segment to be reconstructed
%   pixelspersite: Lattice spacing, in pixels
%   PSFwidth: float, width of Gaussian PSF fit

%Produce 3 vectors for x, y and z input to fitting
[ydim, xdim] = size(inputpicture);
xc = linspace(1,xdim,xdim);
yc = linspace(1,ydim,ydim);

%x and y coordinates for fit
[X,Y] = ndgrid(xc,yc);
xcol = reshape(X,xdim*ydim,1);
ycol = reshape(Y,xdim*ydim,1);
zcol = reshape(inputpicture',xdim*ydim,1);

%Fit bounds and starting values
lowervalues = zeros(9,1);
uppervalues = repmat(2*max(max(inputpicture)),9,1);
startvalues = repmat(max(max(inputpicture))/2,9,1);

%Fit options
options = fitoptions('Method','NonlinearLeastSquares','Lower',lowervalues,...
                        'Upper',uppervalues,'StartPoint',startvalues);                  
ft = fittype(@(a1,a2,a3,a4,a5,a6,a7,a8,a9, x, y) myfun(a1,a2,a3,a4,a5,a6,a7,a8,a9, x, y), 'independent', {'x', 'y'}, 'dependent', 'z', 'coefficients', {'a1','a2','a3','a4','a5','a6','a7','a8','a9'});
[cfun,~,~] = fit([xcol, ycol], zcol, ft, options);

coefs = coeffvalues(cfun);
confs = confint(cfun);

%Do fit
FitZ = reshape(myfun(coefs(1), coefs(2),coefs(3),coefs(4),coefs(5),coefs(6),coefs(7),coefs(8),coefs(9), xcol, ycol),xdim,ydim);

function zt = myfun(a1,a2,a3,a4,a5,a6,a7,a8,a9, x, y)
    %Sum of 2D Gaussians with fixed width and positions
    zt = a1 .* exp(-(x - 0.5*pixelspersite) .^2 ./ (2 * PSFwidth ^ 2) - (y - 0.5*pixelspersite) .^2 / (2 * PSFwidth ^2))...
        + a2 .* exp(-(x - 1.5*pixelspersite) .^2 ./ (2 * PSFwidth ^ 2) - (y - 0.5*pixelspersite) .^2 / (2 * PSFwidth ^2))...
        + a3 .* exp(-(x - 2.5*pixelspersite) .^2 ./ (2 * PSFwidth ^ 2) - (y - 0.5*pixelspersite) .^2 / (2 * PSFwidth ^2))...
        + a4 .* exp(-(x - 0.5*pixelspersite) .^2 ./ (2 * PSFwidth ^ 2) - (y - 1.5*pixelspersite) .^2 / (2 * PSFwidth ^2))...
        + a5 .* exp(-(x - 1.5*pixelspersite) .^2 ./ (2 * PSFwidth ^ 2) - (y - 1.5*pixelspersite) .^2 / (2 * PSFwidth ^2))...
        + a6 .* exp(-(x - 2.5*pixelspersite) .^2 ./ (2 * PSFwidth ^ 2) - (y - 1.5*pixelspersite) .^2 / (2 * PSFwidth ^2))...
        + a7 .* exp(-(x - 0.5*pixelspersite) .^2 ./ (2 * PSFwidth ^ 2) - (y - 2.5*pixelspersite) .^2 / (2 * PSFwidth ^2))...
        + a8 .* exp(-(x - 1.5*pixelspersite) .^2 ./ (2 * PSFwidth ^ 2) - (y - 2.5*pixelspersite) .^2 / (2 * PSFwidth ^2))...
        + a9 .* exp(-(x - 2.5*pixelspersite) .^2 ./ (2 * PSFwidth ^ 2) - (y - 2.5*pixelspersite) .^2 / (2 * PSFwidth ^2));
end

end

