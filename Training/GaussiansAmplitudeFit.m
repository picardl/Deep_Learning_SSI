function [coefs, confs, FitZ] = GaussiansAmplitudeFit(inputpicture,pixelspersite, PSFwidth)
%Fits a sum of 9 Gaussians to an input image

%Produce 3 vectors for x, y and z input to fitting
[ydim, xdim] = size(inputpicture);
xc = linspace(1,xdim,xdim);
yc = linspace(1,ydim,ydim);

[X,Y] = ndgrid(xc,yc);

xcol = reshape(X,xdim*ydim,1);
ycol = reshape(Y,xdim*ydim,1);
zcol = reshape(inputpicture',xdim*ydim,1);

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

%Plot input and fit
FitZ = reshape(myfun(coefs(1), coefs(2),coefs(3),coefs(4),coefs(5),coefs(6),coefs(7),coefs(8),coefs(9), xcol, ycol),xdim,ydim);
% figure(3)
% surf(xc,yc,FitZ);
% figure(4)
% image(inputpicture*63)


function zt = myfun(a1,a2,a3,a4,a5,a6,a7,a8,a9, x, y)
    %2D Gaussian with fixed width and position
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

