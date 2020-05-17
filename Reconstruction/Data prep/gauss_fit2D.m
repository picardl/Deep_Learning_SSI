function [coefs] = gauss_fit2D(inputpicture, x_cent, y_cent, spacing)
%Fits a 2D gaussian to a portion of an input picture, with fixed centre
%coordinates. Returns amplitude of fit.
    %The area of the input which is used for fitting is a square around the
    %centre, with sides of length 2*spacing

%Produce 3 vectors for x, y and z input to fitting
[ydim, xdim] = size(inputpicture);
xc = linspace(1,xdim,xdim);
yc = linspace(1,ydim,ydim);

[X,Y] = ndgrid(xc,yc);

xcol = reshape(X,xdim*ydim,1);
ycol = reshape(Y,xdim*ydim,1);
zcol = reshape(inputpicture',xdim*ydim,1);

lowervalues = 0;
uppervalues = [Inf 50];
startvalues = [0.01 8];

%Fit options
options = fitoptions('Method','NonlinearLeastSquares','Lower',lowervalues,...
                        'Upper',uppervalues,'StartPoint',startvalues, 'Exclude', abs(xcol - x_cent) > spacing | abs(ycol - y_cent) > spacing);                  
ft = fittype(@(a, PSFwidth, x, y) myfun(a, PSFwidth, x, y), 'independent', {'x', 'y'}, 'dependent', 'z', 'coefficients', {'a', 'PSFwidth'});
     
[cfun,~,~] = fit([xcol, ycol], zcol, ft, options);

coefs = coeffvalues(cfun);

% %Plot input and fit
% FitZ = reshape(myfun(coefs(1), coefs(2), xcol, ycol),xdim,ydim);
% figure(3)
% surf(xc,yc,FitZ);
% hold on
% %surf(xc,yc,inputpicture);
% hold off

            
function zt = myfun(a, PSFwidth, x, y)
    %2D Gaussian with fixed width and position
    zt = a * exp(-(x - x_cent) .^2 / (2 * PSFwidth ^ 2) - (y - y_cent) .^2 / (2 * PSFwidth ^2));
end

end