function R = PSF2D(x,y,PSFwidth)
    x = 2*pi*sqrt(x.^2+y.^2)*0.6098/PSFwidth*2;
    R = (2*besselj(1,x)./(x)).^2;
    R(isnan(R)) = 1;
    R = R/sum(sum(R));
end
