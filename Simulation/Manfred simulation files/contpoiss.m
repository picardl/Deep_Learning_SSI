function [ p ] = contpoiss( x, mean )
%Continuous analogue of the Poisson distribution
if (x < 0); p = 0; return; end

p = exp(-mean).*(mean.^x)./gamma(x + 1);

end

