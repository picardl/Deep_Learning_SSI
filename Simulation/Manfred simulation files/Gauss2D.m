function R = Gauss2D(posx,posy,sigma,A,x,y)
    
    R = A*exp(-((x-posx).^2)/(2*sigma.^2)-((y-posy).^2)/(2*sigma.^2));

end