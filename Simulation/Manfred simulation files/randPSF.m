function [Rx,Ry] = randPSF(x,y,PSFwidth)
    L = length(x);
    lengthpix = L*L;
    Pno = PSF2D(x,y,PSFwidth);
    Pn = reshape(Pno,lengthpix,1);
    Rx = L;
    Ry = L;
    prob = rand;%*sum(Pn);
    cs = 0;
    for ind = 1:lengthpix
        cs = cs+Pn(ind);
        if (cs>prob)
            [Rx,Ry] = ind2sub(L,ind);
            break
        end
    end
end
