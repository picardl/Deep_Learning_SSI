function R = randexp(varlambda)
    R = -1*length(varlambda);
    while (R<0)
        R=-1./varlambda.*log(rand(length(varlambda), 1));
    end
end