function R = randpoisson(varlambda)
    if (varlambda<100)
        nr=1; produ = 1;
        produ = produ*rand;
        while produ >= exp(-varlambda)
            produ = produ*rand;
            nr = nr+1;
        end
        R = nr-1;
    else
        R = round(sqrt(varlambda)*randn+varlambda);
    end
end