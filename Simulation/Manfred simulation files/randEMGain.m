function R = randEMGain(input,gain)
    
    Maxoutput = 100000;
    Pn = zeros(Maxoutput,1);
    n = [input:Maxoutput];
    Pn(input:end) = ((n-input+1).^(input-1))./(factorial(input-1).*((gain-1+1/input).^input)).*exp(-(n-input+1)/(gain-1+1/input));
    
    prob = rand*sum(Pn);
    Pn_sum = cumsum(Pn);
    R = find(Pn_sum > prob, 1);
    
    if (Pn_sum(end) < prob)
        R = Maxoutput;
    end
    
end