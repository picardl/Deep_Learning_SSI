function [ left_bound, right_bound, top_bound, bottom_bound ] = lattice_bounds( x_proj, y_proj, SNR)
%lattice_bounds returns bounds of image region occupied by lattice
%   Given projections of a CCD image along the x and y axes and an estimate
%   of the signal-to-noise ratio, this function provides the boundaries of the
%   region of the image occupied by the lattice.


for i = 1:length(x_proj)/2
    mean_test = mean(x_proj(i : i + 5));
    if mean_test > max(x_proj)*(1/SNR)
        break
    end
end

left_bound = i;

for i = 1:length(x_proj)/2
    mean_test = mean(x_proj(end - 5 - i : end - i));
    if mean_test > max(x_proj)*(1/SNR)
        break
    end
end

right_bound = length(x_proj) + 1 - i;

for i = 1:length(y_proj)/2
    mean_test = mean(y_proj(i : i + 5));
    if mean_test > max(y_proj)*(1/SNR)
        break
    end
end

top_bound = i;

for i = 1:length(y_proj)/2
    mean_test = mean(y_proj(end - 5 - i : end - i));
    if mean_test > max(y_proj)*(1/SNR)
        break
    end
end

bottom_bound = length(y_proj) + 1 - i;

end

