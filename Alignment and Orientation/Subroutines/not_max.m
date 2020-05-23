function [ max_false ] = not_max( image_input, index )
%Returns true if the point in an image specified by the vector 'index' is not a
%local maximum.

test_point = image_input(index(2), index(1));
max_false = 0;
neighbours = zeros(8, 1);
im_size = size(image_input);

if (index(2)>1)&&(index(1)>1)
   neighbours(1) = image_input(index(2) - 1, index(1) - 1);
end
if index(1)>1
   neighbours(2) = image_input(index(2), index(1) - 1);
end
if (index(2)<im_size(2))&&(index(1)>1)
   neighbours(3) = image_input(index(2) + 1, index(1) - 1);
end
if index(2) > 1
   neighbours(4) = image_input(index(2) - 1, index(1));
end
if index(2)<im_size(2)
   neighbours(5) = image_input(index(2) + 1, index(1));
end
if (index(2)>1)&&(index(1) < im_size(1))
   neighbours(6) = image_input(index(2) - 1, index(1) + 1);
end
if index(1) < im_size(1)
   neighbours(7) = image_input(index(2), index(1) + 1);
end
if (index(2)<im_size(2))&&(index(1) < im_size(1))
   neighbours(8) = image_input(index(2) + 1, index(1) + 1);
end

if max(neighbours) > test_point || index(1) == 1 || index(1) == im_size(1) || index(2) == 1 || index(2) == im_size(2)
    max_false = 1;
end

end

