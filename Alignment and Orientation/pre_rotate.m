function [ rot_centers ] = pre_rotate( angle, centers, midpoint)
%Rotates a set of cartesian coordinates clockwise by a specified angle
%   Used for lattice alignment testing

centers(:,1) = centers(:,1) - midpoint(2);
centers(:,2) = midpoint(1) - centers(:,2);

radii = sqrt(centers(:,1).^2 + centers(:,2).^2);
init_angles = atan2(centers(:,2),centers(:,1));
new_angles = init_angles - angle;
rot_centers = [radii.*cos(new_angles), radii.*sin(new_angles)];


end

