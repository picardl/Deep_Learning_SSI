sites = 10; %simulated pattern has dimensions sites x sites
pixels_per_site = 20;
fraction_filled = 0.1;
n_pics = 100;

maximum_accuracy = zeros([10 1]);
binary_accuracy = zeros([10 1]);
avg_accuracy = zeros([10 1]);

for i = 1:10
    lattice_angle = 2*pi*(rand - 0.5);
    
    sim_pics = simulate_npictures(sites, pixels_per_site, fraction_filled, lattice_angle, n_pics)  ; %Cell array containing simulated pictures

    [maximum_centers_store, ~] = calculate_atom_centers(sim_pics, sites, pixels_per_site, fraction_filled);
    [binary_centers_store, midpoint] = OLDcalculate_atom_centers(sim_pics, sites, pixels_per_site, fraction_filled);
    
    
    f = @(angle)rotate_centers(angle, maximum_centers_store, midpoint); %anonymous function for optimisation
    maximum_output_angle = -patternsearch(f, 0.1, [], [], [], [], -lattice_angle - 0.2, -lattice_angle + 0.2); %Optimise in 0.4 radian region around expected angle
    maximum_accuracy(i) = abs(lattice_angle - maximum_output_angle);
    
    f = @(angle)rotate_centers(angle, binary_centers_store, midpoint); %anonymous function for optimisation
    binary_output_angle = -patternsearch(f, 0.1, [], [], [], [], -lattice_angle - 0.2, -lattice_angle + 0.2); %Optimise in 0.4 radian region around expected angle
    binary_accuracy(i) = abs(lattice_angle - binary_output_angle);
    
    avg_angle = (maximum_output_angle + binary_output_angle)/2;
    avg_accuracy(i) = abs(lattice_angle - avg_angle);
end

plot(maximum_accuracy);
hold on
plot(binary_accuracy);
