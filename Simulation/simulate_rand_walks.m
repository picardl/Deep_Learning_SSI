%Simulate random walks and photon detection positions of atoms

%Simulation parameters
n = 5; %number of particle walks to simulate
pixelpersite = 20;

%Define simulation constants
NA = 0.85;
latticespacing = 256; %in nm
imagingpulse = 3e-6;
Lossrate = 1e-4;
recoilvel = 5.9e-3;
CollectedPhotonsratio = (sin(asin(NA)/2)^2)*0.74*0.8; %QE of Andor and 20%loss
lambda = 401;
scatteringrate = 2 * pi *30e6;
realpixelsize = latticespacing/pixelpersite;  %in nm
sensorsize = 2*pixelpersite; %#ok<*NASGU>
PSFwidth = (0.43 * lambda / 2 / NA)/latticespacing*pixelpersite;
wavepacketwidth = 1 / 10 * pixelpersite;
EMGain = 1000;
AddNoise = 0.005; %ixon897
photons = scatteringrate * imagingpulse / 2 * CollectedPhotonsratio * exp(-scatteringrate * imagingpulse / 2 * Lossrate);
u = 1.66053892e-27;
m = 168 * u;
c = 299792458;
hbar = 1.05457173e-34;

all_pictures = zeros(sensorsize, sensorsize, n); %CCD picture array
all_walks = cell(n, 1); %Cell array to hold random walks
all_detections = cell(n, 1);
all_detect_steps = cell(n, 1);
all_vel_steps = cell(n, 1);

for i = 1:n
    PhotonCount = [];
    ElectronCount = [];
    AtomCount = [];
    walk = [];
    detect_pos = [];
    detect_step = [];
    vel_step = [];

    CCD = zeros(sensorsize,sensorsize);

    inda = 1.5; indb = 1.5; %Define single occupied site for simulation

    %Simulate photon scattering
    t = 0;
    nrphotons = 0;
    nrphotonsscat = 0;
    pos = [0 0 0];
    vel = [0 0 0];
    
    while (t<imagingpulse)
        %Absorption
        dt = randexp(scatteringrate);
        pos = pos+vel*dt;
        walk = [walk; pos];
        vel_step = [vel_step; vel];
        t = t+dt;
        A = (floor(rand*3.9999));
        switch A
            case 0
                velkick = [0 1 0]*recoilvel;
            case 1
                velkick = [0 -1 0]*recoilvel;
            case 2
                velkick = [1 0 0]*recoilvel;
            case 3
                velkick = [-1 0 0]*recoilvel;
        end
        vel = vel+velkick;
        nrphotonsscat = nrphotonsscat+1;
        %Emmission
        dt = randexp(scatteringrate);
        pos = pos+vel*dt;
        t = t+dt;
        phi = 2*pi*rand;
        theta = asin((rand-0.5)*2)+pi/2;
        velkick = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]*recoilvel;
        vel = vel+velkick;
        %Probability of loss
        if (rand<Lossrate)
            break
        end
        %Probability of Detection
        if (rand<CollectedPhotonsratio)
            r1 = randn;
            r2 = randn;
            posx = round(PSFwidth*r1+(inda-1/2)*pixelpersite+pos(1)*1e9/realpixelsize);
            posy = round(PSFwidth*r2+(indb-1/2)*pixelpersite+pos(2)*1e9/realpixelsize);
            if ((posx>0)&&(posx<=sensorsize)&&(posy>0)&&(posy<=sensorsize))
                CCD(posx,posy) = CCD(posx,posy)+1;
                detect_pos = [detect_pos; [(PSFwidth*r1*realpixelsize+pos(1)*1e9), (PSFwidth*r2*realpixelsize+pos(2)*1e9)]];
                detect_step = [detect_step; size(walk, 1)];
            end                
        end
    end
    

%noise
for ind = 1:(sensorsize^2)
    CCD(ind) = CCD(ind)+randpoisson(AddNoise);
    if (CCD(ind)>0)
        CCD(ind) = randEMGain(CCD(ind),EMGain);
    end
end

all_pictures( :, :, i) = CCD;
all_walks{i} = walk;
all_detections{i} = detect_pos;
all_detect_steps{i} = detect_step;
all_vel_steps{i} = vel_step;

end

col = [rand rand rand];
walk_i = all_walks{1}.*(10^9);
detect_i = all_detections{1};
detect_step_i = all_detect_steps{1};

%Dynamically plot random walk of one atom
figure(1)
subplot(2, 1, 1)
xlabel('x displacement from origin / nm')
ylabel('y displacement from origin / nm')
hold on
plot(walk_i(1, 1), walk_i(1, 2), 'Color', col)
xlim([-300 300])
ylim([-300 300])
for j = 2:size(walk_i, 1)
    plot(walk_i(1:j, 1), walk_i(1:j, 2), 'Color', col)
    if any(detect_step_i == j)
        detect_ind = find(detect_step_i == j);
        scatter(detect_i(detect_ind, 1), detect_i(detect_ind, 2), 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none')
    end
    pause(0.001)
end
scatter(detect_i(:, 1), detect_i(:, 2), 'MarkerEdgeColor', col)
hold off

%Dynamically plot atom velocity over time
subplot(2, 1, 2)
hold on
ylabel('velocity / m s^{-1}')
xlabel('time / s')
vel_i = all_vel_steps{1};
vel_i_tot = sqrt(sum(vel_i.^2, 2));
plot(linspace(0, 3e-6, length(vel_i)), vel_i_tot, 'Color', col);

hold off

%Plot all n random walks
for i = 1:n
    col = [rand rand rand];
    detect_i = all_detections{i};
    walk_i = all_walks{i}.*(10^9);
    figure(2)
    hold on
    for j = 1:size(walk_i, 1)
        plot(walk_i(:, 1), walk_i(:, 2), 'Color', col)
    end
    scatter(detect_i(:, 1), detect_i(:, 2), 'MarkerEdgeColor', col)
    figure(3)
    hold on
    vel_i = all_vel_steps{i};
    vel_i_tot = sqrt(sum(vel_i.^2, 2));
    plot(linspace(0, 3e-6, length(vel_i)), vel_i_tot, 'Color', col);
    ylabel('velocity / m s^{-1}')
    xlabel('time / s')
end
hold off
