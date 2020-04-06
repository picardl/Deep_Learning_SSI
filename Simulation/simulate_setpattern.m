 function [ all_pictures ] = simulate_setpattern( pattern, pixelspersite, n, NA, latticespacing, imagingpulse, lossrate, recoilvel, scatteringrate, atomicmass, lambda, addnoise, collectedphotonsratio, latticedepth)
%Simulates n lattice images for a given pattern of occupied / unoccupied
%sites, returns images in a 3-dimensional array of height n
% 
% 	Inputs: pattern – binary matrix representing pattern of occupation, where ‘1’ is occupied and ‘0’ is unoccupied
% 			pixelspersite – width of site in pixels.
% 			n – number of pictures to simulate.
%           NA - numerical aperture (between 0 and 1)
%           latticespacing - distance between neighbouring lattice sites,
%           in nm
%           imagingpulse - length of imaging laser pulse, in s
%           scatteringrate - rate at which atoms scatter photons, in s^-1
%           lossrate - probability of atom loss in every scattering event
%           recoilvel - atom recoil velocity from scattering in m/s
%           collectedphotonsratio - fraction of scattered photons which are
%               collected on CCD
%           latticedepth - Depth of pinning lattice, in Hz
% 
% 	Outputs:all_pictures – (pixelspersite*sites) squared by n matrix, containing all simulated pictures.

%Define remaining simulation parameters and physical constants
realpixelsize = latticespacing/pixelspersite;  %in nm
sensorsize = round(pixelspersite * size(pattern,2)); %#ok<*NASGU>
PSFwidth = (0.43 * lambda / 2 / NA)/latticespacing*pixelspersite;
wavepacketwidth = 1 / 10 * pixelspersite;
EMGain = 1000;
photons = scatteringrate * imagingpulse / 2 * collectedphotonsratio * exp(-scatteringrate * imagingpulse / 2 * lossrate);
u = 1.66053892e-27;
m = atomicmass * u;
c = 299792458;
hbar = 1.05457173e-34;
latticespacing = latticespacing*1e-9;
latticedepth_J = latticedepth*hbar*2*pi; %lattice depth in J

rng('shuffle'); %Seed random number generator

all_pictures = zeros(sensorsize, sensorsize, n);

for i = 1:n
    PhotonCount = [];
    ElectronCount = [];
    AtomCount = [];

    CCD = zeros(sensorsize,sensorsize);

for inda = 1:size(pattern,1)
for indb = 1:size(pattern,1)
    if (pattern(inda,indb)<1)
        continue
    end
    %Simulate photon scattering
    t = 0;
    nrphotons = 0;
    nrphotonsscat = 0;
    pos = [0 0 0];
    vel = [0 0 0];
    acc = [0 0 0];
    while (t<imagingpulse)
        dt = randexp(scatteringrate);
        acc = -(2*pi*latticedepth_J/latticespacing)*sin(pi.*pos./latticespacing)...
            .*cos(pi.*pos./latticespacing)./m; %Initial acceleration
        pos = pos+vel*dt + 0.5*acc*dt^2; %update position
        t = t+dt;
        
        acc_new = -(2*pi*latticedepth_J/latticespacing)*sin(pi.*pos./latticespacing)...
            .*cos(pi.*pos./latticespacing)./m; %next step acceleration
        
        %Absorption
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
        vel = vel+velkick + 0.5*(acc + acc_new)*dt; %update velocity
        
        nrphotonsscat = nrphotonsscat+1;
        acc = acc_new; %update acceleration
        
        %Emmission
        dt = randexp(scatteringrate);
        pos = pos+vel*dt + 0.5*acc*dt^2; %update position
        acc_new = -(2*pi*latticedepth_J/latticespacing)*sin(pi.*pos./latticespacing)...
            .*cos(pi.*pos./latticespacing)./m; %next step acceleration
        t = t+dt; %update time
        phi = 2*pi*rand; %random emission direction
        theta = asin((rand-0.5)*2)+pi/2;
        velkick = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(phi)]*recoilvel;
        vel = vel+velkick + 0.5*(acc + acc_new)*dt; %update velocity
        
        acc = acc_new; %update acceleration
        
        %Probability of loss
        if (rand<lossrate)
            break
        end
        %Probability of Detection
        if (rand<collectedphotonsratio)
            posx = round(PSFwidth*randn+(inda-1/2)*pixelspersite+pos(1)*1e9/realpixelsize);
            posy = round(PSFwidth*randn+(indb-1/2)*pixelspersite+pos(2)*1e9/realpixelsize);
            if ((posx>0)&&(posx<=sensorsize)&&(posy>0)&&(posy<=sensorsize))
                CCD(posx,posy) = CCD(posx,posy)+1;
            end                
        end
    end
end
end

%noise
for ind = 1:(sensorsize^2)
    CCD(ind) = CCD(ind)+randpoisson(addnoise);
    if (CCD(ind)>0)
        CCD(ind) = randEMGain(CCD(ind),EMGain);
    end
end

all_pictures( :, :, i) = CCD;

end
end

