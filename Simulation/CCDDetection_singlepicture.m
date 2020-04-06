%Simulaten a single CCD picture of a lattice, with occupation defined by
%the matrix 'pattern' and all simulation parameters set in the calling
%function or script

PhotonCount = [];
ElectronCount = [];
AtomCount = [];

CCD = zeros(sensorsize,sensorsize);
pattern = round( rand(sites,sites) + (fractionfilled - 0.5));

pos_diff = [0 0];

for inda = 1:sites
for indb = 1:sites
    if (pattern(inda,indb)<1)
        continue
    end
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
        velkick = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(phi)]*recoilvel;
        vel = vel+velkick;
        %Probability of loss
        if (rand<Lossrate)
            break
        end
        %Probability of Detection
        if (rand<CollectedPhotonsratio)
            posx = round(PSFwidth*randn+(inda-1/2)*pixelpersite+pos(1)*1e9/realpixelsize);
            posy = round(PSFwidth*randn+(indb-1/2)*pixelpersite+pos(2)*1e9/realpixelsize);
            if ((posx>0)&&(posx<sensorsize)&&(posy>0)&&(posy<sensorsize))
                CCD(posx,posy) = CCD(posx,posy)+1;
            end                
        end
    end
    pos_diff = [pos_diff; [(inda-1/2)*pixelpersite, (indb-1/2)*pixelpersite] - pos(1:2)];
end
end

%noise
for ind = 1:(sensorsize^2)
    CCD(ind) = CCD(ind)+randpoisson(AddNoise);
    if (CCD(ind)>0)
        CCD(ind) = randEMGain(CCD(ind),EMGain);
    end
end


% %single picture analysis
% figure(2)
% image(CCD/max(max(CCD))*63);
%view(0,90);
%figure(2)
%image(pattern*63)




