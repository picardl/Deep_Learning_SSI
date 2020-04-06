clear all
u = 1.66053892e-27;
m = 168*u;
c = 299792458;
hbar = 1.05457173e-34;
recoilvel = 5.9e-3;
scatteringrate = 2*pi*30e6;
EMGain = 1000;
lambda = 401;
iterations = 40;

paramA = [11];
paramB = [532];
paramC = [0.75];
paramD = [5]*1e-6;
  
NA = paramC;
CollectedPhotonsratio = (sin(asin(NA)/2)^2)*0.74*0.8; %QE of Andor and 10%loss
imagingpulse = paramD;
photons = scatteringrate*imagingpulse/2*CollectedPhotonsratio;
pixelpersite = paramA;
latticespacing = paramB/2;
realpixelsize = latticespacing/pixelpersite; %in nm
sensorsize = pixelpersite*7;
PSFwidth = (0.43*lambda/2/NA)/latticespacing*pixelpersite;
PSFwidth2 = (2.44*lambda/2/NA)/latticespacing*pixelpersite;
wavepacketwidth = 1/10*pixelpersite;
Lossrate = 1e-4;
koord = [0:sensorsize-1]*latticespacing/pixelpersite;
PhotonCount = zeros(iterations,1);
AtomCount = zeros(iterations,1);
CCDsum = zeros(sensorsize,sensorsize);
nrphotons = 0;
[X,Y] = ndgrid(-sensorsize/2:1:sensorsize/2,-sensorsize/2:1:sensorsize/2);

while(true)

tic
parfor iter = 1:iterations

    CCD = zeros(sensorsize,sensorsize);
    pattern = zeros(5,5);
    AtompatternX = wavepacketwidth*randn;%+sensorsize/2;
    AtompatternY = wavepacketwidth*randn;%+sensorsize/2;
    pattern(13) = 1;
    %Simulate photon scattering
    t = 0;

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
            nrphotons = nrphotons+1;
            [Rx,Ry] = randPSF(X,Y,PSFwidth2);
            posx = round(Rx+AtompatternX+pos(1)*1e9/realpixelsize);
            posy = round(Ry+AtompatternY+pos(2)*1e9/realpixelsize);
%             posx = round(PSFwidth*randn+AtompatternX+pos(1)*1e9/realpixelsize);
%             posy = round(PSFwidth*randn+AtompatternY+pos(2)*1e9/realpixelsize);
            if ((posx>0)&&(posx<sensorsize)&&(posy>0)&&(posy<sensorsize))
                CCD(posx,posy) = CCD(posx,posy)+1;
            end                
        end
    end
    CCDsum = CCDsum+CCD;
end
toc

figure(1)
surf(CCDsum);
title(int2str(nrphotons));
drawnow

end


