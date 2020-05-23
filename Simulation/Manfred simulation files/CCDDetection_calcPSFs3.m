clear all
warning off
u = 1.66053892e-27;
m = 168*u;
c = 299792458;
hbar = 1.05457173e-34;
recoilvel = 5.9e-3;
scatteringrate = 2*pi*30e6;
EMGain = 1000;
lambda = 401;
iterations = 200000;

paramA = [5];
paramB = [532 1064 500:50:1100];
paramC = [0.5:0.05:0.9];
paramD = [1:0.5:5]*1e-6;

for indA = 1:length(paramA)
    
PSFstruct = zeros(length(paramB),length(paramC),length(paramD),paramA(indA)*5,paramA(indA)*5);
    
for indB = 1:length(paramB)
for indC = 1:length(paramC)
for indD = 1:length(paramD)

NA = paramC(indC);
CollectedPhotonsratio = (sin(asin(NA)/2)^2)*0.74*0.8; %QE of Andor and 10%loss
imagingpulse = paramD(indD);
photons = scatteringrate*imagingpulse/2*CollectedPhotonsratio;
pixelpersite = paramA(indA);
latticespacing = paramB(indB)/2;
realpixelsize = latticespacing/pixelpersite; %in nm
sensorsize = pixelpersite*5;
PSFwidth = (0.43*lambda/2/NA)/latticespacing*pixelpersite;
PSFwidth2 = (2.44*lambda/2/NA)/latticespacing*pixelpersite;
wavepacketwidth = 1/10*pixelpersite;
Lossrate = 1e-4;
koord = [0:sensorsize-1]*latticespacing/pixelpersite;
PhotonCount = zeros(iterations,1);
AtomCount = zeros(iterations,1);
CCDsum = zeros(sensorsize,sensorsize);
[X,Y] = ndgrid(-sensorsize/2:1:sensorsize/2,-sensorsize/2:1:sensorsize/2);

name = strcat('PSF\PSF',int2str(pixelpersite),'pixel',int2str(latticespacing*2),'nm',int2str(NA*100),'NA',int2str(imagingpulse*1e7),'us401nm_real')
try
    load(name);
catch
iter2 = 0;
tic
while (max(max(CCDsum))<10000)
    iter2 = iter2+1000;
for iter = 1:1000

    CCD = zeros(sensorsize,sensorsize);
    pattern = zeros(5,5);
    %AtompatternX = wavepacketwidth*randn+sensorsize/2;
    %AtompatternY = wavepacketwidth*randn+sensorsize/2;
    AtompatternX = wavepacketwidth*randn;%
    AtompatternY = wavepacketwidth*randn;%
    pattern(13) = 1;
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
            nrphotons = nrphotons+1;
            [Rx,Ry] = randPSF(X,Y,PSFwidth2);
            posx = round(Rx+AtompatternX+pos(1)*1e9/realpixelsize);
            posy = round(Ry+AtompatternY+pos(2)*1e9/realpixelsize);
            %posx = round(PSFwidth*randn+AtompatternX+pos(1)*1e9/realpixelsize);
            %posy = round(PSFwidth*randn+AtompatternY+pos(2)*1e9/realpixelsize);
            if ((posx>0)&&(posx<sensorsize)&&(posy>0)&&(posy<sensorsize))
                CCD(posx,posy) = CCD(posx,posy)+1;
            end                
        end
    end
    CCDsum = CCDsum+CCD;
end

figure(1)
surf(CCDsum);
%image(PSFFit/max(max(PSFFit))*63);
title(strcat(name,int2str(iter2)));
drawnow
end
toc

PSFFit = CCDsum/sum(sum(CCDsum));
save(name,'PSFFit');


end

PSFstruct(indB,indC,indD,:,:) = PSFFit;
name2 = strcat('PSF',int2str(pixelpersite),'pixel401nm_real');
save(name2,'PSFstruct');

end
end
end
end

