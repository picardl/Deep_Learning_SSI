clear all
u = 1.66053892e-27;
m = 168*u;
c = 299792458;
hbar = 1.05457173e-34;
recoilvel = 5.9e-3;
scatteringrate = 2*pi*30e6;
EMGain = 1000;
lambda = 401;

paramA = [5]*1e-6;

for indA = 1:length(paramA)

NA = 0.75;
CollectedPhotonsratio = (sin(asin(NA)/2)^2)*0.74*0.8; %QE of Andor and 10%loss
imagingpulse = paramA(indA);
pixelpersite = 15;
latticespacing = 532/2;
realpixelsize = latticespacing/pixelpersite; %in nm
sensorsize = pixelpersite*5;
PSFwidth = (0.43*lambda/2/NA)/realpixelsize;
PSFwidth2 = (2.44*lambda/2/NA)/realpixelsize;
wavepacketwidth = 0*1/10*pixelpersite;
AddNoise = 0.000; %ixon897
Lossrate = 1e-4;
koord = [0:sensorsize-1]*realpixelsize;
[X,Y] = ndgrid(-sensorsize/2:1:sensorsize/2,-sensorsize/2:1:sensorsize/2);

iterations = 10000;
nrphotonsdistdetected = zeros(iterations,1);
nrphotonsdistscattered = zeros(iterations,1);

PhotonCount = zeros(iterations,1);
AtomCount = zeros(iterations,1);

CCDsum = zeros(sensorsize,sensorsize);

tic
for iter = 1:iterations
    CCD = zeros(sensorsize,sensorsize);
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
        %A = (round(rand)*2-1);
        %velkick = [A 0 0]*recoilvel;
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
            nrphotons = nrphotons+1;
            [Rx,Ry] = randPSF(X,Y,PSFwidth2);
            posx = round(Rx+pos(1)*1e9/realpixelsize);
            posy = round(Ry+pos(2)*1e9/realpixelsize);
            if ((posx>0)&&(posx<sensorsize)&&(posy>0)&&(posy<sensorsize))
                CCD(posx,posy) = CCD(posx,posy)+1;
            end                
        end
    end
    endvel(iter) = sqrt(vel(1)^2+vel(2)^2+vel(3)^2);
    endpos(iter) = sqrt(pos(1)^2+pos(2)^2+pos(3)^2);
    nrphotonsdistdetected(iter) = nrphotons;
    nrphotonsdistscattered(iter) = nrphotonsscat;
    PhotonCount(iter) = sum(sum(CCD((2*pixelpersite+1):(3*pixelpersite),(2*pixelpersite+1):(3*pixelpersite))));
    CCDsum = CCDsum+CCD;
    if ((iter/100) == floor(iter/100))
        figure(4)
        image(CCDsum/max(max(CCDsum))*63);
        title(strcat('iter:',num2str(iter)));
        drawnow
    end
end
toc
[X,Y] = find(CCDsum>(max(max(CCDsum))/2));
FWHM = (max(X)-min(X)+max(Y)-min(Y))/2*realpixelsize;
PSFreal(indA) = FWHM/2.35482*2.8372;

figure(4)
surf(koord,koord,CCDsum,'FaceColor','interp','EdgeColor','none');


options = fitoptions('Method','NonlinearLeastSquares',...
                        'Lower',[0 0 0 0],...
                        'Upper',[1500 1500 1000 100000],...
                        'StartPoint',[600,600,100,500]);
ft = fittype(@(posx,posy,sigma,A, x, y) Gauss2D(posx,posy,sigma,A, x, y), 'numindep', 2);
[X,Y] = ndgrid(koord,koord);
xcol = reshape(X,length(koord)^2,1);
ycol = reshape(Y,length(koord)^2,1);
zcol = reshape(CCDsum,length(koord)^2,1);
[cfun,gof,output] = fit([xcol, ycol], zcol, ft,options);
fitcoeff = coeffvalues(cfun);
PSFrealFit(indA) = fitcoeff(3)*2.8372;

title(strcat('FWHM-Width:',num2str(PSFreal(indA)),'nm Fit-Width:',num2str(PSFrealFit(indA)),'nm'));

figure(2)
histogram(endvel)
figure(3)
histogram(endpos)


end

% figure
% plot(paramA,PSFreal,paramA,PSFrealFit);
