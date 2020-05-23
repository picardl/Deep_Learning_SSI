clear all
%load('PSFFitER20pixel532nm2us08NA');

paramA = [4:20];
paramB = [500:50:1100];

for indA = 1:length(paramA)
for indB = 1:length(paramB)

NA = 0.75;
CollectedPhotonsratio = (sin(asin(NA)/2)^2)*0.74*0.8; %QE of Andor and 20%loss
lambda = 401;
scatteringrate = 2*pi*30e6;
imagingpulse = 2e-6;%paramA(indA);
pixelpersite = paramA(indA);
latticespacing = paramB(indB)/2;
realpixelsize = latticespacing/pixelpersite; %in nm
sites = 20;
sensorsize = pixelpersite*sites;
PSFwidth = (0.43*lambda/2/NA)/latticespacing*pixelpersite;
wavepacketwidth = 1/10*pixelpersite;
EMGain = 1000;
AddNoise = 0.005; %ixon897
Lossrate = 1e-4;
photons = scatteringrate*imagingpulse/2*CollectedPhotonsratio*exp(-scatteringrate*imagingpulse/2*Lossrate);

u = 1.66053892e-27;
m = 168*u;
c = 299792458;
hbar = 1.05457173e-34;
recoilvel = 5.9e-3;

iterations = 60;
nrphotonsdistdetected = zeros(iterations,1);
nrphotonsdistscattered = zeros(iterations,1);
PhotonCount = [];
ElectronCount = [];
AtomCount = [];
tic
parfor iter = 1:iterations

CCD = zeros(sensorsize,sensorsize);
pattern = round(rand(sites,sites));
% pattern = zeros(sites,sites);
% pattern(2:end-1,2:end-1) = round(rand(sites-2,sites-2)-0.2);

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
        %velkick = ([rand rand rand]-0.5)*2;
        %velkick = velkick/norm(velkick)*recoilvel;
        vel = vel+velkick;
        %Probability of loss
        if (rand<Lossrate)
            break
        end
        %Probability of Detection
        if (rand<CollectedPhotonsratio)
            nrphotons = nrphotons+1;
            posx = round(PSFwidth*randn+(inda-1/2)*pixelpersite+pos(1)*1e9/realpixelsize);
            posy = round(PSFwidth*randn+(indb-1/2)*pixelpersite+pos(2)*1e9/realpixelsize);
            if ((posx>0)&&(posx<sensorsize)&&(posy>0)&&(posy<sensorsize))
                CCD(posx,posy) = CCD(posx,posy)+1;
            end                
        end
    end
    nrphotonsdistdetected(iter) = nrphotons;
    nrphotonsdistscattered(iter) = nrphotonsscat;
end
end
%sum(sum(CCD))

%noise
for ind = 1:(sensorsize^2)
    CCD(ind) = CCD(ind)+randpoisson(AddNoise);
    if (CCD(ind)>0)
        CCD(ind) = randEMGain(CCD(ind),EMGain);
    end
end

for ind1 = 1:sites
    for ind2 = 1:sites
        PhotonCount = [PhotonCount sum(sum(CCD(((ind1-1)*pixelpersite+1):(ind1*pixelpersite),((ind2-1)*pixelpersite+1):(ind2*pixelpersite))>0))];
        ElectronCount = [ElectronCount sum(sum(CCD(((ind1-1)*pixelpersite+1):(ind1*pixelpersite),((ind2-1)*pixelpersite+1):(ind2*pixelpersite))))];
        AtomCount = [AtomCount pattern(ind1,ind2)];
    end
end
             
%PhotonCount(iter) = sum(sum(CCD((pixelpersite+1):(2*pixelpersite),(pixelpersite+1):(2*pixelpersite))));
%AtomCount(iter) = pattern(5);

% %single picture analysis
% figure(1)
% image(CCD/max(max(CCD))*63);
% view(0,90);
% figure(2)
% image(pattern*63)
% figure(3)
% PhotonCountimg = (reshape(PhotonCount,sites,sites))';
% ElectronCountimg = (reshape(ElectronCount,sites,sites))';
% b = bar3(ElectronCountimg);
% for k = 1:length(b)
%     zdata = b(k).ZData;
%     b(k).CData = zdata;
%     b(k).FaceColor = 'interp';
% end
% %image(ElectronCountimg)
% drawnow
% patternreconstruct = zeros(sites,sites);
% patternreconstruct(ElectronCountimg<1) = -1;
% patternreconstruct(ElectronCountimg>((photons-1)*EMGain)) = +1;
% for ind1 = 1:sites
%     for ind2 = 1:sites
% 
%     end
% end

end
toc
% figure(4)
% histogram(nrphotonsdistdetected);
% hold on
% histogram(nrphotonsdistscattered);
% hold off
% 
ZeroAtoms = sort(PhotonCount(AtomCount < 1));
OneAtom = sort(PhotonCount(AtomCount > 0));
%bins = [0:EMGain:max(OneAtom)*1.05];
bins = [0:max(OneAtom)+2];
figure(3)
histogram(ZeroAtoms,bins);
hold on
histogram(OneAtom,bins);
%boundary = (EMGain*(photons+2*AddNoise*pixelpersite^2)/2);
boundary = (1*(photons+2*AddNoise*pixelpersite^2)/2);
%Boundary search
bound = 0;
while (bound < length(OneAtom))
    bound = bound+1;
    
    if (length(ZeroAtoms(ZeroAtoms > OneAtom(bound))) < bound)
        boundary = OneAtom(bound-1);
        break
    end
end
FidelityPhotons(indA,indB) = length(OneAtom(OneAtom>boundary))/length(OneAtom);
title(strcat('Treshold>',int2str(boundary),' Fidelity:',num2str(FidelityPhotons(indA,indB))));
hold off

ZeroAtoms = sort(ElectronCount(AtomCount < 1));
OneAtom = sort(ElectronCount(AtomCount > 0));
bins = [0:EMGain:max(OneAtom)*1.05];
%bins = [0:max(OneAtom)+2];
figure(4)
histogram(ZeroAtoms,bins);
hold on
histogram(OneAtom,bins);
boundary = (EMGain*(photons+2*AddNoise*pixelpersite^2-1)/2);
%Boundary search
bound = 0;
while (bound < length(OneAtom))
    bound = bound+1;
    if (length(ZeroAtoms(ZeroAtoms > OneAtom(bound))) < bound)
        boundary = OneAtom(bound-1);
        break
    end
end
%boundary = (1*(photons+2*AddNoise*pixelpersite^2)/2);
FidelityElectrons(indA,indB) = length(OneAtom(OneAtom>boundary))/length(OneAtom);
title(strcat(num2str(paramA(indA)),' ',num2str(paramB(indB)),'Treshold>',int2str(boundary),' Fidelity:',num2str(FidelityElectrons(indA,indB))));
hold off

% ZeroAtoms = PhotonCountDecon(AtomCount < 1);
% OneAtom = PhotonCountDecon(AtomCount > 0);
% bins = linspace(0,max(OneAtom),100);
% figure(4)
% histogram(ZeroAtoms,bins);
% hold on
% histogram(OneAtom,bins);
% boundary = (EMGain*(photons+2*AddNoise*pixelpersite^2)/2);
% Fidelity(indA,indB) = length(find(OneAtom>boundary))/length(OneAtom);
% title(strcat('Treshold:',int2str(boundary),' Fidelity:',num2str(Fidelity(indA,indB))));
% hold off

drawnow
end
end
% 
figure(5)
surf(paramB,paramA,FidelityPhotons,'FaceColor','interp');
view(0,90);
figure(6)
surf(paramB,paramA,FidelityElectrons,'FaceColor','interp');
view(0,90);
%figure(5)
%[decovpic,P1] = deconvblind(CCD,PSFFit);
%surf(abs(decovpic)/max(max(abs(decovpic)))*63);
%figure(4)
%image(abs(P1)/max(max(abs(decovpic)))*63);
% figure(5)
% surf(paramB,paramA,1-FidelityElectrons,'FaceColor','interp');
% view(0,90);
