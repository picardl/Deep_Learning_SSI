clear all
%load('PSFFitER20pixel532nm2us08NA');

paramA = [20];
paramB = [50:50:500];

for indA = 1:length(paramA)
for indB = 1:length(paramB)

NA = 0.75;
CollectedPhotonsratio = (sin(asin(NA)/2)^2)*0.74*0.8; %QE of Andor and 20%loss
lambda = 583;
scatteringrate = 2*pi*0.2e6;
imagingpulse = paramB(indB);
pixelpersite = paramA(indA);
latticespacing = 532/2;
realpixelsize = latticespacing/pixelpersite; %in nm
sites = 10;
sensorsize = pixelpersite*sites;
PSFwidth = (0.43*lambda/2/NA)/latticespacing*pixelpersite;
wavepacketwidth = 1/10*pixelpersite;
EMGain = 1000;
AddNoise = 0.005; %ixon897
Lossrate = 0*1e-4;
photons = paramB(indB);%scatteringrate*imagingpulse/4*CollectedPhotonsratio*exp(-scatteringrate*imagingpulse/2*Lossrate);
[x,y] = ndgrid([-pixelpersite*1.45:pixelpersite*1.45],[-pixelpersite*1.45:pixelpersite*1.45]);
PSF = PSF2D(x,y,1,PSFwidth/0.43*1.22);%Gauss2D(0,0,PSFwidth,1,x,y);

u = 1.66053892e-27;
m = 168*u;
c = 299792458;
hbar = 1.05457173e-34;
recoilvel = 0*5.9e-3;

iterations = 120;
PhotonCount = [];
ElectronCount = [];
ElectronCountDeconv = [];
AtomCount = [];
tic
parfor iter = 1:iterations
    CCD = zeros(sensorsize,sensorsize);
    CCDreconstruction = zeros(sensorsize,sensorsize);
    %pattern = round(rand(sites,sites));
    pattern = zeros(sites,sites);
    pattern(2:end-1,2:end-1) = round(rand(sites-2,sites-2)-0);
    for inda = 2:sites-1
        for indb = 2:sites-1
            if (pattern(inda,indb)>0)
            %Simulate photon scattering
                NrPhotons = poissrnd(photons);
                posx = max(min(round(PSFwidth*randn(NrPhotons,1)+(inda-1/2)*pixelpersite),sensorsize),1);
                posy = max(min(round(PSFwidth*randn(NrPhotons,1)+(indb-1/2)*pixelpersite),sensorsize),1);
                for id1 = 1:NrPhotons
                    %if ((posx(id1)>0)&&(posx(id1)<sensorsize)&&(posy(id1)>0)&&(posy(id1)<sensorsize))
                        CCD(posx(id1),posy(id1)) = CCD(posx(id1),posy(id1))+1;
                    %end
                end
            end
        end
    end
    CCDPhotons = CCD;
    %noise
    for ind = 1:(sensorsize^2)
        CCD(ind) = CCD(ind)+poissrnd(AddNoise);
        if (CCD(ind)>0)
            CCD(ind) = randEMGain(CCD(ind),EMGain);
        end
    end
    CCDdeconv = deconvlucy(CCD, PSF);
    for ind1 = 2:sites-1
        for ind2 = 2:sites-1
            PhotonCount = [PhotonCount sum(sum(CCDPhotons(((ind1-1)*pixelpersite+1):(ind1*pixelpersite),((ind2-1)*pixelpersite+1):(ind2*pixelpersite))))];
            ElectronCount = [ElectronCount sum(sum(CCD(((ind1-1)*pixelpersite+1):(ind1*pixelpersite),((ind2-1)*pixelpersite+1):(ind2*pixelpersite))))];
            ElectronCountDeconv = [ElectronCountDeconv sum(sum(CCDdeconv(((ind1-1)*pixelpersite+1):(ind1*pixelpersite),((ind2-1)*pixelpersite+1):(ind2*pixelpersite))))];
            AtomCount = [AtomCount pattern(ind1,ind2)];
        end
    end

% figure(1)
% image(CCD/max(max(CCD))*63);
% view(0,90);
% 
% figure(2)
% ElectronCountimg = (reshape(ElectronCount,sites-2,sites-2))';
% b = bar3(ElectronCountimg);
% view(0,90);
% for i = 1:length(b)
%     zdata = get(b(i),'ZData');
%     set(b(i),'CData',zdata)
%     set(b,'EdgeColor','k') 
% end
% 
% figure(3)
% image(pattern*63)
% 
% figure(4)
% image(CCDdeconv/max(max(CCDdeconv))*63);
% view(0,90);
% 
% figure(5)
% ElectronCountDeconvimg = (reshape(ElectronCountDeconv,sites-2,sites-2))';
% b = bar3(ElectronCountDeconvimg);
% view(0,90);
% for i = 1:length(b)
%     zdata = get(b(i),'ZData');
%     set(b(i),'CData',zdata)
%     set(b,'EdgeColor','k') 
% end
    
    
end
toc

ZeroAtoms = sort(PhotonCount(AtomCount < 1));
OneAtom = sort(PhotonCount(AtomCount > 0));
bins = [0:max(OneAtom)+2];
figure(1)
% [n1,xout1] = hist(ZeroAtoms,(bins(1:end-1)+((bins(2)-bins(1))/2)));
% [n2,xout2] = hist(OneAtom,(bins(1:end-1)+((bins(2)-bins(1))/2)));
% bar(xout1,n1,'b');
% hold on
% bar(xout2,n2,'r');
% hold off
histogram(ZeroAtoms,bins);
hold on
histogram(OneAtom,bins);

%Boundary search
bound = 1;
while (bound < length(OneAtom))
    bound = bound+1;
    if (length(ZeroAtoms(ZeroAtoms > OneAtom(bound))) < bound)
        boundary = OneAtom(bound-1)-1;
        break
    end
end
FidelityPhotons(indA,indB) = length(OneAtom(OneAtom>boundary))/length(OneAtom);
title(strcat('Treshold>',int2str(boundary),' Fidelity:',num2str(FidelityPhotons(indA,indB))));
hold off

ZeroAtoms = sort(ElectronCount(AtomCount < 1));
OneAtom = sort(ElectronCount(AtomCount > 0));
bins = [0:EMGain:max(OneAtom)*1.05];
figure(2)
% [n1,xout1] = hist(ZeroAtoms,(bins(1:end-1)+((bins(2)-bins(1))/2)));
% [n2,xout2] = hist(OneAtom,(bins(1:end-1)+((bins(2)-bins(1))/2)));
% bar(xout1,n1,'b');
% hold on
% bar(xout2,n2,'r');
% hold off
histogram(ZeroAtoms,bins);
hold on
histogram(OneAtom,bins);

boundary = (EMGain*(photons+2*AddNoise*pixelpersite^2-1)/2);
%Boundary search
bound = 1;
while (bound < length(OneAtom))
    bound = bound+1;
    if (length(ZeroAtoms(ZeroAtoms > OneAtom(bound))) < bound)
        boundary = OneAtom(bound-1)-1;
        break
    end
end
%boundary = (1*(photons+2*AddNoise*pixelpersite^2)/2);
FidelityElectrons(indA,indB) = length(OneAtom(OneAtom>boundary))/length(OneAtom);
title(strcat(num2str(paramA(indA)),' ',num2str(paramB(indB)),'Treshold>',int2str(boundary),' Fidelity:',num2str(FidelityElectrons(indA,indB))));
hold off

ZeroAtoms = sort(ElectronCountDeconv(AtomCount < 1));
OneAtom = sort(ElectronCountDeconv(AtomCount > 0));
bins = linspace(0,max(OneAtom),100);
figure(3)
% [n1,xout1] = hist(ZeroAtoms,(bins(1:end-1)+((bins(2)-bins(1))/2)));
% [n2,xout2] = hist(OneAtom,(bins(1:end-1)+((bins(2)-bins(1))/2)));
% bar(xout1,n1,'b');
% hold on
% bar(xout2,n2,'r');
% hold off
histogram(ZeroAtoms,bins);
hold on
histogram(OneAtom,bins);

%Boundary search
bound = 1;
while (bound < length(OneAtom))
    bound = bound+1;
    if (length(ZeroAtoms(ZeroAtoms > OneAtom(bound))) < bound)
        boundary = OneAtom(bound-1)-1;
        break
    end
end
FidelityElectronsDeconv(indA,indB) = length(find(OneAtom>boundary))/length(OneAtom);
title(strcat(num2str(paramA(indA)),' ',num2str(paramB(indB)),'Treshold>',int2str(boundary),' Fidelity:',num2str(FidelityElectronsDeconv(indA,indB))));
hold off

drawnow
end
end
% 
figure(4)
surf(paramB,paramA,FidelityPhotons,'FaceColor','interp');
view(0,90);
figure(5)
surf(paramB,paramA,FidelityElectrons,'FaceColor','interp');
view(0,90);
figure(6)
surf(paramB,paramA,FidelityElectronsDeconv,'FaceColor','interp');
view(0,90);
%figure(5)
%[decovpic,P1] = deconvblind(CCD,PSFFit);
%surf(abs(decovpic)/max(max(abs(decovpic)))*63);
%figure(4)
%image(abs(P1)/max(max(abs(decovpic)))*63);
% figure(5)
% surf(paramB,paramA,1-FidelityElectrons,'FaceColor','interp');
% view(0,90);
