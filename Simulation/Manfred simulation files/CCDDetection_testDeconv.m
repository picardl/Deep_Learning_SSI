clear all
%load('PSFFitER20pixel532nm2us08NA');

paramA = [10];
paramB = [532];

for indA = 1:length(paramA)
for indB = 1:length(paramB)

NA = 0.75;
CollectedPhotonsratio = (sin(asin(NA)/2)^2)*0.74*0.8; %QE of Andor and 20%loss
lambda = 401;
scatteringrate = 2*pi*30e6;
imagingpulse = 2e-6;%paramA(indA);
pixelpersite = 10;
latticespacing = 532/2;
realpixelsize = latticespacing/pixelpersite; %in nm
sites = 7;
sensorsize = pixelpersite*sites;
PSFwidth = (0.43*lambda/2/NA)/latticespacing*pixelpersite;
wavepacketwidth = 1/10*pixelpersite;
EMGain = 1000;
AddNoise = 0.005; %ixon897
Lossrate = 1*1e-4;
photons = scatteringrate*imagingpulse/2*CollectedPhotonsratio*exp(-scatteringrate*imagingpulse/2*Lossrate);
[x,y] = ndgrid([-pixelpersite*1.45:pixelpersite*1.45],[-pixelpersite*1.45:pixelpersite*1.45]);
%PSF = Gauss2D(0,0,PSFwidth,1,x,y);

name = strcat('PSF10pixel532nm075NA',int2str(imagingpulse*1e7),'us');
load(name,'PSFFit');
PSF = PSFFit;

u = 1.66053892e-27;
m = 168*u;
c = 299792458;
hbar = 1.05457173e-34;
recoilvel = 1*5.9e-3;

iterations = 1;
nrphotonsdistdetected = zeros(iterations,1);
nrphotonsdistscattered = zeros(iterations,1);
PhotonCount = [];
ElectronCount = [];
ElectronCountDeconv = [];
AtomCount = [];
tic
for iter = 1:iterations

CCD = zeros(sensorsize,sensorsize);
CCDreconstruction = zeros(sensorsize,sensorsize);
%pattern = round(rand(sites,sites));
pattern = zeros(sites,sites);
pattern(2:end-1,2:end-1) = round(rand(sites-2,sites-2));

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
CCDPhotons = CCD;
%noise
for ind = 1:(sensorsize^2)
    CCD(ind) = CCD(ind)+randpoisson(AddNoise);
    if (CCD(ind)>0)
        CCD(ind) = randEMGain(CCD(ind),EMGain);
    end
end

CCDdeconv = deconvlucy(CCD, PSF);

for ind1 = 1:sites
    for ind2 = 1:sites
        PhotonCount = [PhotonCount sum(sum(CCDPhotons(((ind1-1)*pixelpersite+1):(ind1*pixelpersite),((ind2-1)*pixelpersite+1):(ind2*pixelpersite))))];
        ElectronCount = [ElectronCount sum(sum(CCD(((ind1-1)*pixelpersite+1):(ind1*pixelpersite),((ind2-1)*pixelpersite+1):(ind2*pixelpersite))))];
        ElectronCountDeconv = [ElectronCountDeconv sum(sum(CCDdeconv(((ind1-1)*pixelpersite+1):(ind1*pixelpersite),((ind2-1)*pixelpersite+1):(ind2*pixelpersite))))];
        AtomCount = [AtomCount pattern(ind1,ind2)];
    end
end
             
%single picture analysis
figure(1)
image(CCD/max(max(CCD))*63);
view(0,90);

figure(2)
ElectronCountimg = (reshape(ElectronCount,sites,sites))';
b = bar3(ElectronCountimg);
view(0,90);
for i = 1:length(b)
    zdata = get(b(i),'ZData');
    set(b(i),'CData',zdata)
    set(b,'EdgeColor','k') 
end

figure(3)
image(pattern*63)

figure(4)
image(CCDdeconv/max(max(CCDdeconv))*63);
view(0,90);

figure(5)
ElectronCountDeconvimg = (reshape(ElectronCountDeconv,sites,sites))';
b = bar3(ElectronCountDeconvimg);
view(0,90);
for i = 1:length(b)
    zdata = get(b(i),'ZData');
    set(b(i),'CData',zdata)
    set(b,'EdgeColor','k') 
end


drawnow

% patternreconstruct = zeros(sites,sites);
% patternreconstruct(ElectronCountimg>(max(max(ElectronCountimg))/2)) = +1;
% offset = 4;
% for ind1 = 2:sites-1
%     for ind2 = 2:sites-1
%         CCDreconstruction(((ind1-2)*pixelpersite+1):((ind1+1)*pixelpersite),((ind2-2)*pixelpersite+1):((ind2+1)*pixelpersite)) = ...
%         CCDreconstruction(((ind1-2)*pixelpersite+1):((ind1+1)*pixelpersite),((ind2-2)*pixelpersite+1):((ind2+1)*pixelpersite)) + patternreconstruct(ind1,ind2)*PSF;
%     end
% end
% probability = sum(sum(CCDreconstruction.*(CCD-offset)));
% 
% figure(6)
% image(CCDreconstruction/max(max(CCDreconstruction))*63);
% view(0,90);
% title(strcat('prob:',num2str(probability)));
% 
% changes = 1;
% while ( changes > 0);
%     changes = 0;
%     for ind1 = 2:sites-1
%         for ind2 = 2:sites-1
%             if (patternreconstruct(ind1,ind2) > 0)
%                 CCDreconstruction(((ind1-2)*pixelpersite+1):((ind1+1)*pixelpersite),((ind2-2)*pixelpersite+1):((ind2+1)*pixelpersite)) = ...
%                 CCDreconstruction(((ind1-2)*pixelpersite+1):((ind1+1)*pixelpersite),((ind2-2)*pixelpersite+1):((ind2+1)*pixelpersite)) - PSF;
% 
%             else
%                 CCDreconstruction(((ind1-2)*pixelpersite+1):((ind1+1)*pixelpersite),((ind2-2)*pixelpersite+1):((ind2+1)*pixelpersite)) = ...
%                 CCDreconstruction(((ind1-2)*pixelpersite+1):((ind1+1)*pixelpersite),((ind2-2)*pixelpersite+1):((ind2+1)*pixelpersite)) + PSF;
%             end
%             tempprob = sum(sum(CCDreconstruction.*CCD-offset));
%             figure(6)
%             image(CCDreconstruction/max(max(CCDreconstruction))*63);
%             view(0,90);
%             title(strcat('max:',num2str(probability),' prob:',num2str(tempprob),' ind1:',int2str(ind1),' ind2:',int2str(ind2)));
%             drawnow
%             pause(0.1);
%             if (tempprob > probability)
%                     patternreconstruct(ind1,ind2) = abs(patternreconstruct(ind1,ind2)-1);
%                     probability = tempprob;
%                     changes = changes+1;
%             else
%                 if (patternreconstruct(ind1,ind2) > 0)
%                     CCDreconstruction(((ind1-2)*pixelpersite+1):((ind1+1)*pixelpersite),((ind2-2)*pixelpersite+1):((ind2+1)*pixelpersite)) = ...
%                     CCDreconstruction(((ind1-2)*pixelpersite+1):((ind1+1)*pixelpersite),((ind2-2)*pixelpersite+1):((ind2+1)*pixelpersite)) + PSF;
% 
%                 else
%                     CCDreconstruction(((ind1-2)*pixelpersite+1):((ind1+1)*pixelpersite),((ind2-2)*pixelpersite+1):((ind2+1)*pixelpersite)) = ...
%                     CCDreconstruction(((ind1-2)*pixelpersite+1):((ind1+1)*pixelpersite),((ind2-2)*pixelpersite+1):((ind2+1)*pixelpersite)) - PSF;
%                 end                
%             end
%         end
%     end
% end
% 
% figure(6)
% image(CCDreconstruction/max(max(CCDreconstruction))*63);
% view(0,90);
% title(strcat('prob:',num2str(probability)));
% 

end
toc
% figure(4)
% histogram(nrphotonsdistdetected);
% hold on
% histogram(nrphotonsdistscattered);
% hold off
% 
% ZeroAtoms = sort(PhotonCount(AtomCount < 1));
% OneAtom = sort(PhotonCount(AtomCount > 0));
% %bins = [0:EMGain:max(OneAtom)*1.05];
% bins = [0:max(OneAtom)+2];
% figure(3)
% histogram(ZeroAtoms,bins);
% hold on
% histogram(OneAtom,bins);
% %boundary = (EMGain*(photons+2*AddNoise*pixelpersite^2)/2);
% boundary = (1*(photons+2*AddNoise*pixelpersite^2)/2);
% %Boundary search
% bound = 0;
% while (bound < length(OneAtom))
%     bound = bound+1;
%     
%     if (length(ZeroAtoms(ZeroAtoms > OneAtom(bound))) < bound)
%         boundary = OneAtom(bound-1);
%         break
%     end
% end
% FidelityPhotons(indA,indB) = length(OneAtom(OneAtom>boundary))/length(OneAtom);
% title(strcat('Treshold>',int2str(boundary),' Fidelity:',num2str(FidelityPhotons(indA,indB))));
% hold off

% ZeroAtoms = sort(ElectronCount(AtomCount < 1));
% OneAtom = sort(ElectronCount(AtomCount > 0));
% bins = [0:EMGain:max(OneAtom)*1.05];
% %bins = [0:max(OneAtom)+2];
% figure(4)
% histogram(ZeroAtoms,bins);
% hold on
% histogram(OneAtom,bins);
% boundary = (EMGain*(photons+2*AddNoise*pixelpersite^2-1)/2);
% %Boundary search
% bound = 0;
% while (bound < length(OneAtom))
%     bound = bound+1;
%     if (length(ZeroAtoms(ZeroAtoms > OneAtom(bound))) < bound)
%         boundary = OneAtom(bound-1);
%         break
%     end
% end
% %boundary = (1*(photons+2*AddNoise*pixelpersite^2)/2);
% FidelityElectrons(indA,indB) = length(OneAtom(OneAtom>boundary))/length(OneAtom);
% title(strcat(num2str(paramA(indA)),' ',num2str(paramB(indB)),'Treshold>',int2str(boundary),' Fidelity:',num2str(FidelityElectrons(indA,indB))));
% hold off

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
% figure(5)
% surf(paramB,paramA,FidelityPhotons,'FaceColor','interp');
% view(0,90);
% figure(6)
% surf(paramB,paramA,FidelityElectrons,'FaceColor','interp');
% view(0,90);
%figure(5)
%[decovpic,P1] = deconvblind(CCD,PSFFit);
%surf(abs(decovpic)/max(max(abs(decovpic)))*63);
%figure(4)
%image(abs(P1)/max(max(abs(decovpic)))*63);
% figure(5)
% surf(paramB,paramA,1-FidelityElectrons,'FaceColor','interp');
% view(0,90);
