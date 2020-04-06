%Set imaging simulation paramters

pixelspersite = 20;

NA = 0.9;
CollectedPhotonsratio = (sin(asin(NA)/2)^2)*0.74*0.8; %QE of Andor and 20%loss
lambda = 401;
scatteringrate = 2*pi*30e6;
imagingpulse = 3e-6;
pixelpersite = 10;
latticespacing = 532/2; %in nm
realpixelsize = latticespacing/pixelpersite; %in nm
sites = 20;
sensorsize = pixelpersite*sites;
PSFwidth = (0.43 * lambda * pixelpersite) / (2 * NA * latticespacing);
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
init_angle = 0; %(0.2/pi)*180;

PhotonCount = [];
ElectronCount = [];
AtomCount = [];

