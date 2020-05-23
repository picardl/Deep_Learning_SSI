function test()

recoilvel = 5.9e-3;
scatteringrate = 2*pi*30e6;
imagingpulse = 2e-6;
NA = 0.8;
CollectedPhotonsratio = (sin(asin(NA)/2)^2)*0.74*0.9;

iterations = 500;
finalpos = zeros(iterations,3);

for ind = 1:iterations
t = 0;
nrphotons(ind) = 0;
pos = [0 0 0];
vel = [0 0 0];
posx = [];
posy = [];
while (t<imagingpulse)
    %Absorption
    posx = [posx pos(1)];
    posy = [posy pos(2)];

    dt = randexp(scatteringrate);
    pos = pos+vel*dt+[0 0 -9.81]/2*dt^2;
    t = t+dt;
    %One beam config
    %velkick = [0 0 1]*recoilvel;
    
    %Two beam config
    %A = (round(rand)*2-1);
    %velkick = [0 A 0]*recoilvel;
    
    %Four beam config
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
    posx = [posx pos(1)];
    posy = [posy pos(2)];
    
    vel = vel+velkick;
    %Emmission
    dt = randexp(scatteringrate);
    pos = pos+vel*dt+[0 0 -9.81]/2*dt^2;
    t = t+dt;
    phi = 2*pi*rand;
    theta = asin((rand-0.5)*2)+pi/2;
    velkick = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(phi)]*recoilvel;
    vel = vel+velkick;
    nrphotons(ind) = nrphotons(ind)+1;
end
figure(4)
hold on
plot(posx/1e-6,posy/1e-6)
hold off
axis equal
finalpos(ind,:) = pos(1:3)/1e-6;

end

Distancesx = sort(abs(finalpos(:,1)));
RadSigmax = Distancesx(round(iterations*0.682689492137));
Distancesy = sort(abs(finalpos(:,2)));
RadSigmay = Distancesy(round(iterations*0.682689492137));

figure(1)
plot(finalpos(:,1),finalpos(:,2),'o');
axis equal
% hold on
% plot(RadSigmax*cos([0:0.01:2*pi]),RadSigmay*sin([0:0.01:2*pi]),'r');
% hold off
% axis equal
% figure(2)
% histogram(nrphotons*CollectedPhotonsratio);
% figure(3)
% histogram(finalpos(:,3));
% 
% function R = randexp(varlambda)
%     R = -1;
%     while (R<0)
%         R=-1/varlambda*log(rand);
%     end
% end

end