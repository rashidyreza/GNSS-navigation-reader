function [RSatellite , RS , RGeodetic , Vsatellite , VS , Accel] = GroundTrack(Date , Ephemerides)
%
% clear all
% clc
% [Date , Ephemerides] = RinexData('data.18n');
%
%For Whole Of Rinex Data:
% PRN = Ephemerides(:,1);
% Year = 2000+Date(:,1);
% Month = Date(:,2);
% Day = Date(:,3);
% Hour = Date(:,4);
% Minute = Date(:,5);
% Second = 0;
%
format long g
Mio = 3.986004418e+14;
V = 7.292115e-5;
pi = 3.14159265;
t = transpose(0:600:86400);
% 
figure;
axesm( 'MapProjection' , 'mercator' , 'MapLatLimit' , [-75 75]);
geoshow( 'landareas.shp' , 'FaceColor' , [0.5 1.0 0.5]);
setm( gca , 'grid' , 'on' , 'meridianlabel' , 'on' , 'mlabellocation' , 30 , 'parallellabel' , 'on' , 'plabellocation' , 15);
wgs84 = wgs84Ellipsoid('meters');
%
%For Every Suggested Hour & Minute By Us:
p = find(Date(:,4) == 0 & Date(:,5) == 0);
PRN = Ephemerides(p,1);
Year = 2000+Date(p,1);
Month = Date(p,2);
Day = Date(p,3);
Hour = Date(p,4);
Minute = Date(p,5);
Second = 0;
% 
for j = 1:size(PRN,1)
% 
Crs(j,1) = Ephemerides(j,2);
Crc(j,1) = Ephemerides(j,3);
Cuc(j,1) = Ephemerides(j,4);
Cus(j,1) = Ephemerides(j,5);
Cic(j,1) = Ephemerides(j,6);
Cis(j,1) = Ephemerides(j,7);
deltan(j,1) = Ephemerides(j,8);
M0(j,1) = Ephemerides(j,9);
a(j,1) = Ephemerides(j,10);
e(j,1) = Ephemerides(j,11);
Omega0(j,1) = Ephemerides(j,12);
Omegadot(j,1) = Ephemerides(j,13);
i0(j,1) = Ephemerides(j,14);
idot(j,1) = Ephemerides(j,15);
omega(j,1) = Ephemerides(j,16);
t0(j,1) = Ephemerides(j,17);
A(j) = a(j).^2;
n0(j) = sqrt((Mio)./((A(j)).^3));
n(j) = n0(j) + deltan(j);
%
for k = 1:size(t,1)
    %
    T(k) = t(k) - t0(j);
    M(k) = M0(j) + (n(j)).*T(k);
    Mdot(k) = n(j);
    E0(k) = M(k);
    E(k) = E0(k) - ((E0(k) - (e(j)).*sin(E0(k)) - M(k))./(1 - (e(j)).*cos(E0(k))));
    dE(k) = 1;
    loopnumber = 0;
    %
    while dE(k) > 1e-14
        E0(k) = E(k);
        E(k) = E0(k) - ((E0(k) - e(j).*sin(E0(k)) - M(k))./(1 - e(j).*cos(E0(k))));
        dE(k) = abs(E(k) - E0(k));
        loopnumber = loopnumber +1;
    end
    %
    Edot(k) = (Mdot(k))./(1 - (e(j)).*cos(E(k)));
    %
    %Teta Is True Anomaly...
    C = (cos(E(k)) - e(j))./(1 - e(j).*cos(E(k)));
    S = (sqrt(1-(e(j)).^2).*sin(E(k)))./(1 - e(j).*cos(E(k)));
    Teta(k) = atan2(S , C);
    Tetadot(k) = ((sin(E(k)))*(1 + (e(j)).*(cos(Teta(k))))*Edot(k))./((1-(e(j)).*cos(E(k)))*(sin(Teta(k))));
    %Fi Is Argument Of Latitude & omega Is Arument Of Perigee Point...
    Fi(k) = Teta(k) + omega(j);
    Fidot(k) = Tetadot(k);
    %Corrections Of Radius , Inclination Angle & Latitude Angular...
    deltaU(k) = (Cuc(j)).*cos(2*Fi(k)) + (Cus(j)).*sin(2*Fi(k));
    deltaR(k) = (Crc(j)).*cos(2*Fi(k)) + (Crs(j)).*sin(2*Fi(k));
    deltaI(k) = (Cic(j)).*cos(2*Fi(k)) + (Cis(j)).*sin(2*Fi(k));
    deltaUdot(k) = 2*(-(Cuc(j)).*sin(2*Fi(k)) + (Cus(j)).*cos(2*Fi(k)))*Fidot(k);
    deltaRdot(k) = 2*(-(Crc(j)).*sin(2*Fi(k)) + (Crs(j)).*cos(2*Fi(k)))*Fidot(k);
    deltaIdot(k) = 2*(-(Cic(j)).*sin(2*Fi(k)) + (Cis(j)).*cos(2*Fi(k)))*Fidot(k);
    %Corrected Latitude Angular...
    U(k) = Fi(k) + deltaU(k);
    Udot(k) = Fidot(k) + deltaUdot(k);
    %Corrected Radius...
    R(k) = (A(j)).*(1 - (e(j)).*cos(E(k))) + deltaR(k);
    Rdot(k) = (A(j)).*(e(j)).*(sin(E(k)))*(Edot(k)) + deltaRdot(k);
    %Corrected Inclination Angle...
    I(k) = i0(j) + (idot(j)).*T(k) + deltaI(k);
    Idot(k) = idot(j) + deltaIdot(k);
    %x & y Are Position Coordinates In Orbital Plane...
    x(k) = (R(k))*cos(U(k));
    xdot(k) = (Rdot(k))*cos(U(k)) - (R(k))*(sin(U(k)))*Udot(k);
    y(k) = (R(k))*sin(U(k));
    ydot(k) = (Rdot(k))*sin(U(k)) + (R(k))*(cos(U(k)))*Udot(k);
    %Omega Is Corrected Longitude Of Ascending Node. V Is Angular Velocity Of Earth's Rotation...
    OMEGA(k) = Omega0(j) + (Omegadot(j) - V).*T(k) - V*t0(j);
    OMEGAdot(k) = Omegadot(j) - V;
    %X ,Y ,Z Are Geocentric Satellite Coordinates In Earth Fixed Coordinate System Or ECEF...
    X(k) = (x(k))*cos(OMEGA(k)) - (y(k))*(sin(OMEGA(k)))*cos(I(k));
    Y(k) = (x(k))*sin(OMEGA(k)) + (y(k))*(cos(OMEGA(k)))*cos(I(k));
    Z(k) = (y(k))*sin(I(k));
    Xdot(k) = (xdot(k))*cos(OMEGA(k)) - (x(k))*(sin(OMEGA(k)))*OMEGAdot(k) - (ydot(k))*(cos(I(k)))*(sin(OMEGA(k))) - (y(k))*(-(sin(I(k)))*(sin(OMEGA(k)))*Idot(k) + (cos(I(k)))*(cos(OMEGA(k)))*OMEGAdot(k));
    Ydot(k) = (xdot(k))*sin(OMEGA(k)) + (x(k))*(cos(OMEGA(k)))*OMEGAdot(k) + (ydot(k))*(cos(I(k)))*(cos(OMEGA(k))) + (y(k))*(-(sin(I(k)))*(cos(OMEGA(k)))*Idot(k) - (cos(I(k)))*(sin(OMEGA(k)))*OMEGAdot(k));
    Zdot(k) = (ydot(k))*sin(I(k)) + (y(k))*(cos(I(k)))*Idot(k);
    %
    RS(k) = sqrt((X(k))^2 + (Y(k))^2 + (Z(k))^2);
    VS(k) = sqrt((Xdot(k))^2 + (Ydot(k))^2 + (Zdot(k))^2);
    %
    AccelVec(k,:) = (Mio./(RS(k)).^3).*[X(k) ; Y(k) ; Z(k)];
    %
    hVec(k,:) = cross([X(k) ; Y(k) ; Z(k)] , [Xdot(k) ; Ydot(k) ; Zdot(k)]);
end
% 
RSatellite = [X ; Y ; Z]';
RS = RS';
Vsatellite = [Xdot ; Ydot ; Zdot]';
VS = VS';
Accel = sqrt((AccelVec(:,1)).^2 + (AccelVec(:,2)).^2 + (AccelVec(:,3)).^2);
% [fi landa height] = ecef2geodetic(wgs84 , X , Y , Z);
[fi , landa , height] = cart2curve(X , Y , Z);
RGeodetic = [fi ; landa ; height]';
fileID = fopen('SatellitePositionECEF.txt','w');
fprintf(fileID , '%13s             %13s             %13s             %13s\r\n' , 'Epoc' , 'X' , 'Y' , 'Z');
for z = 1:size(t,1)
fprintf(fileID , '%13d             %13.8f             %13.8f              %13.8f\r\n' , t(z,1) , (RSatellite(z,1:3)));
end
fclose(fileID);
type 'SatellitePositionECEF.txt';
%
fileID = fopen('SatellitePositionGeodetic.txt','w');
fprintf(fileID , '%13s             %13s            %13s             %13s\r\n' , 'Epoc' , 'Fi' , 'Landa' , 'Height');
for z = 1:size(t,1)
fprintf(fileID , '%13d             %2.13f              %2.13f              %13.8f\r\n' , t(z,1) , (RGeodetic(z,1:3)));
end
fclose(fileID);
type 'SatellitePositionGeodetic.txt';
linem(fi , landa ,'linewidth' , 1 , 'color' , [rand rand rand] , 'marker' , 'o' );
%
pause
end
end
%%end