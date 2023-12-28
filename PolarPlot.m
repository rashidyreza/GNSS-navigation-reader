function [Azimuth Elevation] = PolarPlot(Ephemerides , FI , LANDA , HEIGHT)
Mio = 3.986004418e+14;
V = 7.292115e-5;
pi = 3.14159265;
figure;
axes('Visible','off','XGrid','off','YGrid','off');
PolarGrid(2 , 2);
PRN = Ephemerides(:,1);
t = transpose(0:300:86400);
for j =1:size(PRN , 1) 
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
A = a.^2;
n0 = sqrt((Mio)./(A.^3));
n = n0 + deltan;
for k = 1:size(t,1)
    T(k) = t(k) - t0(j);
    M(k) = M0(j) + n(j).*T(k);
    E0(k) = M(k);
    E(k) = E0(k) - ((E0(k) - e(j).*sin(E0(k)) - M(k))/(1 - e(j).*cos(E0(k))));
    dE(k) = 1;
    loopnumber = 0;
    while dE(k) > 1e-14
        E0(k) = E(k);
        E(k) = E0(k) - ((E0(k) - e(j).*sin(E0(k)) - M(k))/(1 - e(j).*cos(E0(k))));
        dE(k) = abs(E(k) - E0(k));
        loopnumber = loopnumber +1;
    end
    %Teta Is True Anomaly...
    C = (cos(E(k)) - e(j))/(1 - e(j).*cos(E(k)));
    S = (sqrt(1-(e(j)).^2)*sin(E(k))) / (1 - e(j).*cos(E(k)));
    if S > 0
    if C > 0
        Teta(k) = abs(atan(S/C));
    else
        Teta(k) = pi - abs(atan(S/C));
    end
    else
        if C > 0
        Teta(k) = 2*pi - abs(atan(S/C));
        else
            Teta(k) = pi + abs(atan(S/C));
        end
    end
    %Fi Is Argument Of Latitude & omega Is Arument Of Perigee Point...
    Fi(k) = Teta(k) + omega(j);
    %Corrections Of Radius , Inclination Angle & Latitude Angular...
    deltaU(k) = (Cuc(j))*cos(2*Fi(k)) + (Cus(j))*sin(2*Fi(k));
    deltaR(k) = (Crc(j))*cos(2*Fi(k)) + (Crs(j))*sin(2*Fi(k));
    deltaI(k) = (Cic(j))*cos(2*Fi(k)) + (Cis(j))*sin(2*Fi(k));
    %Corrected Latitude Angular...
    U(k) = Fi(k) + deltaU(k);
    %Corrected Radius...
    R(k) = (A(j)).*(1 - e(j).*cos(E(k))) + deltaR(k);
    %Corrected Inclination Angle...
    I(k) = i0(j) + idot(j).*T(k) + deltaI(k);
    %x & y Are Position Coordinates In Orbital Plane...
    x(k) = (R(k))*cos(U(k));
    y(k) = (R(k))*sin(U(k));
    %Omega Is Corrected Longitude Of Ascending Node. V Is Angular Velocity Of Earth's Rotation...
    Omega(k) = Omega0(j) + (Omegadot(j) - V).*T(k) - V*t0(j);
    %X ,Y ,Z Are Geocentric Satellite Coordinates In Earth Fixed Coordinate System Or ECEF...
    X(k) = (x(k))*cos(Omega(k)) - (y(k))*(sin(Omega(k)))*cos(I(k));
    Y(k) = (x(k))*sin(Omega(k)) + (y(k))*(cos(Omega(k)))*cos(I(k));
    Z(k) = (y(k))*sin(I(k));
    %Reciever Geodetic Coordinate
    a=6378137;
    b=6356752.314;
    N=a^2/sqrt(a^2*(cos(FI)).^2+b^2.*(sin(FI)).^2);
    Rg=[(N+HEIGHT)*cos(FI)*cos(LANDA) , (N+HEIGHT)*cos(FI)*sin(LANDA) , (N*(b^2/a^2)+HEIGHT)*sin(FI)];
    Xg = Rg(:,1);
    Yg = Rg(:,2);
    Zg = Rg(:,3);
    %Position Of Satellite Relative To The Reciever In Geodetic & Local Geodetic System.
    RgRelative = [X - Xg ; Y - Yg ; Zg - Z];
    P2 = [1 , 0 , 0 ; 0 , -1 , 0 ; 0 , 0 , 1];
    R2 = [cos((pi/2) - FI) , 0 , -sin((pi/2) - FI) ; 0 , 1 , 0 ; sin((pi/2) - FI) , 0 , cos((pi/2) - FI)];
    R3 = [cos(pi - LANDA) , sin(pi - LANDA) , 0 ; -sin(pi - LANDA) , cos(pi - LANDA) , 0 ; 0 , 0 , 1];
    TR=R3*R2*P2;
    RLG = inv(TR)*RgRelative;
    XLG(1,k) = (RLG(1,k))./sqrt((RLG(1,k)).^2 + (RLG(2,k)).^2 + (RLG(3,k)).^2);
    YLG(1,k) = (RLG(2,k))./sqrt((RLG(1,k)).^2 + (RLG(2,k)).^2 + (RLG(3,k)).^2);
    ZLG(1,k) = (RLG(3,k))./sqrt((RLG(1,k)).^2 + (RLG(2,k)).^2 + (RLG(3,k)).^2);
    %Azimuth Of Satellite For GroundTrack
    if YLG(k)>0
        if XLG(k)>0
            Azimuth(1,k) = abs(atan((YLG(k))./(XLG(k))));
        else
            Azimuth(1,k) = pi - abs(atan((YLG(k))./(XLG(k))));
        end
    else
        if XLG(k)>0
            Azimuth(k) = 2*pi - abs(atan((YLG(k))./(XLG(k))));
        else
            Azimuth(k) = pi + abs(atan((YLG(k))./(XLG(k))));
        end
    end
    %Heigth Of Satellite For GroundTack
    Elevation(k) = (ZLG(k))./sqrt((XLG(k)).^2 + (YLG(k)).^2);
        %Polar Coordinates Of Satellite
            if Elevation(k)>0
            Xp(k) = ((pi/2) - Elevation(k)).*sin(Azimuth(k));
            Yp(k) = ((pi/2) - Elevation(k)).*cos(Azimuth(k));
            end
%         hold on
%         plot(Xp , Yp , '*');
end
hold on
plot(Xp , Yp , '*');
end
end
%%end