clear 
clc
[Data , Ephemerides] = RinexData('*.18n');
%
% [RSatellite , RS , RGeodetic , Vsatellite , VS , Accel] = GroundTrack(Date , Ephemerides);
FI = deg2rad(35.696);
LANDA = deg2rad(51.423);
HEIGHT = 1590;
[Azimuth ,  Elevation] = PolarPlot(Ephemerides , FI , LANDA , HEIGHT)
