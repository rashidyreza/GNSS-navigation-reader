function [Date , Ephemerides] = RinexData(filename)
filename = uigetfile('*.18n');
fileID = fopen(filename,'r');
formatSpec = '%f';
endofheader = 0;
while endofheader==0
    currentline = fgetl(fileID);
    if strfind(currentline , 'END OF HEADER')
        endofheader = 1;
    end
end
data = textscan(fileID , formatSpec);
m = size(data{1,1},1);
PRN = data{1,1}(1:38:m , 1);
Year = data{1,1}(2:38:m , 1);
Month = data{1,1}(3:38:m , 1);
Day = data{1,1}(4:38:m , 1);
Hour = data{1,1}(5:38:m , 1);
Minute = data{1,1}(6:38:m , 1);
Date = [Year , Month , Day , Hour , Minute];
Crc = data{1,1}(24:38:m , 1);
Crs = data{1,1}(12:38:m , 1);
Cic = data{1,1}(20:38:m , 1);
Cis = data{1,1}(22:38:m , 1);
Cuc = data{1,1}(15:38:m , 1);
Cus = data{1,1}(17:38:m , 1);
Omega0 = data{1,1}(21:38:m , 1);
Omegadot = data{1,1}(26:38:m , 1);
omega = data{1,1}(25:38:m , 1);
i0 = data{1,1}(23:38:m , 1);
idot = data{1,1}(27:38:m , 1);
t0 = data{1,1}(19:38:m , 1);
deltan = data{1,1}(13:38:m , 1);
M0 = data{1,1}(14:38:m , 1);
a = data{1,1}(18:38:m , 1);
e = data{1,1}(16:38:m , 1);
Ephemerides = [ PRN , Crs , Crc , Cuc , Cus , Cic , Cis , deltan , M0 , a , e , Omega0 , Omegadot , i0 , idot , omega , t0];
end
%%end