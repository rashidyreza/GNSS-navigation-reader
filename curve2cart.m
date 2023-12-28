%function that converts Geodetic coordinates to cartezian:
function [Rg]= curve2cart(fi,landa,h)
format long g
a=6378137;
b=6356752.314;
for i=1:size(fi,1)
N(i,:)=a^2./sqrt(a^2.*(cos(fi(i,:))).^2+b^2.*(sin(fi(i,:))).^2);
r(i,:)=[(N(i,1)+h(i,1)).*cos(fi(i,1)).*cos(landa(i,1)) , (N(i,1)+h(i,1)).*cos(fi(i,1)).*sin(landa(i,1)) , (N(i,1)*(b^2/a^2)+h(i,1)).*sin(fi(i,1))];
end
Xg=r(:,1);
Yg=r(:,2);
Zg=r(:,3);
Rg=[Xg,Yg,Zg];
end
%%end