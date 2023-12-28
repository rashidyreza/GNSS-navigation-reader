%function that converts cartezian coordinates to Geodetic:
function[fi , landa , h]=cart2curve(X , Y , Z)
format long g
a=6378137;
b=6356752.314;
e=(sqrt(a^2-b^2))/a;
P= sqrt(X.^2 + Y.^2);
fi= atan(Z./(P.*(1-e^2)));
fi= rad2deg(fi);
d=0;
dfi=1;
while dfi >1.0e-12
fi1=fi;
 N=a^2./sqrt(a^2.*(cos(fi1)).^2+b^2.*(sin(fi1)).^2);
h= (P./cos(fi1))-N;
 fi=atan(Z./(P.*(1-(e^2.*N)./(N+h))));
dfi=norm(fi1-fi);
d=d+1;
end
% dfi=rad2deg(dfi)
% d
landa =2* atan(Y./(X+P));
landa=rad2deg(landa);
N=a^2./sqrt(a^2.*(cos(fi)).^2+b^2.*(sin(fi)).^2);
h=(P./cos(fi))-N;
fi=rad2deg(fi);
end
%%end