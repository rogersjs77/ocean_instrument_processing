function [k,om]=dispersion_approx(h,T)
g = 9.81;
om = 2*pi/T;
% find solution to k using Hunt (1979), Dean pg 72
 yh = (om.^2.*h)/g;                                    
d1 = 0.6666666666*yh;
d2 = 0.3555555555*yh.^2;
d3 = 0.1608465608*yh.^3;
d4 = 0.0632098765*yh.^4;
d5 = 0.0217540484*yh.^5;
d6 = 0.0065407983*yh.^6;
sumd = d1+d2+d3+d4+d5+d6;
dis = yh.^2+(yh./(1+sumd));
k = sqrt(dis)./h;