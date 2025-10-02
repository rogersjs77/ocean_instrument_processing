function    [WaveStats] = wave_stats_burst(u,v,p,Tm,Tstd,doffp,rho)
g = 9.81;

p_prime = p - nanmean(p);
u_prime = u - nanmean(u);
v_prime = v - nanmean(v);
h = nanmean(p)./(g*rho);

T = Tm +[-Tstd, 0, Tstd];
om = 2*pi./T;
k = get_wavenumber(2*pi./T,h);
Kp = cosh(k.*doffp)./cosh(k.*h);

prms = sqrt(nanmean(p_prime.^2));

H = sqrt(8).*prms./(Kp.*rho*g);

Hm = H(2);
Hstd = nanmean(abs(Hm-[H(1), H(3)]));

% find covariance
Cuu = nanvar(u_prime.*u_prime);
Cuv = nanvar(u_prime.*v_prime);
Cvv = nanvar(v_prime.*v_prime);
Cpu = nanvar(p_prime.*u_prime);
Cpv = nanvar(p_prime.*v_prime);
Cpp = nanvar(p_prime.*p_prime);

a1 = Cpu./sqrt(Cpp.*(Cuu+Cvv));
b1 = Cpv./sqrt(Cpp.*(Cuu+Cvv));
a2 = (Cuu-Cvv)./(Cuu+Cvv);
b2 = 2*Cuv./(Cuu+Cvv);

dir1 = rad2deg ( atan2(b1,a1) );  
spread1 = 2*(1-(a1*cos(deg2rad(dir1))+b1*sin(deg2rad(dir1))));

% this is probably more reilable
dir2 = rad2deg ( atan2(b2,a2)/2 ); 
spread2 = (180/pi).*sqrt(abs(0.5 - ...
    0.5.*( a2.*cos(2.*deg2rad(dir2)) + b2.*sin(2.*deg2rad(dir2)))));

WaveStats.depth = h;
WaveStats.Hsig_ss = Hm;
WaveStats.Hsig_ss_std = Hstd;
% WaveStats.dir_ss2 = theta;
WaveStats.dir_ss = dir2;
WaveStats.spread_ss = spread2;
WaveStats.Tm_ss = Tm;
WaveStats.Tm_ss_std = Tstd;
WaveStats.Energy_ss = (1/8)*rho*g*Hm.^2;
WaveStats.Us = g.*k(2).*Hm.*cos(deg2rad(dir2))./(8*h*om(2));
WaveStats.Vs = g.*k(2).*Hm.*sin(deg2rad(dir2))./(8*h*om(2));
end