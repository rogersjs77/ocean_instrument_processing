function [eta_t,ETAS_filt,s] = pressure_to_eta (P,dt,gage_ht,rho,fscutt)
% convert pressure data to free surface
% Justin Rogers
% Stanford EFML, 2012

% P = pressure time series (Pa)
% dt = sampling rate (s)
% gage_ht = height of pressure gage above bottom (m)
% rho = density

% eta_t = free surface time series
% ETAS_filt = fourier coefficients of recreated spectrum
% s = frequencies

%%
g = 9.81;

if size(P,1)==1
else % make sure P is columns
    P=P';
end

Pd=P-mean(P); % dynamic pressure

z = -mean(P)/(rho*g);
h = gage_ht-z;

N = length(Pd); 
fs = 1/dt; % sampling frequency
% fscutt=0.25;%fs/4; 

H = 1;%0.5*(1-cos(2*pi*jj/(N-1))); % could add hamming filter later

PS = 1/(N)*fftshift(fft2(Pd.*H)); %fourier coefficients of pressure

if mod(N,2)==1 % even length
    s = -(N-1)/2*fs/N : fs/N : (N-1)/2*fs/N; % N odd
else % odd length
    s = -N/2*fs/N : fs/N : (N/2-1)*fs/N; % N even
end

ks=0*s;
for i=1:length(s) % wave numbers
    [ks(i),~]=dispersion_approx(h,1/abs(s(i)));
% ks(i) = get_wavenumber(2*pi*abs(s(i)),h);    
end

kps = cosh(ks*(h+z))./cosh(ks*h); % correction factor

ETAS = PS./(rho*g*kps); % fourier coeffcients of eta

indx = find(abs(s)>fscutt);
filter = 0*s+1;
filter(indx)=0;

ETAS_filt = ETAS.*filter; % remove high frequencies

eta_t = real(ifft2(ifftshift(ETAS_filt.*N))); % transform to time space
