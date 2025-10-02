function [f_filt,F2,s2] = band_pass_filter(f,Ts,Tbot,Ttop)

% f = matrix or vector real signal of MxN, columns are time
% Ts = sampling rate (s)
% Tbot = desired low period cutoff (high freq) (s)
% Ttop = desired high period cutoff (low freq) (s)

% f_filt = filtered signal
% F2 = fourier coefficients with positive frequency (1/2 of spectrum)
% s2 = fourier frequencies (Hz) of + frequency

if size(f,2)==1
    f=f';
end

N = length(f); 
fs = 1/Ts; % sampling frequency

t=0:Ts:(N-1)*Ts;
jj = 0:(N-1);
H = 1; % no window
% H = 0.5*(1-cos(2*pi*jj/(N-1))); % hanning window

F = 1/(N)*fftshift(fft2(f.*H)); %fourier coefficients

if mod(N,2)==1 % check if odd length
    s = -(N-1)/2*fs/N : fs/N : (N-1)/2*fs/N; % N odd
else % even
    s = -N/2*fs/N : fs/N : (N/2-1)*fs/N; % N even
end
[S,~]=meshgrid(s,zeros(size(F,1),1));
indx = find(abs(S)>1/Tbot | abs(S)<1/Ttop);
filter = 0*S+1;
filter(indx)=0;

F_filt = F.*filter;

f_filt = (1./H).*real(ifft2(ifftshift(F_filt.*N)));

% manually remove first/last few points ringing
f_filt(:,1:2)=NaN;
f_filt(:,end-2:end)=NaN;
f_filt(isnan(f_filt)) = interp1(find(~isnan(f_filt)), f_filt(~isnan(f_filt)),...
    find(isnan(f_filt)),'linear','extrap');

indx=find(s>0); % only return positive half of spectrum
s2=s(indx);
indx=find(S>0);
F2=F(indx);

end