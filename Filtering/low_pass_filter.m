function [f_filt,F2,s2] = low_pass_filter(f,Ts,Tcutt)
% Input:
% f = matrix or vector real signal of MxN, columns are time
% Ts = sampling rate (s)
% Tcutt = desired low cutoff (s)
% Output:
% f_filt = filtered signal
% F2 = fourier coefficients with positive frequency (1/2 of spectrum)
% s2 = fourier frequencies (Hz) of + frequency
%
% Justin Rogers, Stanford EFML, 2014
%%

if size(f,2)==1 % make row vector
    f=f';
end

indx2 = isnan(f);

% find first and last good point
pt1 = min(find(~isnan(f)));
if isempty(pt1)
    pt1=1;
end
pt2 = max(find(~isnan(f)));
if isempty(pt2)
    pt1=length(f);
end

if sum(sum(indx2))>0
f= naninterp(f);
end

N = size(f,2);
M = size(f,1);
fs = 1/Ts; % sampling frequency

H = 1; % no window
% hanning window
% t=0:Ts:(N-1)*Ts;
% jj = 0:(N-1);
% for j=1:M
% H(j,:) = 0.5*(1-cos(2*pi*jj/(N-1))); 
% end

F = 1/(N)*fftshift(fft2(f.*H)); %fourier coefficients

% find frequencies corresponding to F
if mod(N,2)==1 % check if odd length
    s = -(N-1)/2*fs/N : fs/N : (N-1)/2*fs/N; % N odd
else % even
    s = -N/2*fs/N : fs/N : (N/2-1)*fs/N; % N even
end
[S,~]=meshgrid(s,zeros(size(F,1),1));
indx = find(abs(S)>1/Tcutt); % frequencies above cutoff
filter = 0*S+1;
filter(indx)=0;

F_filt = F.*filter; % remove high frequencies from F
f_filt = (1./H).*real(ifft2(ifftshift(F_filt.*N))); % transform back to real signal

% manually remove first/last few points to reduce ringing if large dataset
points_drop = ceil(0.50*Tcutt/Ts);

if pt1+points_drop<0.2*length(f)
f_filt(:,1:(pt1+points_drop))=NaN;
f_filt(:,(pt2-points_drop):end)=NaN;
end

f_filt(indx2)= nan; % add nans to original data

indx=find(s>0); % only return positive half of spectrum
s2=s(indx);
indx=find(S>0);
F2=F(indx);

end