function [f_filt,F2,s2] = low_pass_filter_cyclic(f,Ts,Tcutt)
% Input:
% f = matrix or vector real signal of MxN, columns are time
% Ts = sampling rate (s)
% Tcutt = desired low cutoff (s)
% Output:
% f_filt = filtered signal
% F2 = fourier coefficients with positive frequency (1/2 of spectrum)
% s2 = fourier frequencies (Hz) of + frequency
%
% this variant adds half the data before f, and half after f
% to improve ringing at edges of f_filt, low pass filtering signal
% don't use this function to get Fourier coefficients
% since we are messing with the original signal.
%
% Justin Rogers, Stanford EFML, 2014
%%

% f = Press(19,:);

if size(f,2)==1 % make row vector
    f=f';
end

% store nans position in f
indx2 = isnan(f);

ii = length(f);
jj = 0;%round(length(f)/10);
% f(1)=999;
% f(end)=99;

ff = [(nan+ones(1,jj)) f (nan+ones(1,jj))];
t=1:length(ff);

% store location of original f within ff
indx3 = [jj+1:jj+ii];

% find first and last good point
pt1 = min(find(~isnan(ff)));
pt2 = max(find(~isnan(ff)));

% good data in the middle
fgood = naninterp(ff(pt1:pt2));

% now fill nans with cyclic data
fff = ff;

%% pad with  neighbors

f = naninterp(fff);

%% try a cosine fit

% try a high order polynomial fit to end of data

% % fit beginnging of data
% indx = pt1:(pt1+100);
% 
% x = t(indx);
% y = ff(indx);
% 
% fo = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0,0,1E-4],...
%                'Upper',[Inf,2*pi,Inf],...
%                'StartPoint',[.1 1 1]);
% ft = fittype('a*cos(x/c-b)','options',fo);
% 
% [curve,gof] = fit(x,y,ft);
% 
% [p] = polyfit(x,y,3);
% fff(1:pt1-1) = polyval(p,t(1:pt1-1));
% 
% % fit end of data
% indx = (pt2-5):pt2;
% 
% x = t(indx);
% y = ff(indx);
% 
% [p] = polyfit(x,y,3);
% fff(pt2+1:end) = polyval(p,t(pt2+1:end));
% 
% f=fff;


%% try a high order polynomial fit to end of data

% fit beginnging of data
% indx = pt1:(pt1+5);
% 
% x = t(indx);
% y = ff(indx);
% 
% [p] = polyfit(x,y,3);
% fff(1:pt1-1) = polyval(p,t(1:pt1-1));
% 
% % fit end of data
% indx = (pt2-5):pt2;
% 
% x = t(indx);
% y = ff(indx);
% 
% [p] = polyfit(x,y,3);
% fff(pt2+1:end) = polyval(p,t(pt2+1:end));
% 
% f=fff;
%% wrap data

% start on first half
% j=pt1-1;
% k=2; %start at next point
% for i=1:j
%    fff(j) = fgood(k);
%    j=j-1;
%    k=k+1;    
% end
% 
% % do second half
% j=pt2+1;
% k=length(fgood)-1; %start at second to last point
% for i=1:length(ff)-j
%    fff(j) = fgood(k);
%    j=j+1;
%    k=k-1;    
% end

% f = fff;

%% this is the original fft routine
if sum(sum(indx2))>0
f= naninterp(f);
end

% store original size


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

% manually remove first/last few points to reduce ringing
f_filt(:,1:4)=NaN;
f_filt(:,end-4:end)=NaN;



indx=find(s>0); % only return positive half of spectrum
s2=s(indx);
indx=find(S>0);
F2=F(indx);

%% return original data
f_filt = f_filt(indx3);
% add nans to original data
f_filt(indx2)= nan; 

end