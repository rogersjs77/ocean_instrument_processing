% function A = calculate_fft(X,nfft,fs)
%   - Falk Feddersen
%   - returns matrix A(num,nfft)  of fourier coefficients
%   - this uses 75% overlapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% notes Justin Rogers 2/27/13
% 
% Input:
% X     data input
% nfft  length of data to use for fft
% 
% Output:
% A     fourier coefficients A(num,nfft), where num = floor(4*n/nfft)-3;
% 
% 1. revised to allow row vector input for X
%

function A = calculate_fft2(X,nfft)

[n,m]=size(X);

if m>n
   X=X'
   [n,m]=size(X);
end
    
num = floor(4*n/nfft)-3;

%[n nfft num]

X = X-mean(X);
% sumXt = (X'*X)/n;

%WIN = hanningwindow(@hamming,nfft);
jj = [0:nfft-1]';
WIN = 0.5 * ( 1 - cos(2*pi*jj/(nfft-1)));


A = zeros(num,nfft);

% set it up so that SQR(|A(i,:)|^2) = sum(X^2)

varXwindtot = 0;

for i=1:num,
      istart = ceil((i-1)*(nfft/4)+1);
      istop = floor(istart+nfft-1);
      Xwind = X(istart:istop);
      lenX = length(Xwind);
      Xwind = Xwind - mean(Xwind);  % demean.   detrend?
      varXwind =( Xwind'*Xwind)/lenX;            
      Xwind = detrend(Xwind);

      varXwindtot = varXwindtot + varXwind;
      Xwind = Xwind .* WIN;
      tmp = ( Xwind'*Xwind)/lenX;
      if (tmp == 0),
	 Xwind = Xwind * 0.0;
      else
	 Xwind = Xwind*sqrt(varXwind/tmp);
      end;
      A(i,:) = fft(Xwind')/sqrt(nfft);
%       meanA2 = mean( A(i,:) .* conj(A(i,:)));
%      [meanA2 varXwind  meanA2/varXwind]
%      A(i,:) = A(i,:);    %  no longer necessary  * (varXwind / meanA2);   % make so above condition holds
end;

%meanA2 = mean( mean( A .* conj(A) ));

%disp('--here is total---')
%[sumXt varXwindtot/num  meanA2 ]
