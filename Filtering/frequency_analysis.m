
function [C]= frequency_analysis(X,Y,mtime,mtime_filt,nfft)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   - Falk Feddersen
%   - returns Spectrum Sm so that   2*nansum(Sm)*df = var(U)
%     as required by Trowbridge appendix for estimating dissipation 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Modifications to script
% Justin Rogers, Stanford EFML
% 1/5/2015
%
% INPUT:
% x         input signal
% V         output signal
% nfft      points to use for FFT, recommend power of 2
%            =5*60*2/dt_sec;% approximate 5 minute window for ig waves
%            =1*60*2/dt_sec, appoximate 1 minute window for swell waves
% fs        [Hz] sampling frequency
%
% OUTPUT: C structure file
% C.Sxx(f) power spectral density of x
% C.Syy(f) power spectral density of y
% C.Sxy(f) power spectral covariance of x & y
% C.coherence; squared cohereence of x&y
% C.phase; phase of xy in degrees
% C.cospectra =  cospectrum of xy
% C.quadrature = quadrature of xy
% C.f = fm; frequency vector
% C.dof = DOF; degrees of freedom
% C.sig = SIG; 95% significance level
% C.df = df; [Hz] frequency resolution is:  = fs/(nfft-1); fs is sampling frequency
% C.num_avg = num_avg; number of ffts used for averaging = floor(4*n/nfft)-3; n is #
%           of points in series

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = naninterp(X);
Y = naninterp(Y);

fs = 1./(24*(mtime(2)-mtime(1)));% do everything in hours
dt_big_day = ((mtime_filt(2)-mtime_filt(1)));
    
 for i=1:(length(mtime_filt)-2) % only process points with full time series

    tdown = mtime_filt(i+1)-dt_big_day/2;
    tup = mtime_filt(i+1)+dt_big_day/2;
    indx = find(mtime>tdown & mtime<=tup);
        
    if ~isempty(indx)
    x = X(indx);
    y = Y(indx);
    



% check for vector direction, need to be colunm vectors

if size(x,1)<size(x,2)
    x=x';
end
if size(y,1)<size(y,2)
    y=y';
end

% interpolate oxt nans
x = naninterp(x);
y = naninterp(y);

%%%%%%%%%%%%%%%%%%%%%%

Amx = calculate_fft2(x,nfft);    % this is w/ 75% overlapping
Amy = calculate_fft2(y,nfft);

[nA,mA] = size(Amx);

df = fs/(nfft-1);   % This is the frequency resolution
nnyq = nfft/2 +1;
num_avg = nA;

fm = [0:nnyq-1]*df;

% are these G is half spectrum, S is full spectrum
Sxx = nanmean( Amx .* conj(Amx) ) / (nnyq * df);  % need to normalize by freq resolution 
Gxx = Sxx(1:nnyq);

Syy = nanmean( Amy .* conj(Amy) ) / (nnyq * df);  % need to normalize by freq resolution 
Gyy = Syy(1:nnyq);

Sxy = nanmean( Amx .* conj(Amy) ) / (nnyq * df);  % need to normalize by freq resolution 
Gxy = Sxy(1:nnyq);


% COHERENCE & PHASE, etc -------------------------------------------------------
% Definitions from Elgar, 1985, JFM ... gives same result as Wunsch's def.

% Cospectrum & Quadrature:

coXY = real(Gxy);   quXY = imag(Gxy);

% Coherence & Phase at each freq-band
% *** note that it's important to calc this AFTER all merging and ensemble avg.

cohXY = (coXY.^2 + quXY.^2)./(Gxx .* Gyy);
% cohXY = (abs(Gxy)).^2./(Gxx .* Gyy); % this is equivalent

phXY  = (180/pi) .* atan2( quXY , coXY ); % degrees  



% DEGREES OF FREEDOM and level of no significant coherence --------------------
% DOF = 2 * (# independent windows) * (# bands merged)
merge = 1;
DOF = 2 * nA * merge;  
% 95% significance level for zero coherence
SIG = sqrt(6/DOF);

%% try multi-input cross spectral analysis
% emery and thomson

transf_XY = Gxy ./ Gxx;
transf_aXY = sqrt(coXY.^2 + quXY.^2)./(Gxx ); % transfer function ampl

% cohXY = (abs(transf_XY)).^2 .* Gxx ./ Gyy; % this is also equivalent


%% %% SPECIFY OUTPUT TO C structure file %%%%%%%%%%%%%%%%%%%%%%%%%%%
C.Gxx(i,:) = Gxx;
C.Gyy(i,:) = Gyy;
C.Gxy(i,:) = Gxy;
% C.coherence = cohXY;
% C.transf_a = transf_aXY;
% C.phase = phXY;
% C.cospectra = coXY;
% C.quadrature = quXY;

    end
 end
C.f = fm;
C.dof = DOF;
C.sig = SIG;
C.df = df;
C.num_avg = num_avg;
C.mtime = mtime_filt(2:end-1);
end
