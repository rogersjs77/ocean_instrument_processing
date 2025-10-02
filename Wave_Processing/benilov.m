
function [B]= benilov(U,V,W,P,nfft,doffp,fs,fcutoff,rho)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Benilov and Fuyiskin, 1970 method
% for Reynolds Stress computation with strong waves
% Justin Rogers, Stanford EFML
% 6/2014
%
% INPUT:
% U         [m/s] u velocity signal
% V         [m/s] v velocity signal
% P         [dBar] pressure signal
% nfft      points to use for FFT, recommend power of 2
%            =5*60/dt_sec;% approximate 5 minute window for ig waves
%            =1*60/dt_sec, appoximate 1 minute window for swell waves
% doffp     [m] height of p measurement off bottom
% fs        [Hz] sampling frequency
% fcutoff   [Hz] high frequency cutoff, depends on fs and instrument depth
% rho       [kg/m^3] density

%
% OUTPUT: WaveStats structure file
% Structure file B with following variables
% B.See = Seet;
% B.Suu = Suut;
% B.Svv = Svvt;
% B.Sww = Swwt;
% B.Suv = Suvt;
% B.Suw = Suwt;
% B.Svw = Suvt;
% B.Sue = Suet;
% B.Sve = Svet;
% B.Swe = Swet;
% 
% B.S_uwave_wwave = S_uwave_wwave;
% B.S_vwave_wwave = S_vwave_wwave;
% B.Supwp = Supwp;
% B.Svpwp = Svpwp;
% B.upwp_bar = upwp_bar;
% B.vpwp_Bar = vpwp_bar;
% % fmt       (f) [Hz] frequency vector

% df        [Hz] frequency resolution is:  = fs/(nfft-1); fs is sampling frequency
% num_avg   number of ffts used for averaging = floor(4*n/nfft)-3; n is #
%           of points in series

% JSR Edits 2/26/2013
% 1. added fs as an input
% 2. changed some notation
% 3. allow for row or column vector input
% 4. inlcluded more output variables
% 5. modified ouput to single structure file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%
g = 9.81;    % gravity, m/s^2
% rho = 1028;    % density (Kg /m^-3) 

if (all(isnan(U)) | all(isnan(V)) | all(isnan(P)) ),
   dir2_ss = nan;   spread2_ss = nan;

   fcentroid_swell = nan;
   Tm_ss = 1/fcentroid_swell;
   Hsig_ss = nan; 
   Hsig_ig = nan;
   depth=nan;
   dirt = [nan nan];   spreadt = nan;
   dir_ss = nan;   spread_ss = nan;   
   a1t = nan;    a2t = nan;    b1t = nan;    b2t = nan;
   a1_ss = nan;    a2_ss = nan;    b1_ss = nan;    b2_ss = nan;
   Sxx_ss = nan;  Sxy_ss=nan;  Eflux_ss = [nan nan nan ];
   ztest_ss = nan;  Cpu_ss = nan;
   SSEt = nan;
   Suut = nan;   Svvt = nan;  Suvt = nan;  dir2t=nan; spread2t=nan; fmt=nan;
   return;
   % finish this
end;

% check for vector direction, need to be colunm vectors
if size(P,1)<size(P,2)
    P=P';
    U=U';
    V=V';
    W=W';
end
% interpolate out nans
P = naninterp(P);
U = naninterp(U);
V = naninterp(V);
W = naninterp(W);

depth = nanmean(P)*1E4/(rho*g) + doffp;  % This is the averaged water depth

if (depth <= 0.0),
   disp('** Warning!  depth < 0.  Abort')
   dir2_ss = nan;   spread2_ss = nan;   
   fcentroid_swell = nan;
   Tm_ss = 1/fcentroid_swell;
   fspread = nan;
   Hsig_ss = nan; 
   Hsig_ig = nan; 
   depth=nan;
   ztest_ss = nan;  Cpu_ss = nan;   
   dirt = [nan nan];   spreadt = nan;
   dir_ss = nan;   spread_ss = nan;      
   a1t = nan;    a2t = nan;    b1t = nan;    b2t = nan;
   a1_ss = nan;    a2_ss = nan;    b1_ss = nan;    b2_ss = nan;
   Sxx_ss = nan;  Sxy_ss=nan;  Eflux_ss = [nan nan nan ];
   SSEt = nan;
   Suut = nan;   Svvt = nan;  Suvt = nan;  dir2t=nan; spread2t=nan; fmt=nan;
   return;
end;   

%%%%%%%%%%%%%%%%%%%%%%

Amu = calculate_fft2(U,nfft);    % this is w/ 75% overlapping
Amv = calculate_fft2(V,nfft);
Amw = calculate_fft2(W,nfft);
Amp = calculate_fft2(P,nfft);

[nA,mA] = size(Amu);

df = fs/(nfft-1);   % This is the frequency resolution
nnyq = nfft/2 +1;
num_avg = nA;

fm = [0:nnyq-1]*df;

Suu = nanmean( Amu .* conj(Amu) ) / (nnyq * df);  % need to normalize by freq resolution 
Suu = Suu(1:nnyq);

Svv = nanmean( Amv .* conj(Amv) ) / (nnyq * df);  % need to normalize by freq resolution 
Svv = Svv(1:nnyq);

Sww = nanmean( Amw .* conj(Amw) ) / (nnyq * df);  % need to normalize by freq resolution 
Sww = Sww(1:nnyq);

Suv = nanmean( Amu .* conj(Amv) ) / (nnyq * df);  % need to normalize by freq resolution 
Suv = Suv(1:nnyq);

Suw = nanmean( Amu .* conj(Amw) ) / (nnyq * df);  % need to normalize by freq resolution 
Suw = Suw(1:nnyq);

Svw = nanmean( Amv .* conj(Amw) ) / (nnyq * df);  % need to normalize by freq resolution 
Svw = Svw(1:nnyq);

Spp = nanmean( Amp .* conj(Amp) ) / (nnyq * df);  % need to normalize by freq resolution 
Spp = Spp(1:nnyq);

Spu = nanmean( Amp .* conj(Amu) ) / (nnyq * df);  % need to normalize by freq resolution 
Spu = Spu(1:nnyq);

Spv = nanmean( Amp .* conj(Amv) ) / (nnyq * df);  % need to normalize by freq resolution 
Spv = Spv(1:nnyq);

Spw = nanmean( Amp .* conj(Amw) ) / (nnyq * df);  % need to normalize by freq resolution 
Spw = Spw(1:nnyq);

% The result is now that  nansum(Suu)*df = var(U) - actually < var(U) because of detrending

% DEPTH CORRECTION & SPECTRAL WEIGHTED AVERAGES / STATS--------------------
% find correction/conversion rcoefs at each f-band 
% to calc sea surface elevation and convert velocities to pressure units 


omega = 2*pi*fm;   % This is radian frequency

warning off MATLAB:divideByZero
k = get_wavenumber( omega , depth);
warning on MATLAB:divideByZero


correction = zeros(1,nnyq);
% convert = zeros(1,nnyq);

ii = find(fm <= fcutoff);
correction(ii) = 1E4 *cosh(k(ii)*depth)./(rho * g * cosh(k(ii)*doffp));
% convert(ii) = (2*pi*fm(ii)./(g*k(ii))) .* cosh(k(ii)*doffp) ./ cosh(k(ii)*doffu);  % record conversion coefs to check them

See = Spp .* (correction.^2) ; % correct pressure for attentuation to get sea-surface elevation
Sue = Spu .* correction;
Sve = Spv .* correction;
Swe = Spw .* correction;

Seet = See(ii);
Suut = Suu(ii);
Svvt = Svv(ii);
Swwt = Sww(ii);
Suvt = Suv(ii);
Suwt = Suw(ii);
Svwt = Svw(ii);
Suet = Sue(ii);
Svet = Sve(ii);
Swet = Swe(ii);

fmt = fm(ii);

% find indices of freq bands
f_ig_low = 1/250;
f_ig_high = 1/33.33;
f_swell_low = 1/25;
f_swell_high = 1/5;

i_ig = find(fm>f_ig_low & fm<f_ig_high);
i_swell = find(fmt>f_swell_low & fmt<f_swell_high);

% Benilov and Filyushkin 1970 --------------------------------------------------------
% Su_wave_w_wave = Sue(f)(Swe(f)*)/See(f)
% Su’w’ = Suw – Su_wave_w_wave
% <u’w’> = integral(real(Su’w’(f)))df

% wave results
% auto
S_uwave_uwave = Suet .* conj(Suet) ./ Seet;
S_vwave_vwave = Svet .* conj(Svet) ./ Seet;
S_wwave_wwave = Swet .* conj(Swet) ./ Seet;
% cross
S_uwave_wwave = Suet .* conj(Swet) ./ Seet;
S_vwave_wwave = Svet .* conj(Swet) ./ Seet;
S_uwave_vwave = Suet .* conj(Svet) ./ Seet;

% turbulence results
% auto
Supup = Suut - S_uwave_uwave;
Svpvp = Svvt - S_vwave_vwave;
Swpwp = Swwt - S_wwave_wwave;
% cross
Supwp = Suwt - S_uwave_wwave;
Svpwp = Svwt - S_vwave_wwave;
Supvp = Suvt - S_uwave_vwave;

% integrate wave over spectra, only in wave band
% auto
uwave_uwave_bar = nansum( real(S_uwave_uwave(i_swell)) .* df);
vwave_vwave_bar = nansum( real(S_vwave_vwave(i_swell)) .* df);
wwave_wwave_bar = nansum( real(S_wwave_wwave(i_swell)) .* df);
% cross
uwave_wwave_bar = nansum( real(S_uwave_wwave(i_swell)) .* df); % i think this is right (coherence not quadrature)
vwave_wwave_bar = nansum( real(S_vwave_wwave(i_swell)) .* df);
uwave_vwave_bar = nansum( real(S_uwave_vwave(i_swell)) .* df);

% integrate turbulence over entire spectra
% the 2 comes from a double sided spectra
% auto
upup_bar = 2*nansum( real(Supup) .* df);
vpvp_bar = 2*nansum( real(Svpvp) .* df);
wpwp_bar = 2*nansum( real(Swpwp) .* df);
% cross
upwp_bar = 2*nansum( real(Supwp) .* df); % i think this is right (coherence not quadrature)
% upwp_bar = nansum( abs(Supwp) .* df); I think this is incorrect
vpwp_bar = 2*nansum( real(Svpwp) .* df);
upvp_bar = 2*nansum( real(Supvp) .* df);

%% condense results

B.See = Seet;
B.Suu = Suut;
B.Svv = Svvt;
B.Sww = Swwt;
B.Suv = Suvt;
B.Suw = Suwt;
B.Svw = Suvt;
B.Sue = Suet;
B.Sve = Svet;
B.Swe = Swet;

B.S_uwave_wwave = S_uwave_wwave;
B.S_vwave_wwave = S_vwave_wwave;
B.Supwp = Supwp;
B.Svpwp = Svpwp;

% wave tensor
B.uwave_uwave_bar = uwave_uwave_bar;
B.vwave_vwave_bar = vwave_vwave_bar;
B.wwave_wwave_bar = wwave_wwave_bar;

B.uwave_vwave_bar = uwave_vwave_bar;
B.uwave_wwave_bar = uwave_wwave_bar;
B.vwave_wwave_bar = vwave_wwave_bar;

% reynolds stress tensor
B.upup_bar = upup_bar;
B.vpvp_bar = vpvp_bar;
B.wpwp_bar = wpwp_bar;

B.upwp_bar = upwp_bar;
B.vpwp_bar = vpwp_bar;
B.upvp_bar = upvp_bar;

% other info
B.fm = fmt;
B.df = df;
B.nfft = nfft;
B.fcutoff = fcutoff;
B. depth = depth;
B.fs = fs;
B.num_avg = num_avg;

%%
% subplot(1,2,1)
% loglog(fmt,abs(Supwp),fmt,abs(Suwt),fmt,abs(S_uwave_wwave))
% legend('<u''w''>','<uw>','<uwave wwave>','location','ne')
% ylabel('S'), xlabel('Hz')
% 
% subplot(1,2,2)
% loglog(fmt,abs(Svpwp),fmt,abs(Svwt),fmt,abs(S_vwave_wwave))
% legend('<v''w''>','<vw>','<vwave wwave>','location','ne')
% ylabel('S'), xlabel('Hz')



end
