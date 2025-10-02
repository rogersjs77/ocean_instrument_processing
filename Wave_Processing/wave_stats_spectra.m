
function [WaveStats]= wave_stats_spectra(U,V,P,nfft,doffu,doffp,fs,fcutoff,rho,dirmethod)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   - Falk Feddersen
%   - returns Spectrum Sm so that   2*nansum(Sm)*df = var(U)
%     as required by Trowbridge appendix for estimating dissipation 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Modifications to script
% Justin Rogers, Stanford EFML
% June 2017
%
% INPUT:
% U         [m/s] u velocity signal
% V         [m/s] v velocity signal
% P         [dBar] pressure signal
% nfft      points to use for FFT, recommend power of 2
%            =5*60*2/dt_sec;% approximate 5 minute window for ig waves
%            =1*60*2/dt_sec, appoximate 1 minute window for swell waves
% doffu     [m] height of u measurement off bottom
% doffp     [m] height of p measurement off bottom
% fs        [Hz] sampling frequency
% fcutoff   [Hz] high frequency cutoff, depends on fs and instrument depth
% rho       [kg/m^3] density
% dirmethod = 1, use a1, b1 (puv) correlations for direction, can cover full 360 deg
%           but possible errors from reflections, less accurate?
%           = 2, use a2, b2 (uv) correlations for direction, can only cover
%           +/- 90 deg from x axis. Only use this method if you know waves
%           are in + x direction, but this is likely more accurate. See
%           Herbers et al, 1999 JGR for more discussion.
%
% notes:    infragravity band = 0.004 < fm < 0.03 Hz (4.2 < 1/fm < 0.5 min) 
%           swell band = 0.05 < fm < 0.20 Hz (20 < 1/fm < 5 s)

%
% OUTPUT: WaveStats structure file
% SSEt      (f) power spectral density of surface
% Suut      (f) power spectral density of u
% Svvt      (f) power spectral density of v
% Suvt      (f) power covaiance of uv
% fmt       (f) [Hz] frequency vector
% dirt      (f) wave direction, row 1 is from 1st moment, row 2 is 2nd moment
% spreadt   (f) wave spread row 1 is from 1st moment, row 2 is 2nd moment
% a1t       (f) a1 moment 
% a2t       (f) a2 moment
% b1t       (f) b1 moment
% b2t       (f) b2 moment
% Hsig      (f) [m] significant wave height
% Hsig_ss   [m] significant wave height for swell waves
% Hsig_ig   [m] significant wave height for infragravity waves
% Tm_ss     [s] mean wave period for swell
% Eflux_ss  [?] mean energy flux for swell [posX2; posX; negX; posY2; posY; negY];
% Sxy_ss    [?] radiation stress xy for swell
% dir_ss    [deg] mean direction for swell, 0 deg is in u direction, +CCW (math notation)
% spread_ss [deg] mean wave spread for swell
% a1_ss     a1 moment for swell
% a2_ss     a2 moment for swell
% b1_ss     b1 moment for swell
% b2_ss     b2 moment for swell
% ztest_ss  zest for swell
% Cpu_ss    cross correlation of P and U in swell band
% depth     [m] mean depth
% Us        [m/s] stokes drift in u direction
% Vs        [m/s] stokes drift in v direction
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
% does this work for buried paros sensors??

% check for vector direction, need to be colunm vectors
if size(P,1)<size(P,2)
    P=P';
end
if size(U,1)<size(U,2)
    U=U';
end
if size(V,1)<size(V,2)
    V=V';
end

% interpolate out nans
P = naninterp(P);
U = naninterp(U);
V = naninterp(V);


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

Suv = nanmean( Amu .* conj(Amv) ) / (nnyq * df);  % need to normalize by freq resolution 
Suv = Suv(1:nnyq);

Spp = nanmean( Amp .* conj(Amp) ) / (nnyq * df);  % need to normalize by freq resolution 
Spp = Spp(1:nnyq);

Spu = nanmean( Amp .* conj(Amu) ) / (nnyq * df);  % need to normalize by freq resolution 
Spu = Spu(1:nnyq);

Spv = nanmean( Amp .* conj(Amv) ) / (nnyq * df);  % need to normalize by freq resolution 
Spv = Spv(1:nnyq);

% The result is now that  nansum(Suu)*df = var(U) - actually < var(U) because of detrending

% DEPTH CORRECTION & SPECTRAL WEIGHTED AVERAGES / STATS--------------------
% find correction/conversion rcoefs at each f-band 
% to calc sea surface elevation and convert velocities to pressure units 

% find indices of freq bands
f_ig_low = 1/250;
f_ig_high = 1/33.33;
f_swell_low = 1/25;
f_swell_high = 1/5;

i_ig = find(fm>f_ig_low & fm<f_ig_high);
i_swell = find(fm>f_swell_low & fm<f_swell_high);
i_all = find(fm>f_ig_low & fm<f_swell_high);

omega = 2*pi*fm;   % This is radian frequency

warning off MATLAB:divideByZero
k = get_wavenumber( omega , depth);
warning on MATLAB:divideByZero

% Now do correction for the buried paros sensors
if doffp<0
    dzp=doffp;
    doffp=0;
    G= exp(-abs(k)*dzp); %needed because paroses are buried!!!!
    ind=find(fm>.5);
    G(ind)=1;
else
    G=1;
end

Spp=Spp.*G.*G;
Spu=Spu.*G;
Spv=Spv.*G;

correction = zeros(1,nnyq);
convert = zeros(1,nnyq);

ii = find(fm <= fcutoff);
correction(ii) = 1E4 *cosh(k(ii)*depth)./(rho * g * cosh(k(ii)*doffp));
convert(ii) = (2*pi*fm(ii)./(g*k(ii))) .* cosh(k(ii)*doffp) ./ cosh(k(ii)*doffu);  % record conversion coefs to check them

SSE = Spp .* (correction.^2) ; % correct pressure for attentuation to get sea-surface elevation
SSEt = SSE(ii);

Suut = Suu(ii);
Svvt = Svv(ii);
Suvt = Suv(ii);

fmt = fm(ii);

UUpres = Suu .* (convert.^2) ;   % convert velocities into "equivalent pressure" (to compare with PP)
VVpres = Svv .* (convert.^2) ; 
UVpres = Suv .* (convert.^2) ; 

PUpres = Spu .* convert ;   % convert x-spectra as well (but conversion is not squared here)
PVpres = Spv .* convert ; 

% should check that:   PP^2 = UUpres^2 + VVpres^2
% the 'pres' velocities are now the right units to compare with P

SppU = sqrt( UUpres.^2 + VVpres.^2);
SSEU = SppU .* (correction.^2);

Ztest =  Spp(i_swell) ./  SppU(i_swell);
ztest_ss =  nansum ( (Spp(i_swell) ./  SppU(i_swell) ).* SSE(i_swell) )/...
    nansum(SSE(i_swell));

% COHERENCE & PHASE, etc -------------------------------------------------------
% Definitions from Elgar, 1985, JFM ... gives same result as Wunsch's def.
% using pressure converted or regular velocity should give the same coh & phase
% choose to use pressure converted in order to calc dirs later
% Cospectrum & Quadrature:

coPUpres = real(PUpres);   quPUpres = imag(PUpres);
coPVpres = real(PVpres);   quPVpres = imag(PVpres);
coUVpres = real(UVpres);   quUVpres = imag(UVpres);

% Coherence & Phase at each freq-band
% *** note that it's important to calc this AFTER all merging and ensemble avg.

warning off MATLAB:divideByZero
cohPUpres = sqrt( (coPUpres.^2 + quPUpres.^2) ./ (Spp.* UUpres) );
phPUpres  = (180/pi) .* atan2( quPUpres , coPUpres );  
cohPVpres = sqrt((coPVpres.^2 + quPVpres.^2)./ (Spp.* VVpres));
phPVpres  = (180/pi) .* atan2( quPVpres , coPVpres );  
cohUVpres = sqrt((coUVpres.^2 + quUVpres.^2)./(UUpres .* VVpres));
phUVpres  = (180/pi) .* atan2( quUVpres , coUVpres );  

%  WAVE DIRECTION & SPREAD ------------------------------------------------
%  before rotation 0 deg is for waves headed towards positive x (usually EAST)
%  method 1 (puv) covers 360 deg radius but can be prone to errors
%  method 2 (uv) covers only +/- 90 deg from x axis, possibly more accurate
%  see Herbers et al 1999 JGR

a1 = coPUpres ./ sqrt( Spp .* ( UUpres + VVpres ) );
b1 = coPVpres ./ sqrt( Spp .* ( UUpres + VVpres ) );
dir1 = radtodeg ( atan2(b1,a1) );          
spread1 = radtodeg ( sqrt( 2 .* ( 1-( a1.*cosd(dir1) + b1.*sind(dir1)) ) ) );
% spread1 = radtodeg ( sqrt( 2 .* ( 1-sqrt(a1.^2 + b1.^2) ) ) );

a2 = (UUpres - VVpres) ./ (UUpres + VVpres);
b2 = 2 .* coUVpres ./ ( UUpres + VVpres );
dir2 = radtodeg ( atan2(b2,a2)/2 ); 
spread2 = (180/pi) .* sqrt(abs( 0.5 - 0.5 .* ( a2.*cos(2.*degtorad(dir2)) +...
    b2.*sin(2.*degtorad(dir2)) )  ));

% ESTIMATE ENERGY FLUX -------------------------------------------------------
% spectral estimator from T.H.C. Herbers, using group velocity:

C = omega./k;
Cg = get_cg(k,depth);   %0.5 * omega.* k .* ( 1 + (2*k*depth)./sinh(2*k*depth) ); 

const = g * Cg .* ((cosh(k*depth)).^2) ./ ((cosh(k*doffp)).^2);  

om = omega;
CO = const;

% energy flux (by freq) in cartesian coordinates:

posX=const.*(0.5.* ( abs(Spp) + (UUpres - VVpres) ) + real(PUpres)); 
negX=const.*(0.5.* ( abs(Spp) + (UUpres - VVpres) ) - real(PUpres)); 
posY=const.*(0.5.* ( abs(Spp) + (VVpres - UUpres) ) + real(PVpres)); 
negY=const.*(0.5.* ( abs(Spp) + (VVpres - UUpres) ) - real(PVpres)); 

posX2= g * Cg .* a1 .* SSE;
posY2 = g * Cg .* b1 .* SSE;

%posY=const.*(0.5.* ( abs(Spp) + ( ( omega./(g*k) ).^2 ).*( abs(Svv) - abs(Suu) ) ) + (omega./(g*k)).*abs(real(Spv)));
%negY=const.*(0.5.* ( abs(Spp) + ( ( omega./(g*k) ).^2 ).*( abs(Svv) - abs(Suu) ) ) - (omega./(g*k)).*abs(real(Spv)));
%posX=const.*(0.5.* ( abs(Spp) + ( ( omega./(g*k) ).^2 ).*( abs(Suu) - abs(Svv) ) ) + (omega./(g*k)).*abs(real(Spu)));
%negX=const.*(0.5.* ( abs(Spp) + ( ( omega./(g*k) ).^2 ).*( abs(Suu) - abs(Svv) ) ) - (omega./(g*k)).*abs(real(Spu)));

% *** POSSIBLE PROBLEM WITH NEGATIVE ENERGY FLUXES ***
% use 2 lines below to set negative energy fluxes to be very small positives:
%posX( find(posX<10^-7) ) = 10^-7; negX( find(negX<10^-7) ) = 10^-7; 
%posY( find(posY<10^-7) ) = 10^-7; negY( find(negY<10^-7) ) = 10^-7; 
% create energy flux variable (for output)

Eflux = [posX2; posX; negX; posY2; posY; negY];

Eflux_ss = nansum(Eflux(:,i_swell),2) * df;
Eflux_ig = nansum(Eflux(:,i_ig),2) * df;


% Significant wave height (m)---------------------------------------------
Hsigt =  4 * sqrt( SSE(ii) * df ) ;   
Hsig_ss = 4 * sqrt( nansum( SSE(i_swell) * df ) );      
Hsig_ig = 4 * sqrt( nansum( SSE(i_ig) * df ) );      

Hrmst =  sqrt(8* SSE(ii) * df ) ;   
Hrms_ss = sqrt(8* nansum( SSE(i_swell) * df ) );      
Hrms_ig = sqrt(8* nansum( SSE(i_ig) * df ) );      
Hrms_all = sqrt(8* nansum( SSE(i_all) * df ) );   

%Hsig_ig = 4 * sqrt( nansum( SSE(i_ig) * df ) );            
%Hsig = Hsig_swell;

%HsigU_swell = 4 * sqrt( nansum( SSEU(i_swell) * df ) );      
%HsigU_ig = 4 * sqrt( nansum( SSEU(i_ig) * df ) );            
%HsigU = HsigU_swell;
% swell direction:
 
%dir1_swell = nansum( dir1(i_swell).*SSEU(i_swell) ) / nansum(SSEU(i_swell)) ;
%spread1_swell = nansum( spread1(i_swell).* SSEU(i_swell) ) / nansum(SSEU(i_swell)) ;
%dir2_swell = nansum( dir2(i_swell).* SSEU(i_swell) ) / nansum(SSEU(i_swell)) ;
%spread2_swell = nansum( spread2(i_swell).* SSEU(i_swell) ) / nansum(SSEU(i_swell)) ;

dirt  = [dir1(ii); dir2(ii)]  ;
spreadt =  [spread1(ii); spread2(ii)]  ;

dir_calc = dirt(dirmethod,:);

a1t = a1(ii);
a2t = a2(ii);
b1t = b1(ii);
b2t = b2(ii);

a1_ss = nansum(  a1(i_swell).*SSE(i_swell) ) / nansum(SSE(i_swell)) ;
b1_ss = nansum(  b1(i_swell).*SSE(i_swell) ) / nansum(SSE(i_swell)) ;
a2_ss = nansum(  a2(i_swell).*SSE(i_swell) ) / nansum(SSE(i_swell)) ;
b2_ss = nansum(  b2(i_swell).*SSE(i_swell) ) / nansum(SSE(i_swell)) ;

dir_ss1 = radtodeg ( atan2(b1_ss,a1_ss) );          
spread_ss1 = radtodeg ( sqrt( 2 .* ( 1-sqrt(a1_ss.^2 + b1_ss.^2) ) ) );

dir_ss2 = radtodeg ( atan2(b2_ss,a2_ss)/2 );          
spread_ss2 = (180/pi) .* sqrt(abs( 0.5 - 0.5 .* ( a2_ss.*cos(2.*degtorad(dir_ss2)) + b2_ss.*sin(2.*degtorad(dir_ss2)) )  ));

% centriod frequency
fcentroid_swell = nansum ( fm(i_swell).* SSE(i_swell) ) / sum ( SSE(i_swell) ) ;
Tm_ss = 1/fcentroid_swell;

% peak frequency
[~,indx] = max(SSE(i_swell));
temp = fmt(i_swell);
Tp_ss = 1./temp(indx);
if isempty(Tp_ss)
    Tp_ss=nan;
end

% ig dir and spread
a1_ig = nansum(  a1(i_ig).*SSE(i_ig) ) / nansum(SSE(i_ig)) ;
b1_ig = nansum(  b1(i_ig).*SSE(i_ig) ) / nansum(SSE(i_ig)) ;
a2_ig = nansum(  a2(i_ig).*SSE(i_ig) ) / nansum(SSE(i_ig)) ;
b2_ig = nansum(  b2(i_ig).*SSE(i_ig) ) / nansum(SSE(i_ig)) ;

dir_ig1 = radtodeg ( atan2(b1_ig,a1_ig) );          
spread_ig1 = radtodeg ( sqrt( 2 .* ( 1-sqrt(a1_ig.^2 + b1_ig.^2) ) ) );

dir_ig2 = radtodeg ( atan2(b2_ig,a2_ig)/2 );          
spread_ig2 = (180/pi) .* sqrt(abs( 0.5 - 0.5 .* ( a2_ig.*cos(2.*degtorad(dir_ig2)) + b2_ss.*sin(2.*degtorad(dir_ss2)) )  ));

% ig centriod frequency
fcentroid_ig = nansum ( fm(i_ig).* SSE(i_ig) ) / sum ( SSE(i_ig) ) ;
Tm_ig = 1/fcentroid_ig;

% ig peak frequency
[~,indx] = max(SSE(i_ig));
temp = fmt(i_ig);
Tp_ig = 1./temp(indx);

% centriod frequency
fcentroid_all = nansum ( fm(i_all).* SSE(i_all) ) / sum ( SSE(i_all) ) ;
Tm_all = 1/fcentroid_all;

% all peak frequency
[~,indx] = max(SSE(i_all));
temp = fmt(i_all);
Tp_all = 1./temp(indx);

% Radiation Stress Estimates ----------------------------------------------
% original formulation from Falk F.
Sxx = rho*g*( (1.5 + 0.5*a2) .* (Cg./C) - 0.5 ) .* SSE;
Syy = rho*g*( (1.5 - 0.5*a2) .* (Cg./C) - 0.5 ) .* SSE;
Sxy = rho*g*0.5*b2 .* (Cg./C) .* SSE;

Sxx_ss = nansum(Sxx(i_swell)) * df;
Syy_ss = nansum(Syy(i_swell)) * df;
Sxy_ss = nansum(Sxy(i_swell)) * df;

% this definition seems to be the same (Dean and Dalrymple)
% Sxx2 = rho*g*SSE.*( (Cg./C).*(cosd(dir_calc).^2+1) - 0.5 );
% Syy2 = rho*g*SSE.*( (Cg./C).*(sind(dir_calc).^2+1) - 0.5 );
% Sxy2 = 0.5*rho*g*SSE .* (Cg./C).* sind(2*dir_calc);
% 
% Sxx_ss = nansum(Sxx2(i_swell)) * df;
% Syy_ss = nansum(Syy2(i_swell)) * df;
% Sxy_ss = nansum(Sxy2(i_swell)) * df;

Cpu_ss = nansum( cohPUpres(i_swell).*SSE(i_swell) ) / nansum(SSE(i_swell));


% STOKES DRIFT -------------------------------------------------------------
% Us(f) = g*k(f)*H(f)^2/(8*om*h)*cos(theta(f));
% Vs(f) = g*k(f)*H(f)^2/(8*om*h)*cos(sin(f));
kt = k(ii);
omegat = omega(ii);

Ust = g*kt.*Hrmst.^2./(8*omegat.*depth).*cosd(dir_calc);
Vst = g*kt.*Hrmst.^2./(8*omegat.*depth).*sind(dir_calc);

Us = nansum(Ust) ;
Vs = nansum(Vst) ;

% DEGREES OF FREEDOM and level of no significant coherence --------------------
% DOF = 2 * (# independent windows) * (# bands merged)
merge = 1;
DOF = 2 * nA * merge;  
% 95% significance level for zero coherence
SIG = sqrt(6/DOF);
%phPUsig = phPUpres( find(cohPUpres > SIG) );
%freqPUsig = freq( find(cohPUpres > SIG) );
%phPVsig = phPVpres( find(cohPVpres > SIG) ); 
%freqPVsig = freq( find(cohPVpres > SIG) );

%%%% SPECIFY OUTPUT TO WaveStats structure file %%%%%%%%%%%%%%%%%%%%%%%%%%%
WaveStats.SSEt=SSEt;
WaveStats.Suut=Suut;
WaveStats.Svvt=Svvt;
WaveStats.Suvt=Suvt;
WaveStats.fmt=fmt;
WaveStats.Tt = 1./fmt; WaveStats.Tt(1)=nan;
WaveStats.dirt=dirt;
WaveStats.spreadt=spreadt;
WaveStats.Hsigt=Hsigt;
WaveStats.Hsig_ss=Hsig_ss;
WaveStats.Hsig_ig=Hsig_ig;
WaveStats.Hrmst=Hrmst;
WaveStats.Hrms_ss=Hrms_ss;
WaveStats.Hrms_ig=Hrms_ig;
WaveStats.Hrms_all=Hrms_all;
WaveStats.Tm_ss=Tm_ss;
WaveStats.Tm_ig=Tm_ig;
WaveStats.Tm_all=Tm_all;
WaveStats.Tp_ss=Tp_ss;
WaveStats.Tp_ig=Tp_ig;
WaveStats.Tp_all=Tp_all;
WaveStats.Eflux_ss=Eflux_ss;
WaveStats.Eflux_ig = Eflux_ig;
WaveStats.dir_ss1=dir_ss1;
WaveStats.dir_ss2=dir_ss2;
WaveStats.dir_ig1=dir_ig1;
WaveStats.dir_ig2=dir_ig2;
WaveStats.spread_ss1=spread_ss1;
WaveStats.spread_ss2=spread_ss2;
WaveStats.spread_ig1=spread_ig1;
WaveStats.spread_ig2=spread_ig2;
WaveStats.Sxx_ss=Sxx_ss;
WaveStats.Sxy_ss=Sxy_ss;
WaveStats.Syy_ss=Syy_ss;
WaveStats.a1t=a1t;
WaveStats.a2t=a2t;
WaveStats.b1t=b1t;
WaveStats.b2t=b2t;
WaveStats.ztest_ss=ztest_ss;
WaveStats.Cpu_ss=Cpu_ss;
WaveStats.depth=depth;
WaveStats.Us=Us;
WaveStats.Vs=Vs;
WaveStats.Ust=Ust;
WaveStats.Vst=Vst;
WaveStats.df=df;
WaveStats.num_avg = num_avg;


end
