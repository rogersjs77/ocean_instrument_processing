% function [ ] = aquadopp_process_data(bottom_trim,depth_max,filter_time,...
%     P_atm,Waves,WavesBurst,nfft,fcutoff,dirmethod,theta,head_height,extrap,z0,d,...
%     time_start,time_end,time_bad_vel,bins_from_surf)
%
% this file takes ADP data,
% file trims time sequence, depth, rotates to major/minor axis,
%
% Justin Rogers, stanford efml
% June 2017
%
% input.bottom_trim = data to drop from bottom (m), if empty [] uses average
%               =0 as criteria
% depth_max = data to trim from top (m), if empty [], uses Rmax constraint
%               based on angle and depth
% filter_time = low pass filter window (min) for u and waves;
% Waves flag,   if 1= process wave module, 0=don't process
% nfft  = number of points to use for waves fft, affects freq resolution 
%       =5*60/dt_sec;% approximate 5 minute window for ig waves
%       =1*60/dt_sec, appoximate 1 minute window for swell waves
% dirmethod = 1, use a1, b1 (puv) correlations for direction, can cover full 360 deg
%           but possible errors from reflections, less accurate?
%           = 2, use a2, b2 (uv) correlations for direction, can only cover
%           +/- 90 deg from x axis. Only use this method if you know waves
%           are in + x direction, but this is likely more accurate. See
%           Herbers et al, 1999 JGR for more discussion.
% theta = angle in degrees to rotate coordinates, if [] finds maximum U variance
% P_atm = atmospheric pressure in Pa - remember recorded pressure in ADCP
% is in dPa = 10*Pa;
% head_height is height of ADP head above bottom (m), for adp profiler h=0.037m.
% extrap command = 1 extrapolates missing u data in water column using log
% layer,
%
% frequency resolution is: df = fs/(nfft-1); fs is sampling frequency
% number of averages num = floor(4*n/nfft)-3; n is # of points in dt_big
%%

function [ ] = aquadopp_process_data(bottom_trim,depth_max,filter_time,...
    P_atm,Waves,WavesBurst,nfft,fcutoff,dirmethod,theta,head_height,extrap,z0,d,...
    time_start,time_end,time_bad_vel,bins_from_surf)

files=dir('*.hdr');   

for i=1:length(files)

files=dir('*.hdr');
    
[~,filename,~]=fileparts(files(i).name);

if nargin<17
    bins_from_surf=2;
end

%% load data
adp.hdr = importdata(strcat(filename,'.hdr'),'');
adp.cells = cell2mat(textscan(char(adp.hdr(10)),'%*s %*s %*s %f %*s'));
adp.cell_size = cell2mat(textscan(char(adp.hdr(11)),'%*s %*s %f %*s'))/100;
adp.blanking_dist= cell2mat(textscan(char(adp.hdr(15)),'%*s %*s %f %*s'));
adp.sample_time = cell2mat(textscan(char(adp.hdr(9)),'%*s %*s %f %*s'));

adp.vel_east = dlmread(strcat(filename,'.v1'))';
adp.vel_north = dlmread(strcat(filename,'.v2'))';
adp.vel_vert = dlmread(strcat(filename,'.v3'))';
adp.sen =  dlmread(strcat(filename,'.sen'));

adp.amp1 =  dlmread(strcat(filename,'.a1'))';
adp.amp2 =  dlmread(strcat(filename,'.a2'))';
adp.amp3 =  dlmread(strcat(filename,'.a3'))';

 

%% trim mtime to time underwater
adp.mtime = datenum(adp.sen(:,3),adp.sen(:,1),adp.sen(:,2),...
    adp.sen(:,4),adp.sen(:,5),adp.sen(:,6))';

adp.pressure = adp.sen(:,14)';

if isempty(P_atm) % if no manually defined atmospheric P
    % find avg atmospheric pressure within expected range < 1E4
    P_atm = 1E4*nanmedian(adp.pressure(adp.pressure<1 & adp.pressure>0));
    if isempty(P_atm) || isnan(P_atm)
      P_atm = 0;
    end
end

adp.press_Pa = 1E4*adp.pressure-P_atm;
adp.press_Pa(adp.press_Pa>2E6)=nan;% remove unrealistic pressures

% find pressures outside 1 std of middle 1/2 of data set
if isempty(time_start) % if no manual trim defined
indx3 = find(adp.press_Pa>(median(adp.press_Pa(round(0.25*end:0.75*end)))...
    -std(adp.press_Pa(round(0.25*end:0.75*end)))));
time_start=adp.mtime(indx3(1));
time_end=adp.mtime(indx3(end));
end

% adp time is not always constant spacing
mtime=time_start:adp.sample_time/3600/24:time_end;
dt = mtime(2)-mtime(1); % data time step (days)

% always use nearest interp method for raw data
press_Pa = interp1(adp.mtime,adp.press_Pa,mtime,'nearest',NaN);

% find density and depth
S=34.6; g=9.81;
temperature = interp1(adp.mtime,adp.sen(:,15)',mtime,'nearest',NaN);
[~,sigma] = swstate(S,nanmean(temperature),nanmean(press_Pa)*1E-4);
rho = 1000+sigma;
depth = press_Pa./(rho*g)+head_height; % hydrostatic depth above bottom in m
Depth = nanmean(depth);

heading = interp1(adp.mtime,adp.sen(:,11)',mtime,'linear',NaN);
pitch = interp1(adp.mtime,adp.sen(:,12)',mtime,'linear',NaN);
roll = interp1(adp.mtime,adp.sen(:,13)',mtime,'linear',NaN);
voltage = interp1(adp.mtime,adp.sen(:,9)',mtime,'linear',NaN);

% low pass filter
dt_filt = filter_time/60/24/2;
mtime_filt = time_start:dt_filt:time_end;

[pressure_filt,Fourier.pressure,~] = low_pass_filter(press_Pa,dt*24*3600,filter_time*60);
[depth_filt,Fourier.depth,Fourier.s_depth] = low_pass_filter(depth,dt*24*3600,filter_time*60);
[temperature_filt,Fourier.temperature,~] = low_pass_filter(temperature,dt*24*3600,filter_time*60);
[depth_tides,~,~] = low_pass_filter(depth,dt*24*3600,2*3600);

depth_filt = interp1(mtime,depth_filt,mtime_filt,'linear');
pressure_filt = interp1(mtime,pressure_filt,mtime_filt,'linear');
temperature_filt = interp1(mtime,temperature_filt,mtime_filt,'linear');
depth_tides = interp1(mtime,depth_tides,mtime_filt,'linear');

%% trim vertical bins
bin_MAB_orig = adp.blanking_dist+adp.cell_size+adp.cell_size*[0:1:adp.cells-1];
N=length(depth);
% calculate vertical data to drop
if ~isempty(depth_max) % specified constraint
    Rmax = depth_max;
else % automated calculation
%     Rmax = depth_filt*cos(deg2rad(25));
    Rmax = min(depth(round(N/4):round(N*3/4)))*cos(deg2rad(25));
end
% bins_drop_top = length(find(nanmax(Rmax)<bin_MAB_orig));
bins_drop_top = length(find(Rmax<bin_MAB_orig));

if isempty(bottom_trim) % base bottom trim on min adcp.corr and vel_error
bins_drop_bottom=0;
else
bins_drop_bottom = length(find(bottom_trim>bin_MAB_orig));
end

bin_MAB = bin_MAB_orig(1+bins_drop_bottom:end-bins_drop_top);
bin_z = bin_MAB-Depth;
bins_good=length(bin_MAB);

% drop vertical bad bins
vel_east2 = adp.vel_east(1+bins_drop_bottom:end-bins_drop_top,:);
vel_north2 = adp.vel_north(1+bins_drop_bottom:end-bins_drop_top,:);
vel_vert2 = adp.vel_vert(1+bins_drop_bottom:end-bins_drop_top,:);

% remove manually specified bad velocity data
for kk=1:length(time_bad_vel)/2
    indxv = time_bad_vel(kk) < adp.mtime & time_bad_vel(kk+1) > adp.mtime;
    for j=1:size(vel_east2,1)
       vel_east2(j,indxv)=nan;
       vel_north2(j,indxv)=nan;
       vel_vert2(j,indxv)=nan;
    end
end
% interpolate to common mtime & nan velocity below Rmax
% Rmax_full = depth*cos(deg2rad(25));
% always use nearest interp method on raw data
for j=1:size(vel_east2,1)
   vel_east(j,:) = interp1(adp.mtime,vel_east2(j,:),mtime,'nearest');
   vel_north(j,:) = interp1(adp.mtime,vel_north2(j,:),mtime,'nearest');
   vel_vert(j,:) = interp1(adp.mtime,vel_vert2(j,:),mtime,'nearest');
   
%    indx = Rmax_full < bin_MAB(j);
%    vel_east(j,indx) = NaN;
%    vel_north(j,indx) = NaN;
%    vel_vert(j,indx) = NaN;
end
clear vel_east2 vel_north2 vel_vert2

% calculate error, 1% of value +/- 0.5 cm/s
vel_error = 0.01*sqrt(vel_east.^2+vel_north.^2+vel_vert.^2)+0.5/100;




%% filter data

% [vel_east_filt,Fourier.vel_east,~]=low_pass_filter(naninterp(vel_east),dt*24*3600,filter_time*60);
% [vel_north_filt,Fourier.vel_north,~]=low_pass_filter(naninterp(vel_north),dt*24*3600,filter_time*60);
% [vel_vert_filt,Fourier.vel_vert,~]=low_pass_filter(naninterp(vel_vert),dt*24*3600,filter_time*60);
% 
% vel_east_filt2 = interp1(mtime,vel_east_filt',mtime_filt,'linear',nan)';
% vel_north_filt2 = interp1(mtime,vel_north_filt',mtime_filt,'linear',nan)';
% vel_vert_filt2 = interp1(mtime,vel_vert_filt',mtime_filt,'linear',nan)';

 dt_big_day = 1/24/60*(filter_time);   
 for i=1:length(mtime_filt)
    tdown = mtime_filt(i)-dt_big_day/2;
    tup = mtime_filt(i)+dt_big_day/2;
    indx = find(mtime>tdown & mtime<=tup);
        
   % find rms values of all velocity
    u = vel_east(:,indx);
    v = vel_north(:,indx);
    w = vel_vert(:,indx);
    uerr = vel_error(:,indx);
    
    % remove unrealistic values
    uerr(abs(u)>2)=NaN;
    u(abs(u)>2)=NaN;
    v(abs(v)>2)=NaN;
    w(abs(w)>0.5) = NaN;
           
%     indxnan = sum(~isnan(u),2) < 20; % sample size bigger than 20 
    
    vel_east_filt(:,i) = nanmean(u,2);
    vel_north_filt(:,i) = nanmean(v,2);
    vel_vert_filt(:,i) = nanmean(w,2);
    
    % estimate error from measurement, std_b=(1/n)*sqrt(sum(std^2))
    vel_err_filt(:,i) = sqrt(nansum(uerr.^2,2))./sum(~isnan(u),2);
    
    % standard error 95% from variance err=std/sqrt(n)    
%     vel_err_east(:,i) = nanstd(u,0,2)./sqrt(sum(~isnan(u),2));
%     vel_err_north(:,i) = nanstd(v,0,2)./sqrt(sum(~isnan(v),2));
%     vel_err_vert(:,i) = nanstd(w,0,2)./sqrt(sum(~isnan(w),2));
%  
%     vel_east_filt(indxnan,i) = NaN;
%     vel_north_filt(indxnan,i) = NaN;
%     vel_vert_filt(indxnan,i) = NaN;    
%     vel_err_filt(indxnan,i) = NaN;
%     vel_err_east(indxnan,i) = NaN;
%     vel_err_north(indxnan,i) = NaN;
%     vel_err_vert(indxnan,i) = NaN;
    
    clear u v w uerr indxnan indx
 end
 
 
 %% low pass filter 
 [~,Fourier.Vel_East,~] = low_pass_filter(naninterp(nanmean(vel_east,1)),dt*24*3600,filter_time*60);
[~,Fourier.Vel_North,~] = low_pass_filter(naninterp(nanmean(vel_north,1)),dt*24*3600,filter_time*60);

 
 %% nan data with measurement greater than depth

for i=1:length(mtime)
   indx = bin_MAB > (depth(i)-bins_from_surf*(bin_MAB(2)-bin_MAB(1)));
   vel_east(indx,i)=nan;
   vel_north(indx,i)=nan;
   vel_vert(indx,i)=nan;
   vel_error(indx,i)=nan;
end
dtemp = naninterp(depth_filt);
for i=1:length(mtime_filt)
   indx = bin_MAB > (dtemp(i)-bins_from_surf*(bin_MAB(2)-bin_MAB(1)));
   vel_east_filt(indx,i)=nan;
   vel_north_filt(indx,i)=nan;
   vel_vert_filt(indx,i)=nan;
   vel_err_filt(indx,i)=nan;
end
 clear dtemp
%% depth average mean velocities
Vel_East = nanmean(vel_east_filt,1); % depth averaged values
Vel_North=nanmean(vel_north_filt,1);
Vel_Vert=nanmean(vel_vert_filt,1);
Vel_Err = sqrt(nansum(vel_err_filt.^2,1))./size(vel_err_filt,1);


%% compute angle of maximum variance for mean currents

% depth avg current direction of max variance
[theta_meanU]=adcp_paa(Vel_East,Vel_North);

% instantaneous u direction of max variance
bin = round(length(bin_MAB/4));
[theta_u]=adcp_paa(vel_east(bin,:),vel_north(bin,:));

if isempty(theta) % if no input theta is defined
theta=theta_meanU;
end

%% check compass drift
if mtime(end)-mtime(1)>28
mtime_compass = (mtime_filt(1)-14):28:(mtime_filt(end)+14);
for i=1:length(mtime_compass)
   indx = mtime_filt > mtime_compass(i)-14 &...
       mtime_filt < mtime_compass(i)+14;
   [temp(i)]=adcp_paa(Vel_East(indx),Vel_North(indx)); 
end
theta_drift = interp1(mtime_compass,temp,mtime,'linear','extrap');
theta_drift_filt = interp1(mtime_compass,temp,mtime_filt,'linear','extrap');
end
%% rotate coordinate frame

[vel_maj_filt, vel_min_filt]= adcp_rotation(vel_east_filt,vel_north_filt, theta);
[vel_maj, vel_min]= adcp_rotation(vel_east,vel_north, theta);   

% [~,Fourier.Vel_Maj,Fourier.s_vel] = low_pass_filter(naninterp(mean(vel_maj,1)),dt*24*3600,filter_time*60);
% [~,Fourier.Vel_Min,~] = low_pass_filter(naninterp(mean(vel_min,1)),dt*24*3600,filter_time*60);


%% Compute extrapolated profiles

if extrap

% compute overall roughness height z0
vel_min_corr = 0.020;
[roughness] = get_z0_log_layer(vel_maj_filt,vel_min_filt,bin_MAB,Depth,vel_min_corr);
 
%%% !! REVIEW z0 RESULTS BEFORE SELECTING !! %%%
if isempty(z0)
z0_u = roughness.z0.mean;
z0_v = roughness.z0.mean;
d_u = roughness.d.mean;
d_v = roughness.d.mean;
else
z0_u = z0;
z0_v = z0;
d_u = d;
d_v = d;
end

% interpolate profile
bin_MAB_extrap = linspace(0,max(depth_filt),4*round(Depth/adp.cell_size(1,1)));
bin_z_extrap = bin_MAB_extrap-max(depth_filt);

% interpolate u,v,w
[vel_maj_filt_extrap,vel_min_filt_extrap,vel_vert_filt_extrap,ustar,vstar] = ...
    velocity_profile_interp_log(vel_maj_filt,vel_min_filt,vel_vert_filt,z0_u,z0_v,...
    d_u,d_v,bin_MAB,bin_MAB_extrap,depth_filt);

[vel_east_filt_extrap,vel_north_filt_extrap,~,ustar,vstar] = ...
    velocity_profile_interp_log(vel_east_filt,vel_north_filt,vel_vert_filt,z0_u,z0_v,...
    d_u,d_v,bin_MAB,bin_MAB_extrap,depth_filt);

%compute extrapolated profiles
Vel_Maj_Extrap_Prof=nanmean(vel_maj_filt_extrap,2);
Vel_Min_Extrap_Prof=nanmean(vel_min_filt_extrap,2);
Vel_Vert_Extrap_Prof=nanmean(vel_vert_filt_extrap,2);

Vel_East_Extrap_Prof=nanmean(vel_east_filt_extrap,2);
Vel_North_Extrap_Prof=nanmean(vel_north_filt_extrap,2);

Vel_Maj = nanmean(vel_maj_filt_extrap,1);
Vel_Min = nanmean(vel_min_filt_extrap,1);
Vel_Vert = nanmean(vel_vert_filt_extrap,1);

Vel_East = nanmean(vel_east_filt_extrap,1);
Vel_North = nanmean(vel_north_filt_extrap,1);

else % if no extrapolation, use best estimate
Vel_Maj = nanmean(vel_maj_filt,1);
Vel_Min = nanmean(vel_min_filt,1);
Vel_East = nanmean(vel_east_filt,1);
Vel_North = nanmean(vel_north_filt,1);

end

qe_Maj = depth_filt.*Vel_Maj;
qe_Min = depth_filt.*Vel_Min;


%% Compute Avg vertical profiles
indx_pos = naninterp(Vel_Maj) > 0 ;
indx_neg = naninterp(Vel_Maj) < 0 ;
Vel_Maj_Prof(1,:)=sqrt(nanmean(vel_maj_filt(:,indx_pos).^2,2));
Vel_Maj_Prof(2,:)=sqrt(nanmean(vel_maj_filt(:,indx_neg).^2,2));

indx_pos = naninterp(Vel_Min) > 0 ;
indx_neg = naninterp(Vel_Min) < 0 ;
Vel_Min_Prof(1,:)=sqrt(nanmean(vel_min_filt(:,indx_pos).^2,2));
Vel_Min_Prof(2,:)=sqrt(nanmean(vel_min_filt(:,indx_neg).^2,2));

%% tidal analysis
eta = depth-nanmean(depth);
eta_filt = depth_filt-nanmean(depth_filt);
if mtime(end)-mtime(1)>2
[TIDES_d.nameu,TIDES_d.fu,TIDES_d.tidecon,eta_tides]=...
    t_tide(eta_filt,'interval',24*dt_filt); 
eta_notides = eta_filt-eta_tides;   
eta_2hrfilt = depth_tides-nanmean(depth_tides);

%tidal velocity anaylsis
[TIDES_UMaj.nameu,TIDES_UMaj.fu,TIDES_UMaj.tidecon,Vel_Maj_tides]=...
    t_tide(Vel_Maj,'interval',24*dt_filt); 
Vel_Maj_notides = Vel_Maj-Vel_Maj_tides;   

[TIDES_UMin.nameu,TIDES_UMin.fu,TIDES_UMin.tidecon,Vel_Min_tides]=...
    t_tide(Vel_Min,'interval',24*dt_filt); 
Vel_Min_notides = Vel_Min-Vel_Min_tides;
end
%% Process Wave Data

if Waves
    
    if WavesBurst        
        adp.wave.interval = cell2mat(textscan(char(adp.hdr(20)),'%*s %*s %*s %f %*s'));
        adp.wave.samples = cell2mat(textscan(char(adp.hdr(21)),'%*s %*s %*s %*s %*s %f %*s'));
        adp.wave.sample_rate = cell2mat(textscan(char(adp.hdr(22)),'%*s %*s %*s %*s %f %*s'));
        fs = 1./adp.wave.sample_rate;
        
        % load whd general burst info
        DATA = dlmread(strcat(filename,'.whd'));
        % find center of sampling time (add half time of sample period to
        % start)
        adp.wave.mtime = datenum(DATA(:,3),DATA(:,1),DATA(:,2),...
            DATA(:,4),DATA(:,5),DATA(:,6))'+...
            adp.wave.samples/(2*adp.wave.sample_rate*3600*24);
        indx = adp.wave.mtime>time_start&adp.wave.mtime<time_end;
        adp.wave.mtime = adp.wave.mtime(indx);
        adp.wave.burst = DATA(indx,7)';
        doffu = nanmedian(DATA(indx,9));
        doffp = head_height;
        
        % load full wad burst data
        DATA = dlmread(strcat(filename,'.wad'));
        adp.wave.burst_wad = DATA(:,1)';
        Pw = DATA(:,3)'-P_atm*1E-4;
        Ww = DATA(:,8)'; 
        [Uw,Vw] = adcp_rotation(DATA(:,6)',DATA(:,7)',theta);
        adp.wave.Pw = Pw;
        adp.wave.Uw = Uw;
        adp.wave.Vw = Vw;
        clear DATA   

    else % continuous sampling
        bin = round(length(bin_MAB)/8); % use bin near bottom, %find(Error==min(Error)); % select bin with lowest total error
        % bin = floor(length(bin_MAB)/2);% select bin at midpoint
        dt_sec = dt*24*3600;
        doffu = bin_MAB(bin);
        doffp = head_height;
        % dirmethod = 1; % use method 1 if direction unknown, use 2 if waves in +x direction.

        % use rotated coordinates
        Uw = vel_maj(bin,:);
        Vw = vel_min(bin,:);
        Ww = vel_vert(bin,:);
        Pw = press_Pa;
        % use NE coordinates
        % Uw = vel_east(bin,:);
        % Vw = vel_north(bin,:);

        fs = 1./dt_sec;
        dt_big_day = 1/24/60*(filter_time);
        % only process points with full time series
        adp.wave.mtime = mtime_filt(2:end-1);
    end
    
 for i=1:(length(adp.wave.mtime)) 
    if WavesBurst
       indx = adp.wave.burst(i)==adp.wave.burst_wad;         
    else
        tdown = mtime_filt(i)-dt_big_day/2;
        tup = mtime_filt(i)+dt_big_day/2;
        indx = find(mtime>tdown & mtime<=tup);
    end
        
    if ~isempty(indx)
    U = Uw(indx);
    V = Vw(indx);
    W = Ww(indx);
    P = Pw(indx); % pressure is in dbar

    [WaveStats] = wave_stats_spectra(U,V,P,nfft,doffu,doffp,fs,fcutoff,rho,dirmethod);

    WaveData.depth(i) = WaveStats.depth;
    WaveData.Hsig_ss(i) = WaveStats.Hsig_ss;
    WaveData.Hsig_ig(i) = WaveStats.Hsig_ig;
    WaveData.HSIG(i,:) = WaveStats.Hsigt;
    
    WaveData.Hrms_ss(i) = WaveStats.Hrms_ss;
    WaveData.Hrms_ig(i) = WaveStats.Hrms_ig;
    WaveData.Hrms_all(i) = WaveStats.Hrms_all;
    WaveData.HRMS(i,:) = WaveStats.Hrmst;
    
    WaveData.Tm_ss(i) = WaveStats.Tm_ss;
    WaveData.Tm_ig(i) = WaveStats.Tm_ig;
    WaveData.Tm_all(i) = WaveStats.Tm_all;
    
    WaveData.Tp_ss(i) = WaveStats.Tp_ss;
    WaveData.Tp_ig(i) = WaveStats.Tp_ig;
    WaveData.Tp_all(i) = WaveStats.Tp_all;
        
    WaveData.Eflux_ss(i,:) = WaveStats.Eflux_ss';
    WaveData.ztest_ss(i) = WaveStats.ztest_ss;
    WaveData.Cpu_ss(i) = WaveStats.Cpu_ss;
    
    if dirmethod==1
    WaveData.dir_ss(i) = WaveStats.dir_ss1;
    WaveData.dir_ig(i) = WaveStats.dir_ig1;
    WaveData.DIR(i,:) = WaveStats.dirt(1,:);
    else
    WaveData.dir_ss(i) = WaveStats.dir_ss2;
    WaveData.dir_ig(i) = WaveStats.dir_ig2;
    WaveData.DIR(i,:) = WaveStats.dirt(2,:);
    end
    
    WaveData.spread_ss(i) = WaveStats.spread_ss2;
    WaveData.spread_ig(i) = WaveStats.spread_ig2;
    WaveData.SPREAD(i,:) = WaveStats.spreadt(2,:);
    
    WaveData.Eflux_ig(i,:) = WaveStats.Eflux_ig';
    
    WaveData.Sxx_ss(i) = WaveStats.Sxx_ss;
    WaveData.Sxy_ss(i) = WaveStats.Sxy_ss;
    WaveData.Syy_ss(i) = WaveStats.Syy_ss;
    
    WaveData.SUU(i,:) = WaveStats.Suut;
    WaveData.SUV(i,:) = WaveStats.Suvt;
    WaveData.SVV(i,:) = WaveStats.Svvt;
    WaveData.SSE(i,:) = WaveStats.SSEt;
    
    WaveData.USt(i,:) = WaveStats.Ust;
    WaveData.VSt(i,:) = WaveStats.Vst;
    WaveData.Us(i) = WaveStats.Us;
    WaveData.Vs(i) = WaveStats.Vs;
    
    % find rms values, separate waves/turbulence
    if WavesBurst
        u = Uw(indx);
        v = Vw(indx);
        w = Ww(indx);
        p = Pw(indx);
    else
        u = vel_maj(:,indx);
        v = vel_min(:,indx);
        w = vel_vert(:,indx);
        p = press_Pa(indx)*1E-4; % pressure is in dbar
    end
   
    for jj = 1:size(u,1)  % do for each vertical bin 
    up(jj,:) = u(jj,:) - nanmean(u(jj,:));
    vp(jj,:) = v(jj,:) - nanmean(v(jj,:));
    wp(jj,:) = w(jj,:) - nanmean(w(jj,:));
    
    % wave current separation
    [B]= benilov(u(jj,:),v(jj,:),w(jj,:),p,nfft,doffp,fs,fcutoff,rho);
    WaveData.S_uwave_wwave(jj,i,:) = B.S_uwave_wwave;
    WaveData.S_vwave_wwave(jj,i,:) = B.S_vwave_wwave;
    WaveData.Suw(jj,i,:) = B.Suw;
    WaveData.Svw(jj,i,:) = B.Svw;
    WaveData.Supwp(jj,i,:) = B.Supwp;
    WaveData.Svpwp(jj,i,:) = B.Svpwp;
    
    % Reynolds stress tensor
    WaveData.upup_bar(jj,i) = B.upup_bar;
    WaveData.vpvp_bar(jj,i) = B.vpvp_bar;
    WaveData.wpwp_bar(jj,i) = B.wpwp_bar;
    WaveData.upwp_bar(jj,i) = B.upwp_bar;
    WaveData.vpwp_bar(jj,i) = B.vpwp_bar;
    WaveData.upvp_bar(jj,i) = B.upvp_bar;
    WaveData.TKE(jj,i) = 0.5*(B.upup_bar+B.vpvp_bar+B.wpwp_bar);
    
    % Wave velocity tensor
    WaveData.uuwave_rms(jj,i) = sqrt(B.uwave_uwave_bar);
    WaveData.vvwave_rms(jj,i) = sqrt(B.vwave_vwave_bar);
    WaveData.wwwave_rms(jj,i) = sqrt(B.wwave_wwave_bar);
    WaveData.uwwave_rms(jj,i) = sqrt(B.uwave_wwave_bar);
    WaveData.vwwave_rms(jj,i) = sqrt(B.vwave_wwave_bar);
    WaveData.uvwave_rms(jj,i) = sqrt(B.uwave_vwave_bar);

    end
    
    WaveData.uvw_rms(:,i) = sqrt(nanmean(u.^2+v.^2+w.^2,2));
    WaveData.u_rms(:,i) = sqrt(nanmean(u.^2,2));
    WaveData.v_rms(:,i) = sqrt(nanmean(v.^2,2));
    WaveData.w_rms(:,i) = sqrt(nanmean(w.^2,2));
        
    WaveData.uvwp_rms(:,i) = sqrt(nanmean(up.^2+vp.^2+wp.^2,2));
    WaveData.up_rms(:,i) = sqrt(nanmean(up.^2,2));
    WaveData.vp_rms(:,i) = sqrt(nanmean(vp.^2,2));
    WaveData.wp_rms(:,i) = sqrt(nanmean(wp.^2,2));   
   
    % find stress parameters <u|u|>
    WaveData.uabsu_bar(:,i) = nanmean(u.*sqrt(u.^2+v.^2),2);
    WaveData.vabsu_bar(:,i) = nanmean(v.*sqrt(u.^2+v.^2),2);
    
    end
    clear U P V up vp wp u v w
 end
 WaveData.fmt = WaveStats.fmt;
 WaveData.mtime = adp.wave.mtime;
 WaveData.nfft = nfft;
 WaveData.doffp = doffp;
 WaveData.doffu = doffu;
 WaveData.fs = fs;
 WaveData.fcutoff = fcutoff;
 WaveData.df = WaveStats.df;
 WaveData.num_avg = WaveStats.num_avg;
 
 % bottom stress    
 WaveData.Tau_bx = rho*WaveData.upwp_bar;
 WaveData.Tau_by = rho*WaveData.vpwp_bar;
 WaveData.Tau = sqrt(WaveData.Tau_bx.^2+WaveData.Tau_by.^2);
 
 % compute wave heading direction from angle (CCW from x)
 WaveData.dir_heading_ss = -WaveData.dir_ss+90-theta;
 WaveData.dir_heading_ss(WaveData.dir_heading_ss<0) = ...
      WaveData.dir_heading_ss(WaveData.dir_heading_ss<0)+360;
 WaveData.dir_heading_ig = -WaveData.dir_ig+90-theta;
 WaveData.dir_heading_ig(WaveData.dir_heading_ig<0) = ...
      WaveData.dir_heading_ig(WaveData.dir_heading_ig<0)+360;
 WaveData.DIR_heading = -WaveData.DIR+90-theta;
 WaveData.DIR_heading(WaveData.DIR_heading<0) = ...
      WaveData.DIR_heading(WaveData.DIR_heading<0)+360;
   
  clear WaveStats doffp doffu fs bin Error 

 % Lagrangian mean flows in xy direction
 Vel_Maj_Stokes = interp1(WaveData.mtime,WaveData.Us,mtime_filt,'linear',NaN);
 Vel_Min_Stokes = interp1(WaveData.mtime,WaveData.Vs,mtime_filt,'linear',NaN);

 Vel_Maj_Lagr = Vel_Maj+Vel_Maj_Stokes;
 Vel_Min_Lagr = Vel_Min+Vel_Min_Stokes;
% rotate these back to NE coordinates
 [Vel_East_Stokes,Vel_North_Stokes]= adcp_rotation(Vel_Maj_Stokes,Vel_Min_Stokes, -theta);
 [Vel_East_Lagr,Vel_North_Lagr]= adcp_rotation(Vel_Maj_Lagr,Vel_Min_Lagr, -theta);

end

%% clean up data
input.bottom_trim=bottom_trim;
input.depth_max=depth_max;
input.filter_time=filter_time;
input.Waves=Waves;
input.nfft=nfft;

Misc.bins_drop_top = bins_drop_top;
Misc.bins_good = bins_good;
Misc.bins_drop_bottom = bins_drop_bottom;

%% save files
clear i j ii bottom_trim depth_max bad_beam Waves nfft corr_min...
    N indx indx2 indx3 indx_corr jj kk...
    time_start time_end time_id bins_drop_top bins_good ...
    bins_drop_bottom beta 



save('ADP_PostProcess.mat','-v7.3');
end



end