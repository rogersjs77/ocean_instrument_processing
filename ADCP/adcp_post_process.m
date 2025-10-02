% function [ ] = adcp_post_process(bottom_trim,depth_max,filter_time,...
%     bad_beam,P_atm,Waves,nfft,fcutoff,dirmethod,corr_min,theta,head_height,...
%     extrap,z0,d,time_start,time_end,sample_size,time_bad_vel)
%
% this file takes ADCP data from rdradcp process,
% file trims time sequence, depth, rotates to major/minor axis,
%
% Justin Rogers, stanford efml,
% June 2017
%
% input.bottom_trim = data to drop from bottom (m), if empty [] uses average
%               adcp.corr<115 as criteria
% depth_max = data to trim from top (m), if empty [], uses Rmax constraint
%               based on angle and depth
% filter_time = low pass filter window (min) for u and waves;
% bad_beam = adcp beam to remove, leave blank [] if keeping fulll solution
% Waves flag,   if 1= process wave module, 0=don't process
% nfft  = number of points to use for waves fft, affects freq resolution 
%       =5*60/dt_sec;% approximate 5 minute window for ig waves
%       =1*60/dt_sec, appoximate 1 minute window for swell waves
% corr_min = minimum corr value, if empty [] defaults to 0 (no filtering)
% theta = angle in degrees to rotate coordinates, if [] finds maximum U variance
% P_atm = atmospheric pressure in Pa - remember recorded pressure in ADCP
% is in dPa = 10*Pa;
% head_height is height of ADCP head above bottom (m), for WH Monitor = 0.202
% m, for WH Sentinel = 0.396 m.
% extrap command = 1 extrapolates missing u data in water column using log
% layer,
% time_start = optional starting time in matlab format
% time_end = optional end time in matlab format
%
% frequency resolution is: df = fs/(nfft-1); fs is sampling frequency
% number of averages num = floor(4*n/nfft)-3; n is # of points in dt_big
%%

function [ ] = adcp_post_process(bottom_trim,depth_max,filter_time,...
    bad_beam,P_atm,Waves,nfft,fcutoff,dirmethod,corr_min,theta,head_height,...
    extrap,z0,d,time_start,time_end,sample_size,time_bad_vel)


%% Combine data from multiple files
filesrdr=dir('rdrADCP*.mat');
for ii=1:size(filesrdr,1);
filenamerdr=filesrdr(ii).name;  
load(filenamerdr)

%combine rdrADCP files into one
if ii==1
    ADCP=adcp;
    CFG=cfg;
    ENS=ens;
    HDR=hdr;
    fadcp=fieldnames(adcp);
    fcfg=fieldnames(cfg);
    fhdr=fieldnames(hdr);
else % process files > 1
   
    for jj = 1:length(fadcp)
        n=ndims(ADCP.(fadcp{jj}));
        ADCP.(fadcp{jj})=cat(n,ADCP.(fadcp{jj}), adcp.(fadcp{jj}));    
    end
    
    for jj = 1:length(fcfg)
        CFG.(fcfg{jj})=[CFG.(fcfg{jj}) cfg.(fcfg{jj})];
    end
    
    for jj = 1:length(fhdr)
        HDR.(fhdr{jj})=[HDR.(fhdr{jj}) hdr.(fhdr{jj})];
    end    
    
end

end
adcp=ADCP;
cfg=CFG;
ens=ENS;
hdr=HDR;
clear jj ii fadcp fcfg fens fhdr ADCP CFG ENS HDR temp n

% save('./ADCP_Combined_RawData.mat','-v7.3','adcp','cfg','ens','hdr')

%% Post process ADCP data

%%  clean up adcp structure file
indx= adcp.mtime<datenum(1990,1,1) | ~isfinite(adcp.mtime); % clean up zeros and nans from mtime
adcp.mtime(indx)=nan;
adcp.mtime(isnan(adcp.mtime)) = interp1(find(~isnan(adcp.mtime)), ...
    adcp.mtime(~isnan(adcp.mtime)), find(isnan(adcp.mtime)),'linear','extrap');

% set important variables in adcp to nan where mtime was bad
adcp.pressure(indx) = nan;
adcp.pressure(adcp.pressure>1E7)=nan; % remove unrealistic pressures
adcp.temperature(indx) = nan;
adcp.depth(indx) = nan;
adcp.heading(indx) = nan;
adcp.pitch(indx) = nan;
adcp.roll(indx) = nan;
for i=1:size(adcp.east_vel,1)
   adcp.east_vel(i,indx) = nan;
   adcp.north_vel(i,indx) = nan;
   adcp.vert_vel(i,indx) = nan;
   adcp.error_vel(i,indx) = nan;
end


if isempty(P_atm) % if no manually defined atmospheric P
    % find avg atmospheric pressure within expected range < 1E3    
    P_atm = 10*nanmedian(adcp.pressure(adcp.pressure<1E3 & adcp.pressure>0));
    if isempty(P_atm) || isnan(P_atm)
      P_atm = 0;
    end
end
% find normalized pressure
adcp.press_Pa = 10*adcp.pressure-P_atm;
adcp.press_Pa(adcp.press_Pa>2E6)=nan; % remove unrealistic pressures

% find pressures outside 1 std of middle 1/2 of data set
if isempty(time_start) % if no manual trim defined
indx3 = find(adcp.press_Pa>(median(adcp.press_Pa(round(0.25*end:0.75*end)))...
    -std(adcp.press_Pa(round(0.25*end:0.75*end)))));
time_start=adcp.mtime(indx3(1));
time_end=adcp.mtime(indx3(end));
end

% adcp time is not always constant spacing
dt = nanmedian(diff(adcp.mtime));% data time step (days)
mtime=time_start:dt:time_end;
time_id=find(adcp.mtime>=time_start & adcp.mtime<=time_end); % original mtime index

% get variables on mtime
press_Pa = interp1(adcp.mtime,adcp.press_Pa,mtime,'nearest',NaN);
temperature = interp1(adcp.mtime,adcp.temperature,mtime,'nearest',NaN);
heading = interp1(adcp.mtime,adcp.heading,mtime,'linear',NaN);
pitch = interp1(adcp.mtime,adcp.pitch,mtime,'linear',NaN);
roll = interp1(adcp.mtime,adcp.roll,mtime,'linear',NaN);

%% low pass filter depth, pressure and temp
% large step mtime for filtered data filter_time/2
dt_filt = filter_time/60/24/2;
mtime_filt = time_start:dt_filt:time_end;

% find density and depth
S=34.6; g=9.81;
[~,sigma] = swstate(S,nanmean(temperature),nanmean(press_Pa)*1E-4);
rho = 1000+sigma;
depth = press_Pa./(rho*g)+head_height; % hydrostatic depth above bottom
Depth = nanmean(depth);

% low pass filter
[pressure_filt,Fourier.pressure,~] = low_pass_filter(naninterp(press_Pa),dt*24*3600,filter_time*60);
[depth_filt,Fourier.depth,~] = low_pass_filter(naninterp(depth),dt*24*3600,filter_time*60);
[temperature_filt,Fourier.temperature,~] = low_pass_filter(naninterp(temperature),dt*24*3600,filter_time*60);
[depth_tides,~,~] = low_pass_filter(naninterp(depth),dt*24*3600,2*3600);

depth_filt = interp1(mtime,depth_filt,mtime_filt,'linear',nan);
pressure_filt = interp1(mtime,pressure_filt,mtime_filt,'linear',nan);
temperature_filt = interp1(mtime,temperature_filt,mtime_filt,'linear',nan);
depth_tides = interp1(mtime,depth_tides,mtime_filt,'linear',nan);

%% three beam transformation
if ~isempty(bad_beam) && length(bad_beam)==1
[adcp]=adcp_three_beam_solution(adcp,bad_beam,time_id);
end

%% find average error profile
if ~isempty(bad_beam)
Vel_Error = sqrt(nanmean(adcp.error_vel_orig(:,time_id).^2,2));

else
Vel_Error = sqrt(nanmean(adcp.error_vel(:,time_id).^2,2));    
end

%% trim vertical extent
%calculate bin height
bin_MAB_orig = cfg.ranges(:,1)'+head_height;
bin_z_orig = bin_MAB_orig-Depth;
N=length(depth);

% calculate vertical data to drop
if ~isempty(depth_max) % specified constraint
    Rmax = depth_max;
else % automated calculation
    Rmax = min(depth(round(N/4):round(N*3/4)))*cos(deg2rad(20));
end
bins_drop_top = length(find(bin_MAB_orig>Rmax));

for i = 1:4;
    corr_avg(i,:)=mean(adcp.corr(:,i,time_id),3);
end
corr_avg_min = min(corr_avg,[],1);

if isempty(bottom_trim) % base bottom trim on min adcp.corr and vel_error
bd1 = length(find(corr_avg_min'<115 & bin_MAB_orig<Rmax));
bd2 = 0;%length(find(Vel_Error>0.05 & bin_h_orig<Rmax));
bins_drop_bottom = max(bd1,bd2);
bottom_trim_actual = bin_MAB_orig(bins_drop_bottom);
clear bd1 bd2
else
bins_drop_bottom = length(find(bottom_trim>bin_MAB_orig));
end

bin_MAB = bin_MAB_orig(1+bins_drop_bottom:end-bins_drop_top);
bin_z = bin_MAB-Depth;
bins_good=length(bin_MAB);

%% velocities

% drop vertical bad bins
vel_east2 = adcp.east_vel(1+bins_drop_bottom:end-bins_drop_top,:);
vel_north2 = adcp.north_vel(1+bins_drop_bottom:end-bins_drop_top,:);
vel_vert2 = adcp.vert_vel(1+bins_drop_bottom:end-bins_drop_top,:);
if ~isempty(bad_beam)
    vel_error2 = adcp.error_vel_orig(1+bins_drop_bottom:end-bins_drop_top,:);
else
    vel_error2 = adcp.error_vel(1+bins_drop_bottom:end-bins_drop_top,:);
end

% remove bad correlation data
if bad_beam==1
    kk=[2 3 4];
elseif bad_beam==2
    kk=[1 3 4];
elseif bad_beam==3
    kk=[1 2 4];
elseif bad_beam==4
    kk=[1 2 3];
else
    kk=[1 2 3 4];
end
ii = 1:length(adcp.mtime); % condense tensor
for jj= (1+bins_drop_bottom):(length(bin_MAB_orig)-bins_drop_top)
    corr_max(jj,ii) = max(adcp.corr(jj,kk,ii),[],2);
end

if isempty(corr_min) % if not defined
    corr_min=0;
end
% set low corr values to nan
vel_east2(corr_max < corr_min)=NaN;
vel_north2(corr_max < corr_min)=NaN;
vel_vert2(corr_max < corr_min)=NaN;
vel_error2(corr_max < corr_min)=NaN;

% remove manually specified bad velocity data
for kk=1:length(time_bad_vel)/2;
    indxv = time_bad_vel(kk) < adcp.mtime & time_bad_vel(kk+1) > adcp.mtime;
    for j=1:size(vel_east2,1)
       vel_east2(j,indxv)=nan;
       vel_north2(j,indxv)=nan;
       vel_vert2(j,indxv)=nan;
    end
end

% interpolate to common mtime
for j=1:size(vel_east2,1)
   vel_east(j,:) = interp1(adcp.mtime,vel_east2(j,:),mtime,'nearest',nan);
   vel_north(j,:) = interp1(adcp.mtime,vel_north2(j,:),mtime,'nearest',nan);
   vel_vert(j,:) = interp1(adcp.mtime,vel_vert2(j,:),mtime,'nearest',nan);
   vel_error(j,:) = interp1(adcp.mtime,vel_error2(j,:),mtime,'nearest',nan);
end
clear vel_east2 vel_north2 vel_vert2 vel_error2

%% low pass filter velocity data
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
           
    indxnan = sum(~isnan(u),2) < sample_size; % sample size bigger than 20 
    
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
    vel_east_filt(indxnan,i) = NaN;
    vel_north_filt(indxnan,i) = NaN;
    vel_vert_filt(indxnan,i) = NaN;    
    vel_err_filt(indxnan,i) = NaN;
%     vel_err_east(indxnan,i) = NaN;
%     vel_err_north(indxnan,i) = NaN;
%     vel_err_vert(indxnan,i) = NaN;
    
    clear u v w uerr indxnan indx
 end

%% depth average mean velocities
Vel_East = nanmean(vel_east_filt,1); % depth averaged values
Vel_North=nanmean(vel_north_filt,1);
Vel_Vert=nanmean(vel_vert_filt,1);

Vel_Err = sqrt(nansum(vel_err_filt.^2,1))./size(vel_err_filt,1);
% Vel_East_Err = sqrt(nansum(vel_err_east.^2,1))./size(vel_err_east,1);
% Vel_North_Err = sqrt(nansum(vel_err_north.^2,1))./size(vel_err_north,1);
% Vel_Vert_Err = sqrt(nansum(vel_err_vert.^2,1))./size(vel_err_vert,1);

[~,Fourier.Vel_East,~] = low_pass_filter(naninterp(nanmean(vel_east,1)),dt*24*3600,filter_time*60);
[~,Fourier.Vel_North,~] = low_pass_filter(naninterp(nanmean(vel_north,1)),dt*24*3600,filter_time*60);

%% compute angle of maximum variance for mean currents

% depth avg current direction of max variance
[theta_meanU]=adcp_paa(Vel_East,Vel_North);

% instantaneous u direction of max variance
bin = round(length(bin_MAB/4));
[theta_u]=adcp_paa(vel_east(bin,:),vel_north(bin,:));

if isempty(theta) % if no input theta is defined
theta=theta_meanU;
end
%% rotate coordinate frame 
[vel_maj_filt, vel_min_filt]= adcp_rotation(vel_east_filt,vel_north_filt, theta);
[vel_maj, vel_min]= adcp_rotation(vel_east,vel_north, theta);
[~,Fourier.Vel_Maj,Fourier.s_vel] = low_pass_filter(naninterp(nanmean(vel_maj,1)),dt*24*3600,filter_time*60);
[~,Fourier.Vel_Min,~] = low_pass_filter(naninterp(nanmean(vel_min,1)),dt*24*3600,filter_time*60);

%% Compute extrapolated profiles

if extrap
% compute overall roughness height z0
vel_min_corr = 0.020;
max_depth_fract = 0.5;
zref = 1;
[roughness] = get_z0_log_layer(vel_maj_filt,vel_min_filt,bin_MAB,Depth,max_depth_fract,vel_min_corr,zref);
% [roughness] = get_z0_log_layer(vel_x,vel_y,bin_MAB,Depth,vel_min_corr) 
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
bin_MAB_extrap = linspace(0,max(depth_filt),4*round(Depth/cfg.cell_size(1,1)));
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
if mtime(end)-mtime(1) > 2
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
%% velocity correlations
% umaj_abs_u = naninterp(vel_maj.*sqrt(vel_maj.^2+vel_min.^2));
% umin_abs_u = naninterp(vel_min.*sqrt(vel_maj.^2+vel_min.^2));
% 
% [umaj_abs_u_filt,Fourier.umaj_abs_u,~]=low_pass_filter(umaj_abs_u,dt*24*3600,filter_time*60);
% [umin_abs_u_filt,Fourier.umin_abs_u,~]=low_pass_filter(umin_abs_u,dt*24*3600,filter_time*60);
% 
% umaj_abs_u_filt = interp1(mtime,umaj_abs_u_filt',mtime_filt,'linear',nan)';
% umin_abs_u_filt = interp1(mtime,umin_abs_u_filt',mtime_filt,'linear',nan)';
%  
%% Process Wave Data

if Waves
% Error = Vel_Error(1+bins_drop_bottom:end-bins_drop_top);
bin = round(length(bin_MAB)/8); % use bin near bottom, %find(Error==min(Error)); % select bin with lowest total error
% bin = floor(length(bin_MAB)/2);% select bin at midpoint
% bin = floor(length(bin_MAB)*2/3);% select bin at top
dt_sec = dt*24*3600;
doffu = bin_MAB(bin);
doffp = head_height;
% dirmethod = 1; % use method 1 if direction unknown, use 2 if waves in +x direction.

% use rotated coordinates
Uw = vel_maj(bin,:);
Vw = vel_min(bin,:);
Ww = vel_vert(bin,:);
% use NE coordinates
% Uw = vel_east(bin,:);
% Vw = vel_north(bin,:);

fs = 1./dt_sec;
dt_big_day = 1/24/60*(filter_time);
    
 for i=1:(length(mtime_filt)-2) % only process points with full time series
    tdown = mtime_filt(i+1)-dt_big_day/2;
    tup = mtime_filt(i+1)+dt_big_day/2;
    indx = find(mtime>tdown & mtime<=tup);
        
    if ~isempty(indx)
    U = Uw(indx);
    V = Vw(indx);
    W = Ww(indx);
    P = press_Pa(indx)*1E-4; % pressure is in dbar

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
    u = vel_maj(:,indx);
    v = vel_min(:,indx);
    w = vel_vert(:,indx);
    p = press_Pa(indx)*1E-4; % pressure is in dbar
    
    for jj = 1:size(u,1)  % do for each vertical bin 
        if sum(~isnan(u(jj,:)))<3
           u(jj,:) = zeros(1,length(u(jj,:)));
           v(jj,:) = zeros(1,length(v(jj,:)));
           w(jj,:) = zeros(1,length(w(jj,:)));
         end
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
    
    indxnan = sum(~isnan(u),2) < 30; % sample size bigger than 100 for 10% err
    
    WaveData.uvw_rms(:,i) = sqrt(nanmean(u.^2+v.^2+w.^2,2));
    WaveData.u_rms(:,i) = sqrt(nanmean(u.^2,2));
    WaveData.v_rms(:,i) = sqrt(nanmean(v.^2,2));
    WaveData.w_rms(:,i) = sqrt(nanmean(w.^2,2));
    
    WaveData.uvw_rms(indxnan,i) = NaN;
    WaveData.u_rms(indxnan,i) = NaN;
    WaveData.v_rms(indxnan,i) = NaN;
    WaveData.w_rms(indxnan,i) = NaN;
        
    WaveData.uvwp_rms(:,i) = sqrt(nanmean(up.^2+vp.^2+wp.^2,2));
    WaveData.up_rms(:,i) = sqrt(nanmean(up.^2,2));
    WaveData.vp_rms(:,i) = sqrt(nanmean(vp.^2,2));
    WaveData.wp_rms(:,i) = sqrt(nanmean(wp.^2,2));
    
    WaveData.uvwp_rms(indxnan,i) = NaN;
    WaveData.up_rms(indxnan,i) = NaN;
    WaveData.vp_rms(indxnan,i) = NaN;
    WaveData.wp_rms(indxnan,i) = NaN;
    
    % find stress parameters <u|u|>
    WaveData.uabsu_bar(:,i) = nanmean(u.*sqrt(u.^2+v.^2),2);
    WaveData.vabsu_bar(:,i) = nanmean(v.*sqrt(u.^2+v.^2),2);
    
  
    end
    clear U P V up vp wp u v w
 end
 WaveData.fmt = WaveStats.fmt;
 WaveData.mtime = mtime_filt(2:end-1);
 WaveData.nfft = nfft;
 WaveData.doffp = doffp;
 WaveData.doffu = doffu;
 WaveData.fs = fs;
 WaveData.fcutoff = fcutoff;
 WaveData.df = WaveStats.df;
 WaveData.num_avg = WaveStats.num_avg;
 
 % bottom stress    
 WaveData.Tau_bx = -rho*WaveData.upwp_bar;
 WaveData.Tau_by = -rho*WaveData.vpwp_bar;
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
input.bad_beam=bad_beam;
input.Waves=Waves;
input.nfft=nfft;
input.corr_min=corr_min;

Misc.bins_drop_top = bins_drop_top;
Misc.bins_good = bins_good;
Misc.bins_drop_bottom = bins_drop_bottom;
% Misc.beta=beta;
Misc.ens = ens;
Misc.hdr = hdr;

%% save files
clear i j ii bottom_trim depth_max bad_beam Waves corr_min...
    N filename filenamerdr files filesrdr indx indx2 indx3 indx_corr jj kk...
    time_start time_end time_id bins_drop_top bins_good ens hdr...
    bins_drop_bottom beta

%%
save('ADCP_PostProcess.mat','-v7.3');


