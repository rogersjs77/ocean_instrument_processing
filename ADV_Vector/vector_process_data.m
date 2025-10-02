% function [ ] = vector_process_data(filter_time,P_atm,Waves_Spectra,Waves_Burst,nfft,...
%     fcutoff,theta,dirmethod,head_height,case_height,vertical_orientation,...
%     x_heading,time_start,time_end,extrap,sample_size,corr_min,time_bad_vel)
%
% this file takes ADV data,
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
% theta = angle in degrees to rotate coordinates, if [] finds maximum U variance
% P_atm = atmospheric pressure in Pa - remember recorded pressure in ADCP
% is in dPa = 10*Pa;
% head_height is height of ADV head above bottom (m)
% case_height is height of ADV body above bottom (m)
% vertical_orientation = 'up' or 'down', looking upward 'up', looking downward 'down'
% x_heading =  heading of x axis CW from N (N=0, E=90, S=180,W=270)
% vel_min_corr = minimum velocity to use for roughness calc, typ 1-5cm/s
%
% frequency resolution is: df = fs/(nfft-1); fs is sampling frequency
% number of averages num = floor(4*n/nfft)-3; n is # of points in dt_big
%%

function [ ] = vector_process_data(filter_time,P_atm,Waves_Spectra,Waves_Burst,nfft,...
    fcutoff,theta,dirmethod,head_height,case_height,vertical_orientation,...
    x_heading,time_start,time_end,extrap,sample_size,corr_min,time_bad_vel)

files=dir('*.hdr');   

for i=1:length(files)

files=dir('*.hdr');
    
[~,filename,~]=fileparts(files(i).name);

%% load data
adv.hdr = importdata(strcat(filename,'.hdr'),'');
adv.sample_rate = cell2mat(textscan(char(adv.hdr(11)),'Sampling rate   %f %*s'));
adv.burst_interval = cell2mat(textscan(char(adv.hdr(13)),'Burst interval   %f %*s'));
[m] = textscan(char(adv.hdr(7)),'Time of first measurement %s %s %s');
adv.time_first = datenum([cell2mat(m{1}) ' ' cell2mat(m{2}) ' ' cell2mat(m{3})]);
[m] = textscan(char(adv.hdr(8)),'Time of last measurement %s %s %s');
adv.time_last = datenum([cell2mat(m{1}) ' ' cell2mat(m{2}) ' ' cell2mat(m{3})]);

DATA = dlmread(strcat(filename,'.dat'));
adv.data_burst_no = DATA(:,1)';
data_ens_no = DATA(:,2)';
adv.data_vel_x = DATA(:,3)';
adv.data_vel_y = DATA(:,4)';
adv.data_vel_z = DATA(:,5)';
data_amp1 =  DATA(:,6)';
data_amp2 =  DATA(:,7)';
data_amp3 =  DATA(:,8)';
data_corr1 =  DATA(:,12)';
data_corr2 =  DATA(:,13)';
data_corr3 =  DATA(:,14)';
adv.data_pressure = DATA(:,15)';
clear DATA

DATA =  dlmread(strcat(filename,'.vhd',' '));
adv.burst_time_start = datenum(DATA(:,3),DATA(:,1),DATA(:,2),...
    DATA(:,4),DATA(:,5),DATA(:,6))';
adv.burst_no = DATA(:,7)';
adv.burst_samples = DATA(:,8)';

if ~isempty(adv.burst_interval) % burst sampling
% burst time should be middle of burst
adv.mtime = adv.burst_time_start+0.5*nanmean(adv.burst_samples)/adv.sample_rate/(24*3600);
% process burst data to mean
adv.data_mtime = NaN+adv.data_vel_x;
else % continuous sampling
adv.time_last2 = adv.time_first+(1/(adv.sample_rate*24*3600))*(length(adv.data_vel_x)-1);    
adv.mtime = adv.time_first:1/(adv.sample_rate*24*3600):adv.time_last2;
% process burst data to mean
adv.data_mtime = adv.mtime;    
end

if ~isempty(adv.burst_interval) % burst sampling
    for j = 1:length(adv.burst_time_start)
        indx = find(adv.data_burst_no == adv.burst_no(j));
        if ~isempty(indx)
        adv.data_mtime(indx) = adv.burst_time_start(j) + ...
            (data_ens_no(indx)-1)./adv.sample_rate/(24*3600);
        adv.vel_x(j) = nanmean(adv.data_vel_x(indx)); 
        adv.vel_y(j) = nanmean(adv.data_vel_y(indx));
        adv.vel_z(j) = nanmean(adv.data_vel_z(indx));
        adv.amp1(j) = nanmean(data_amp1(indx));
        adv.amp2(j) = nanmean(data_amp2(indx));
        adv.amp3(j) = nanmean(data_amp3(indx));
        adv.corr1(j) = nanmean(data_corr1(indx));
        adv.corr2(j) = nanmean(data_corr2(indx));
        adv.corr3(j) = nanmean(data_corr3(indx));
        adv.pressure(j) = nanmean(adv.data_pressure(indx));
        adv.pressure_diff(j) = max(adv.data_pressure(indx))-min(adv.data_pressure(indx));
        else % if not defined use NaN
        adv.vel_x(j) = NaN; 
        adv.vel_y(j) = NaN;
        adv.vel_z(j) = NaN;
        adv.amp1(j) = NaN;
        adv.amp2(j) = NaN;
        adv.amp3(j) = NaN;
        adv.corr1(j) = NaN;
        adv.corr2(j) = NaN;
        adv.corr3(j) = NaN;
        adv.pressure(j) = NaN;
        adv.pressure_diff(j) = NaN;    
        end

    end
else % continuous sampling
    adv.vel_x = adv.data_vel_x;
    adv.vel_y = adv.data_vel_y;
    adv.vel_z = adv.data_vel_z;
    adv.amp1 = data_amp1;
    adv.amp2 = data_amp2;
    adv.amp3 = data_amp3;
    adv.corr1 = data_corr1;
    adv.corr2 = data_corr2;
    adv.corr3 = data_corr3;
    adv.pressure = adv.data_pressure;
    adv.pressure_diff = [NaN diff(adv.data_pressure)];   
end

clear DATA
DATA =  dlmread(strcat(filename,'.sen',' '));
sen_time = datenum(DATA(:,3),DATA(:,1),DATA(:,2),...
    DATA(:,4),DATA(:,5),DATA(:,6))';
% downsample temp and other less important variables to adv.mtime
adv.temp = interp1(sen_time,DATA(:,14)',adv.mtime,'linear',NaN);
adv.heading_body = interp1(sen_time,DATA(:,11)',adv.mtime,'linear',NaN);
adv.pitch_body = interp1(sen_time,DATA(:,12)',adv.mtime,'linear',NaN);
adv.roll_body = interp1(sen_time,DATA(:,13)',adv.mtime,'linear',NaN);
adv.voltage = interp1(sen_time,DATA(:,9)',adv.mtime,'linear',NaN);

% delete high res corr and amp data
clear DATA DATA2 DATA3 data_corr1 data_corr2 data_corr3...
    data_amp1 data_amp2 data_amp3 sen_time j indx ...
    data_ens_no

%% Rotate xyz coordinate to earth
% heading is degrees CW from North of Y axis
% pitch is typically 0
% roll is 0 for downward facing, 180 for upward facing

if strcmp(vertical_orientation,'up')
    roll = 180;
    pitch= 0;
    heading = x_heading+90;
elseif strcmp(vertical_orientation,'down')
    roll = 0;
    pitch= 0;
    heading = x_heading-90;
end
adv=adv_XYZ_to_ENU(adv,heading,pitch,roll); 

%% trim mtime to time underwater

if isempty(P_atm) % if no manually defined atmospheric P
    % find avg atmospheric pressure within expected range < 1.2 dbar
    P_atm = 1E4*nanmedian(adv.pressure(adv.pressure<1.2 & adv.pressure>0));
    if isempty(P_atm) || isnan(P_atm)
      P_atm = 0;
    end
end

% find trim time from pressure within 1std of mean
adv.press_Pa = 1E4*adv.pressure-P_atm;
indx2 = adv.press_Pa>1E7; % remove unrealistic pressures
adv.press_Pa(indx2)=nan;
adv.press_Pa=naninterp(adv.press_Pa); % this is pressure at head sensor

if isempty(time_start) % if no manual trim defined
% find pressures outside 1 std of middle 1/2 of data set
indx3 = find(adv.press_Pa>(median(adv.press_Pa(0.25*end:0.75*end))...
    -std(adv.press_Pa(0.25*end:0.75*end))));
time_start=adv.mtime(indx3(1));
time_end=adv.mtime(indx3(end));
end

% adp time is not always constant spacing
if ~isempty(adv.burst_interval) % burst sampling
mtime=time_start:adv.burst_interval/3600/24:time_end;
else
mtime=time_start:1/adv.sample_rate/3600/24:time_end;    
end
dt = mtime(2)-mtime(1); % data time step (days)

% nearest interpolation method
press_Pa = interp1(adv.mtime,adv.press_Pa,mtime,'nearst',NaN);
temperature = interp1(adv.mtime,adv.temp,mtime,'nearest',NaN);

% find density and depth
S=34.6; g=9.81;
[~,sigma] = swstate(S,nanmean(temperature),nanmean(press_Pa)*1E-4);
rho = 1000+sigma;
depth = press_Pa./(rho*g)+case_height; % hydrostatic depth above bottom in m
Depth = nanmean(depth);

% low pass filter pressure and velocity
dt_filt = filter_time/60/24/2;
mtime_filt = time_start:dt_filt:time_end;

[pressure_filt,Fourier.pressure,Fourier.s] = low_pass_filter(press_Pa,dt*24*3600,filter_time*60);
[depth_filt,Fourier.depth,~] = low_pass_filter(depth,dt*24*3600,filter_time*60);
[temperature_filt,Fourier.temperature,~] = low_pass_filter(temperature,dt*24*3600,filter_time*60);
[depth_tides,~,~] = low_pass_filter(depth,dt*24*3600,2*3600);

depth_filt = interp1(mtime,depth_filt,mtime_filt,'linear');
pressure_filt = interp1(mtime,pressure_filt,mtime_filt,'linear');
temperature_filt = interp1(mtime,temperature_filt,mtime_filt,'linear');
depth_tides = interp1(mtime,depth_tides,mtime_filt,'linear');

%% clean up velocity data

vel_east = interp1(adv.mtime,adv.vel_east,mtime,'nearest',nan);
vel_north = interp1(adv.mtime,adv.vel_north,mtime,'nearest',nan);
vel_vert = interp1(adv.mtime,adv.vel_vert,mtime,'nearest',nan);

% remove low correlaton data, and 
% corr_min = 60;
corr_avg = interp1(adv.mtime,nanmin([adv.corr1; adv.corr2; adv.corr3]),...
    mtime,'nearest',nan);
indx = corr_avg < corr_min;
vel_east(indx)=NaN;
vel_north(indx)=NaN;
vel_vert(indx)=NaN;

% remove large outlyers 8 std
indx = vel_east < nanmedian(vel_east)-5*nanstd(vel_east) |...
    vel_east > nanmedian(vel_east)+5*nanstd(vel_east);
vel_east(indx)=NaN;
indx = vel_north < nanmedian(vel_north)-5*nanstd(vel_north) |...
    vel_north > nanmedian(vel_north)+5*nanstd(vel_north);
vel_north(indx)=NaN;
indx = vel_vert < nanmedian(vel_vert)-8*nanstd(vel_vert) |...
    vel_vert > nanmedian(vel_vert)+8*nanstd(vel_vert);
vel_vert(indx)=NaN;

% remove manually specified bad velocity data
for kk=1:length(time_bad_vel)/2;
    indxv = time_bad_vel(kk) < mtime & time_bad_vel(kk+1) > mtime;
       vel_east(indxv)=nan;
       vel_north(indxv)=nan;
       vel_vert(indxv)=nan;
end

% calculate error, 0.5% of value +/- 1 mm/s
vel_error = 0.005*sqrt(vel_east.^2+vel_north.^2+vel_vert.^2)+1/1000;

%% low pass filter velocity data
% [vel_east_filt,Fourier.Vel_East,~]=low_pass_filter(vel_east,dt*24*3600,filter_time*60);
% [vel_north_filt,Fourier.Vel_North,~]=low_pass_filter(vel_north,dt*24*3600,filter_time*60);
% [vel_vert_filt,Fourier.Vel_Vert,~]=low_pass_filter(vel_vert,dt*24*3600,filter_time*60);
[corr_avg_filt,~,~]=low_pass_filter(corr_avg,dt*24*3600,filter_time*60);
corr_avg_filt = interp1(mtime,corr_avg_filt',mtime_filt,'linear')';

% vel_east_filt = interp1(mtime,vel_east_filt',mtime_filt,'linear')';
% vel_north_filt = interp1(mtime,vel_north_filt',mtime_filt,'linear')';
% vel_vert_filt = interp1(mtime,vel_vert_filt',mtime_filt,'linear')';


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

    clear u v w uerr indxnan indx
 end

%% depth average mean velocities (keep for consistency w ADCP/ADP data)
Vel_East = vel_east_filt; % depth averaged values
Vel_North= vel_north_filt;
Vel_Vert = vel_vert_filt;
Vel_Err = vel_err_filt;

%% compute angle of maximum variance for mean currents
% depth avg current direction of max variance
[theta_meanU]=adcp_paa(Vel_East,Vel_North);

% instantaneous u direction of max variance
[theta_u]=adcp_paa(vel_east,vel_north);

if isempty(theta) % if no input theta is defined
theta=theta_meanU;
end

%% rotate coordinate frame
[vel_maj_filt, vel_min_filt]= adcp_rotation(vel_east_filt,vel_north_filt, theta);
[vel_maj, vel_min]= adcp_rotation(vel_east,vel_north, theta);

[~,Fourier.Vel_Maj,Fourier.s_vel] = low_pass_filter(naninterp(vel_maj),dt*24*3600,filter_time*60);
[~,Fourier.Vel_Min,~] = low_pass_filter(naninterp(vel_min),dt*24*3600,filter_time*60);

%% estimate roughness height and depth avg flow
if extrap
vel_min_corr = 0.02;
[roughness] = get_z0_reynolds_stress(vel_maj,vel_min,vel_vert,...
    mtime,mtime_filt,head_height,Depth,vel_min_corr);

%%% !! REVIEW z0 RESULTS BEFORE SELECTING !! %%%
z0_u = roughness.z0.med;
z0_v = roughness.z0.med;

[Vel_East,Vel_North,ustar,vstar] = ...
    get_Ubar_point(vel_east_filt,vel_north_filt,z0_u,z0_v,head_height,depth_filt);

[Vel_Maj,Vel_Min,ustar,vstar] = ...
    get_Ubar_point(vel_maj_filt,vel_min_filt,z0_u,z0_v,head_height,depth_filt);

else
   Vel_East = vel_east_filt;
   Vel_North = vel_north_filt;
   Vel_Maj = vel_maj_filt;
   Vel_Min = vel_min_filt;    
end

%% tidal analysis
eta = depth-nanmean(depth);
eta_filt = depth_filt-nanmean(depth_filt);
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

%% Process Wave Data

if Waves_Spectra
bin = 1; % adv is only one    
dt_sec = dt*24*3600;
doffu = head_height;
doffp = case_height;

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
    
    % remove nans (dont interpolate) since we only want p-u correlation
%     indx2 = ~isnan(u); 
%     u = u(indx2);
%     v = v(indx2);
%     w = w(indx2);
%     p = p(indx2);
%     
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
    
    WaveData.uvw_rms(:,i) = sqrt(nanmean(u.^2+v.^2+w.^2,2));
    WaveData.u_rms(:,i) = sqrt(nanmean(u.^2,2));
    WaveData.v_rms(:,i) = sqrt(nanmean(v.^2,2));
    WaveData.w_rms(:,i) = sqrt(nanmean(w.^2,2));
        
    WaveData.uvwp_rms(:,i) = sqrt(nanmean(up.^2+vp.^2+wp.^2,2));
    WaveData.up_rms(:,i) = sqrt(nanmean(up.^2,2));
    WaveData.vp_rms(:,i) = sqrt(nanmean(vp.^2,2));
    WaveData.wp_rms(:,i) = sqrt(nanmean(wp.^2,2));
    
    % find stress parameters <u|u|>
    WaveData.uabsu_bar(:,i) = nanmean(u.*sqrt(u.^2+v.^2));
    WaveData.vabsu_bar(:,i) = nanmean(v.*sqrt(u.^2+v.^2));
 
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
   
 clear WaveStats doffp doffu fs fcutoff bin Error U P V W

 % Lagrangian mean flows in xy direction
 Vel_Maj_Stokes = interp1(WaveData.mtime,WaveData.Us,mtime_filt,'linear',NaN);
 Vel_Min_Stokes = interp1(WaveData.mtime,WaveData.Vs,mtime_filt,'linear',NaN);

 Vel_Maj_Lagr = Vel_Maj+Vel_Maj_Stokes;
 Vel_Min_Lagr = Vel_Min+Vel_Min_Stokes;
% rotate these back to NE coordinates
 [Vel_East_Stokes,Vel_North_Stokes]= adcp_rotation(Vel_Maj_Stokes,Vel_Min_Stokes, -theta);
 [Vel_East_Lagr,Vel_North_Lagr]= adcp_rotation(Vel_Maj_Lagr,Vel_Min_Lagr, -theta);

end

%% try burst waves

if Waves_Burst
bin = 1; % adv is only one    
% dt_sec = dt*24*3600;
% doffu = head_height;
% doffp = case_height;
    
    
doffu = head_height;
doffp = case_height;% assumed wh sentinel casing;
dt_big_day = 1/24/60*(filter_time);
% Tm = 10;
% Tstd = 2;

fs = adv.sample_rate;

 for i=1:(length(mtime_filt)-2) % only process points with full time series   
    tup = mtime_filt(i+1)+dt_big_day/2;
    tdown = mtime_filt(i+1)-dt_big_day/2;
    indx = find(adv.data_mtime>tdown & adv.data_mtime<=tup);
     
    p = adv.data_pressure(indx);
    u = adv.data_vel_east(indx);
    v = adv.data_vel_north(indx);
    w = adv.data_vel_vert(indx);
    if ~isempty(u)
      [WaveStats] = wave_stats_spectra(u,v,p,nfft,doffu,doffp,fs,fcutoff,rho,dirmethod);

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
        
        
%     [WaveStats] = wave_stats_burst(u,v,p,Tm,Tstd,doffp,rho);            
%     WaveData.depth(i) = WaveStats.depth;
%     WaveData.Hsig_ss(i) = WaveStats.Hsig_ss;
%     WaveData.Hsig_ss_std(i) = WaveStats.Hsig_ss_std;
%     WaveData.dir_ss(i) = WaveStats.dir_ss;
%     WaveData.spread_ss(i) = WaveStats.spread_ss; 
%     WaveData.Tm_ss(i) = WaveStats.Tm_ss;
%     WaveData.Tm_ss_std(i) = WaveStats.Tm_ss_std;
%     WaveData.Energy_ss(i,:) = WaveStats.Energy_ss;
%     WaveData.Us(i) = WaveStats.Us;
%     WaveData.Vs(i) = WaveStats.Vs;
    
    for jj = 1  % do for each vertical bin 
         if sum(~isnan(u(jj,:)))<3
           u(jj,:) = zeros(1,length(u(jj,:)));
           v(jj,:) = zeros(1,length(v(jj,:)));
           w(jj,:) = zeros(1,length(w(jj,:)));
         end
        
    up = u - nanmean(u);
    vp = v - nanmean(v);
    wp = w - nanmean(w);
    
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
    WaveData.uabsu_bar(:,i) = nanmean(u.*sqrt(u.^2+v.^2));
    WaveData.vabsu_bar(:,i) = nanmean(v.*sqrt(u.^2+v.^2));
    end    
    WaveData.mtime(i) = mtime_filt(i);    
 end
 WaveData.fmt = WaveStats.fmt;
 WaveData.mtime = mtime_filt(2:end-1);
 WaveData.nfft = nfft;

 WaveData.fs = fs;
 WaveData.fcutoff = fcutoff;
 WaveData.df = WaveStats.df;
 WaveData.num_avg = WaveStats.num_avg;
 WaveData.doffp = doffp;
 WaveData.doffu = doffu;
 
 % bottom stress    
 WaveData.Tau_bx = -rho*WaveData.upwp_bar;
 WaveData.Tau_by = -rho*WaveData.vpwp_bar;
 WaveData.Tau = sqrt(WaveData.Tau_bx.^2+WaveData.Tau_by.^2);
 
 % wave dir relative to minor axis (typical beach scenario)
 WaveData.dir_Min_ss = -WaveData.dir_ss+90+theta; 
 
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

 clear WaveStats doffp doffu fs fcutoff bin Error U P V u v p   indx tup tdown 
 
 % Lagrangian mean flows
 Vel_East_Stokes = interp1(WaveData.mtime,WaveData.Us,mtime_filt,'linear',NaN);
 Vel_North_Stokes = interp1(WaveData.mtime,WaveData.Vs,mtime_filt,'linear',NaN);
 % depth average
 Vel_East_Lagr = Vel_East+Vel_East_Stokes;
 Vel_North_Lagr = Vel_North+Vel_North_Stokes;
 % point measurement - needs work
 vel_east_lagr = vel_east_filt+Vel_East_Stokes;
 vel_north_lagr = vel_north_filt+Vel_North_Stokes;

 % rotate to major/minor axes
 [Vel_Maj_Stokes,Vel_Min_Stokes]= adcp_rotation(Vel_East_Stokes,Vel_North_Stokes, theta);
 [Vel_Maj_Lagr,Vel_Min_Lagr]= adcp_rotation(Vel_East_Lagr,Vel_North_Lagr, theta);
end

%% clean up data
input.P_atm = P_atm;
input.filter_time=filter_time;
input.case_height=case_height;
input.head_height=head_height;
input.Waves_Spectra=Waves_Spectra;
input.Waves_Burst = Waves_Burst;
input.nfft=nfft;
input.theta = theta;
input.vertical_orientation = vertical_orientation;
input.x_heading=x_heading;

%% save files
clear i j ii bottom_trim depth_max bad_beam Waves nfft corr_min...
    N indx indx2 indx3 indx_corr jj kk...
    time_start time_end bins_drop_top bins_good ...
    bins_drop_bottom beta 

save('ADV_PostProcess.mat','-v7.3');
end

end