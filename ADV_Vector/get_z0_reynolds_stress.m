function [roughness] = get_z0_reynolds_stress(vel_maj,vel_min,vel_vert,...
    mtime,mtime_filt,ref_height,Depth,vel_min_corr)
% This function finds roughnes height z0, and drag coefficient Cd,
% using reynolds stress method given instantaneous velocity
% removing low velocity and unrealistic values
%
% Justin Rogers, Stanford EFML, Semptember 2013
%

k = 0.41; % von karman const
max_depth_fract = 0.6; % fraction of depth allowed for z0;
min_z0 = 1E-4; % minumum allowed z0

%% calculations

% find low pass filtered u
dt = mtime(2)-mtime(1);
filter_time = 2*(mtime_filt(2)-mtime_filt(1))*24*60;
vel_maj_filt = low_pass_filter(naninterp(vel_maj),dt*24*3600,filter_time*60);
vel_min_filt = low_pass_filter(naninterp(vel_min),dt*24*3600,filter_time*60);
vel_vert_filt = low_pass_filter(naninterp(vel_vert),dt*24*3600,filter_time*60);

% find reynolds stresses
up = vel_maj - vel_maj_filt;
vp = vel_min - vel_min_filt;
wp = vel_vert - vel_vert_filt;
roughness.up_wp = low_pass_filter(naninterp(up.*wp),dt*24*3600,filter_time*60);
roughness.up_wp = interp1(mtime,roughness.up_wp,mtime_filt,'linear',NaN);
roughness.vp_wp = low_pass_filter(naninterp(vp.*wp),dt*24*3600,filter_time*60);
roughness.vp_wp = interp1(mtime,roughness.vp_wp,mtime_filt,'linear',NaN);

% find ustar
roughness.ustar_cov = sqrt(abs(roughness.up_wp));
roughness.vstar_cov = sqrt(abs(roughness.vp_wp));

% find roughness height
vel_maj_filt = interp1(mtime,vel_maj_filt,mtime_filt,'linear',NaN);
vel_min_filt = interp1(mtime,vel_min_filt,mtime_filt,'linear',NaN);
roughness.z0.ut = ref_height ./ exp(vel_maj_filt*k./roughness.ustar_cov);
roughness.z0.vt = ref_height ./ exp(vel_min_filt*k./roughness.vstar_cov);

% remove low velocity values
indx = abs(vel_maj_filt) < vel_min_corr;
roughness.z0.ut(indx) = NaN;
indx = abs(vel_min_filt) < vel_min_corr;
roughness.z0.vt(indx) = NaN;

% remove unrealistic results
indx = roughness.z0.ut > max_depth_fract*Depth/30 | roughness.z0.ut < min_z0;
roughness.z0.ut(indx) = NaN;
indx = roughness.z0.vt > max_depth_fract*Depth/30 | roughness.z0.vt < min_z0;
roughness.z0.vt(indx) = NaN;

% compute drag coeff Cd
roughness.Cd.ut = k^2./(log(ref_height./roughness.z0.ut)).^2;
roughness.Cd.vt = k^2./(log(ref_height./roughness.z0.vt)).^2;

% compute statistics for z0
roughness.z0.mean = nanmean([roughness.z0.ut, roughness.z0.vt]);
roughness.z0.u_mean = nanmean(roughness.z0.ut);
roughness.z0.v_mean = nanmean(roughness.z0.vt);

roughness.z0.med = nanmedian([roughness.z0.ut, roughness.z0.vt]);
roughness.z0.u_med = nanmedian(roughness.z0.ut);
roughness.z0.v_med = nanmedian(roughness.z0.vt);

roughness.z0.std = nanstd([roughness.z0.ut, roughness.z0.vt]);
roughness.z0.u_std = nanstd(roughness.z0.ut);
roughness.z0.v_std = nanstd(roughness.z0.vt);

roughness.z0.perc_bad_u = sum(isnan(roughness.z0.ut))/length(roughness.z0.ut)*100;
roughness.z0.perc_bad_v = sum(isnan(roughness.z0.vt))/length(roughness.z0.vt)*100;
roughness.z0.perc_bad_mag = sum(isnan([roughness.z0.ut, roughness.z0.vt]))/...
    length([roughness.z0.ut, roughness.z0.vt])*100;

% compute statistics for Cd
roughness.Cd.mean = nanmean([roughness.Cd.ut, roughness.Cd.vt]);
roughness.Cd.u_mean = nanmean(roughness.Cd.ut);
roughness.Cd.v_mean = nanmean(roughness.Cd.vt);

roughness.Cd.med = nanmedian([roughness.Cd.ut, roughness.Cd.vt]);
roughness.Cd.u_med = nanmedian(roughness.Cd.ut);
roughness.Cd.v_med = nanmedian(roughness.Cd.vt);

roughness.Cd.std = nanstd([roughness.Cd.ut, roughness.Cd.vt]);
roughness.Cd.u_std = nanstd(roughness.Cd.ut);
roughness.Cd.v_std = nanstd(roughness.Cd.vt);

roughness.Cd.perc_bad_u = sum(isnan(roughness.Cd.ut))/length(roughness.Cd.ut)*100;
roughness.Cd.perc_bad_v = sum(isnan(roughness.Cd.vt))/length(roughness.Cd.vt)*100;
roughness.Cd.perc_bad_mag = sum(isnan([roughness.Cd.ut, roughness.Cd.vt]))/...
    length([roughness.Cd.ut, roughness.Cd.vt])*100;

end