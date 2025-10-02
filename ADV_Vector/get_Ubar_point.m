function [Vel_East_Avg,Vel_North_Avg,ustar,vstar] = ...
    get_Ubar_point(vel_east_filt,vel_north_filt,z0_u,z0_v,head_height,depth_filt)
% this function finds Ustar, and depth average flow assuming a log layer profile
% given filtered velocites at a point z and roughness height z0
% Justin Rogers, Stanford EFML, September 2013

k = 0.41; % von karman const
depth_filt = depth_filt';

if head_height < z0_u || head_height < z0_v
    ustar = NaN+depth_filt;
    Vel_East_Avg  = NaN+depth_filt;
    vstar = NaN+depth_filt;
    Vel_North_Avg  = NaN+depth_filt;
    disp('avg depth not valid, z0 > head height')
    return
end

% calculate ustar based on avg roughness
ustar = vel_east_filt*k./log(head_height/z0_u);
vstar = vel_north_filt*k./log(head_height/z0_v);

Vel_East_Avg = ustar.*depth_filt./((depth_filt-z0_u)*k).*...
    (log(depth_filt/z0_u)-1+z0_u./depth_filt);
Vel_North_Avg = vstar.*depth_filt./((depth_filt-z0_v)*k).*...
    (log(depth_filt/z0_v)-1+z0_v./depth_filt);


end