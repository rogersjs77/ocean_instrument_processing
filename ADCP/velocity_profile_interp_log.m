% function [vel_maj_extrap,vel_min_extrap,vel_vert_extrap,ustar,vstar] = ...
%     velocity_profile_interp_log(vel_maj_filt,vel_min_filt,vel_vert_filt,z0_u,z0_v,...
%     bin_MAB,bin_MAB_extrap,depth_filt)
%
% This function extrapolates data from bottom to surface given velocities
% in the middle of the water column and bottom roughness height z0,
% top velocities are extended to surface
% bottom velocity bin is fit to a log profile to the bed for u and v
% and linear fit for w
%
% Justin Rogers, Stanford EFML, September 2013
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vel_maj_extrap,vel_min_extrap,vel_vert_extrap,ustar,vstar] = ...
    velocity_profile_interp_log(vel_maj_filt,vel_min_filt,vel_vert_filt,z0_u,z0_v,...
    d_u,d_v,bin_MAB,bin_MAB_extrap,depth_filt)

k = 0.41; % von karman constant

% find high and low bins
indx = bin_MAB_extrap >  bin_MAB(end);
bin_top = bin_MAB_extrap(indx);
indx = bin_MAB_extrap <  bin_MAB(1);
bin_bott = bin_MAB_extrap(indx);

temph = [bin_bott bin_MAB bin_top];
temph_vert = [z0_u bin_MAB bin_top];

% extend top velocity to surface
[vtop_maj,~] = meshgrid(vel_maj_filt(end,:),bin_top);
[vtop_min,~] = meshgrid(vel_min_filt(end,:),bin_top);
[vtop_vert,~] = meshgrid(vel_vert_filt(end,:),bin_top);

% calculate ustar based on avg roughness
ustar = vel_maj_filt(1,:)*k./log((bin_MAB(1)-d_u)/z0_u);
vstar = vel_min_filt(1,:)*k./log((bin_MAB(1)-d_v)/z0_v);

for i = 1:size(vel_maj_filt,2)
   
    vll_maj = ustar(i)/k*log((bin_bott-d_u)./z0_u)';
    vll_min = vstar(i)/k*log((bin_bott-d_v)./z0_v)';
    vll_vert = 0;    
    
    indx = ~isfinite(vll_maj);
    vll_maj(indx) = NaN;
    indx = ~isfinite(vll_min);
    vll_min(indx) = NaN;
   
    % combine three parts and interpolate
    tempv_maj = [vll_maj; vel_maj_filt(:,i); vtop_maj(:,i)];
    tempv_min = [vll_min; vel_min_filt(:,i); vtop_min(:,i)];
    tempv_vert = [vll_vert; vel_vert_filt(:,i); vtop_vert(:,i)];
    
    vel_maj_extrap(:,i) = interp1(temph,tempv_maj,bin_MAB_extrap,'linear','extrap');
    vel_min_extrap(:,i) = interp1(temph,tempv_min,bin_MAB_extrap,'linear','extrap');
    vel_vert_extrap(:,i) = interp1(temph_vert,tempv_vert,bin_MAB_extrap,'linear','extrap');
end

% remove values above surface
[DEPTH,BIN] = meshgrid(depth_filt,bin_MAB_extrap);

indx = DEPTH < BIN;
vel_maj_extrap(indx)=NaN;
vel_min_extrap(indx)=NaN;
vel_vert_extrap(indx)=NaN;

% remove values below z0 + d
Z0_U = z0_u + d_u + 0*DEPTH;
Z0_V = z0_v + d_v + 0*DEPTH;

indx =  Z0_U > BIN;
vel_maj_extrap(indx)=NaN;
vel_vert_extrap(indx)=NaN;

indx =  Z0_V > BIN;
vel_min_extrap(indx)=NaN;

end
