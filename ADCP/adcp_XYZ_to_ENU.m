% function [adcp] = adcp_XYZ_to_ENU(adcp,time_id)
% 
% Convert ADCP velocities from xyz to earth coordinates
% Justin Rogers & Kara Scheu
% Stanford EFML, 3/2013
%
% INPUT
% adcp = standard adcp structure file from rdradcp in xyz coordinates
% time_id = optional index of mtime to use for heading pitch and roll averaqes
%           if empty [], uses full time series for averages 
% 
% OUTPUT
% adcp = adcp structure file with velocities in earth coordinates
%        velocities (east_vel, north_vel, vert_vel, error_vel)
%
%%
function [adcp] = adcp_XYZ_to_ENU(adcp,time_id)

if isempty(time_id)
    time_id=~isnan(adcp.mtime);
end

CH = cos(deg2rad(nanmean(adcp.heading(time_id))));
SH = sin(deg2rad(nanmean(adcp.heading(time_id))));
CP = cos(deg2rad(nanmean(adcp.pitch(time_id))));
SP = sin(deg2rad(nanmean(adcp.pitch(time_id))));
CR = cos(deg2rad(nanmean(adcp.roll(time_id)+180)));
SR = sin(deg2rad(nanmean(adcp.roll(time_id)+180)));

H = [CH SH 0 ;-SH CH 0 ; 0 0 1];
P = [1 0 0; 0 CP -SP; 0 SP CP];
R = [CR 0 SR; 0 1 0; -SR 0 CR];
M = H*P*R; % this is the transformation matrix
Minv = inv(M);

% compute transformed matrices and assign values
for i = 1:size(adcp.east_vel,1)

    V_temp = M*[adcp.east_vel(i,:);adcp.north_vel(i,:);...
        adcp.vert_vel(i,:)];
    adcp.east_vel(i,:) = V_temp(1,:);
    adcp.north_vel(i,:) = V_temp(2,:);
    adcp.vert_vel(i,:) = V_temp(3,:);
    
end

adcp.config(1,1).coord_sys = 'earth';
adcp.coord_sys = 'earth';
adcp.coord_expl = 'east_vel=east, north_vel=north vert_vel=vert, error_vel = error_up';

end