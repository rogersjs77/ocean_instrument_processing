% function [adv] = adv_XYZ_to_ENU(adv,heading,pitch,roll)
% 
% Convert adv velocities from xyz to earth coordinates
% Justin Rogers & Kara Scheu
% Stanford EFML, 8/2013
%
% INPUT
% adv = standard adcp structure file from vector_process_data in xyz coordinates
% heading = degrees CW from N of y axis
% pitch = standard definition
% roll = 0 deg for doward pointing, 180 deg for upward facing
% 
% OUTPUT
% adv = adv structure file with velocities in earth coordinates
%        velocities (east_vel, north_vel, vert_vel, error_vel)
%
%%
function [adv] = adv_XYZ_to_ENU(adv,heading,pitch,roll)

CH = cos(deg2rad(heading));
SH = sin(deg2rad(heading));
CP = cos(deg2rad(pitch));
SP = sin(deg2rad(pitch));
CR = cos(deg2rad(roll));
SR = sin(deg2rad(roll));

% CR = cos(deg2rad(roll+180));
% SR = sin(deg2rad(roll+180));

H = [CH SH 0 ;-SH CH 0 ; 0 0 1];
P = [1 0 0; 0 CP -SP; 0 SP CP];
R = [CR 0 SR; 0 1 0; -SR 0 CR];
M = H*P*R; % this is the transformation matrix
Minv = inv(M);

% compute transformed matrices and assign values
V_temp = M*[adv.vel_x;adv.vel_y;adv.vel_z];
adv.vel_east = V_temp(1,:);
adv.vel_north = V_temp(2,:);
adv.vel_vert = V_temp(3,:);

if isfield(adv,'data_vel_x')
% compute transformed matrices and assign values
V_temp = M*[adv.data_vel_x;adv.data_vel_y;adv.data_vel_z];
adv.data_vel_east = V_temp(1,:);
adv.data_vel_north = V_temp(2,:);
adv.data_vel_vert = V_temp(3,:);    
end
    
adv.coord_sys = 'earth';
adv.heading = heading;
adv.pitch = pitch;
adv.roll = roll;

end