% function [adcp] = adcp_Beam_to_XYZ(adcp)
% 
% Convert ADCP velocities from beam to xyz coordinates
% Justin Rogers & Kara Scheu
% Stanford EFML, 3/2013
%
% INPUT
% adcp = standard adcp structure file from rdradcp in beam coordinates
% 
% OUTPUT
% adcp = adcp structure file with velocities in xyz coordinates
%        velocities (east_vel, north_vel, vert_vel, error_vel)
%
%%
function [adcp] = adcp_Beam_to_XYZ(adcp)

% find coefficients
% Where: c = +1 for a convex transducer head, -1 for concave
% a = 1/ [2 sin(?)] = 1.4619 for ? = 20°, 1.0000 for ? = 30°
% b = 1/ [4 cos(? )] = 0.2660 for ? = 20°, 0.2887 for ? = 30°
% d = a/ 2^0.5 = 1.0337 for ? = 20°, 0.7071 for ? = 30°

a = 1 / (2* sin(deg2rad(adcp.config(1,1).beam_angle)));
b = 1 / (4* cos(deg2rad(adcp.config(1,1).beam_angle)));
if strcmp(adcp.config(1,1).beam_pattern,'convex')
    c = 1;
elseif strcmp(adcp.config(1,1).beam_pattern,'concave')
    c = -1;
end
d = a /sqrt(2);

N = [c*a -c*a 0 0; 0 0 -c*a c*a; b b b b; d d -d -d];
Ninv = inv(N);

% compute transformed matrices and assign values
for i = 1:size(adcp.east_vel,1)

    Vemp = N*[adcp.east_vel(i,:);adcp.north_vel(i,:);...
        adcp.vert_vel(i,:); adcp.error_vel(i,:)];
    adcp.east_vel(i,:) = Vemp(1,:);
    adcp.north_vel(i,:) = Vemp(2,:);
    adcp.vert_vel(i,:) = Vemp(3,:);
    adcp.error_vel(i,:) = Vemp(4,:);
    
end

adcp.config(1,1).coord_sys = 'xyz';
adcp.coord_sys = 'xyz';
adcp.coord_expl = 'east_vel=x, north_vel=y, vert_vel=z, error_vel = errorz';

end