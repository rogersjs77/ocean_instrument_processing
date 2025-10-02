% function [adcp] = adcp_three_beam_solution(adcp,bad_beam,time_id)
%
% Three Beam Solution to ADCP Data
% developed by Justin Rogers & Kara Scheu
% Stanford EFML, 3/2013
%
% INPUT
% adcp = adcp structure file from rdradcp in earth coordinates
% bad_beam = beam to remove (1,2,3,or 4)
% time_id = optional index of mtime to use for heading pitch and roll averaqes
%           if empty [], uses full time series for averages 
%
% OUTPUT
% adcp = corrected structure file with east_vel, north_vel and vert_vel
%        from 3 beam solution.  error_vel is now meaningless.
%        error_vel_orig is original error estimate
%%
function [adcp] = adcp_three_beam_solution(adcp,bad_beam,time_id)

adcp.error_vel_orig = adcp.error_vel;

% convert ENU to XYZ coord
[adcp] = adcp_ENU_to_XYZ(adcp,time_id);
% convert XYZ to beam coord
[adcp] = adcp_XYZ_to_Beam(adcp);

% set beam velocities
b1 = adcp.east_vel;
b2 = adcp.north_vel;
b3 = adcp.vert_vel;
b4 = adcp.error_vel;
% set error = 0 = b1+b2-b3-b4 , replace bad beam data
if bad_beam==1
    adcp.bad_beam = b1;
    adcp.east_vel = -b2+b3+b4;
elseif bad_beam==2
    adcp.bad_beam = b2;
    adcp.north_vel = -b1+b3+b4;
elseif bad_beam==3
    adcp.bad_beam = b3;    
    adcp.vert_vel = b1+b2-b4;
elseif bad_beam==4
    adcp.bad_beam = b4;    
    adcp.error_vel = b1+b2-b3;
end

% convert beam to XYZ coord
[adcp] = adcp_Beam_to_XYZ(adcp);

% convert XYZ to ENU coord
[adcp] = adcp_XYZ_to_ENU(adcp,time_id);

end