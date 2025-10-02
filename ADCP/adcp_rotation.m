%Ryan Walter, August 2009
%Rotates to principal axes defined by the angle theta

function [vel_maj_paa, vel_min_paa]= adcp_rotation(u, v, theta)
% PRINCAX Principal axis, rotation angle, principal ellipse
%
%   [theta]=princax(east,north) 
%
%   Input:  u = east velocity matrix (depth, time)
%           v = north velocity matrix (depth, time)
%           theta = rotation angle (in degrees) from principal axes (positive defined as
%           CCW)
%   Output: vel_maj_paa = major axis velocity matrix (depth, time)
%            vel_min_paa = minor axis velocity matrix (depth, time)
%    

if length(theta)>1
   [theta,~]=meshgrid(theta,u(:,1)); 
end

w = complex(u,v);
% ind = isfinite(w); % use only the good (finite) points for calculating variance and
% w = w(ind);    %only finite points
wr = w.*exp(-i*theta*pi/180);   %negative indicates a CCW rotation.  multiply by pi/180 for angle in degrees
vel_maj_paa = real(wr);
vel_min_paa = imag(wr);

