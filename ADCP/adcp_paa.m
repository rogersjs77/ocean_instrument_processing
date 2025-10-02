%Ryan Walter, August 2009

%Computes the along channel (u) and cross-channel (v) components using
%principle axes analysis using only points that are finite (not NaN) by
%calculating the variance ellipse and using max/min variance via analysis
%found in Emery and Thompson S4.3.5

function [theta]= adcp_paa(u,v)
% PRINCAX Principal axis, rotation angle, principal ellipse
%
%   [theta]=princax(east,north) 
%
%   Input:  u = east velocity vector (depth-averaged typically)
%           v = north velocity vector (depth-averaged typically)
%   Output: theta = angle of maximum variance, math notation (east == 0, north=90)       
%
% For derivation, see Emery and Thompson, "Data Analysis Methods
%   in Oceanography", 1998, Pergamon, pages 325-327.  ISBN 0 08 0314341
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 1.0 (12/4/1996) Rich Signell (rsignell@usgs.gov)
% Version 1.1 (4/21/1999) Rich Signell (rsignell@usgs.gov) 
%     fixed bug that sometimes caused the imaginary part 
%     of the rotated time series to be aligned with major axis. 
%     Also simplified the code.
% Version 1.2 (3/1/2000) Rich Signell (rsignell@usgs.gov) 
%     Simplified maj and min axis computations and added reference
%     to Emery and Thompson book
% Updated August 2009 by Ryan Walter 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = complex(u,v);   %complex variable w = u + iv
ind = isfinite(w);  % use only the good (finite) points for calculating variance and
w = w(ind);         %only finite points

% find covariance matrix
cv = cov([real(w(:)) imag(w(:))]); 

% find direction of maximum variance (degrees)
theta = 180/pi*(0.5*atan2(2.*cv(2,1),(cv(1,1)-cv(2,2)) ));