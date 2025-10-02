
%   function [nu,nv]=rotateuvfun(u,v,theta);
%
%    theta in degrees, positive for counterclockwise 
%                      negative for clockwise rotation
%
    function [nu,nv]=rotateuvfun(u,v,theta);
    theta = theta*pi/180;
    nu = u.*cos(theta)+v.*sin(theta);
    nv = -u.*sin(theta)+v.*cos(theta);
