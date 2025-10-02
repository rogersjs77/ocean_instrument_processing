%function [C,Cg] = get_c_cg(om,h)
% given depth h, wave period T
% returns phase speed C, group velocity Cg
%
% Justin Rogers
% Stanford EFML, 2013

function [C,Cg] = get_c_cg(om,h)

k = get_wavenumber(om,h);
n = 0.5*(1+2*k.*h./sinh(2*k.*h));

C = om./k;
Cg = n.*C;

end