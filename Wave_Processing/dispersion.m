function [k,om] = dispersion(h,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear dispersion relation solver [om^2=gk tanh(kh)]
% Justin Rogers, Stanford EFML, 2012
% this is an adaptive code which minimizes solving time for k
% by using approximate solution for k as initial guess
%
% h is depth [m], can be matrix of any size
% T is scalar wave period [s]
% k is wavenumber [1/m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g=9.81;
k=0.*h;
om=2*pi./T;
test=0.*h; % tracks method used for calculation
for i = 1:size(h,1)
    for j=1:size(h,2)
        k0 = dispersion_approx(h(i,j),T); % use approximation for starting
        if k0*h(i,j)>10 % deep water
            k(i,j) = om^2/g;
            test(i,j)=2;
        elseif k0*h(i,j)<0.04 % shallow water
            k(i,j)=om./sqrt(g*h(i,j));
            test(i,j)=4;
        else  % mid-range use full dispersion         
            f=@(k)om.^2-g*k*tanh(k*h(i,j));
            k(i,j)=fzero(f,k0); % start with previous k0
            test(i,j)=3;
        end
    end
end