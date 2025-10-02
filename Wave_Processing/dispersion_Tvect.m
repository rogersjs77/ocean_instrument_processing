function [k,om] = dispersion_Tvect(h,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear dispersion relation solver [om^2=gk tanh(kh)]
% Justin Rogers, Stanford EFML, 2012
% this is an adaptive code which minimizes solving time for k
% by using shallow and deep water approximations and k from previous step
%
% h is depth [m], scalar
% direction is columns
% T is vector of wave periods [s]
% k is wavenumber [1/m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g=9.81;
k=0.*h;
om=2*pi./T;
kd = om.^2/g; % deep water approx
test=0.*T; % tracks method used for calculation
for i = 1:size(T,1)
    for j=1:size(T,2)
        if j==1 % calculate first point explicitly
           ks=om(i,j)./sqrt(g*h);% shallow water approx starting point
           k0 = mean([ks kd(i,j)]); % start at midpoint
           f=@(k)om(i,j).^2-g*k*tanh(k*h);
           k(i,j)=fzero(f,k0); % start with intial estimate
           k0=k(i,j);  % update next starting point
           test(i,j)=1;
        elseif k0*h>3.2 % deep water
            k(i,j) = om(i,j)^2/g;
            k0=k(i,j);  % update next starting point
            test(i,j)=2;
        elseif k0*h<0.08 % shallow water
            k(i,j)=om(i,j)./sqrt(g*h);
            k0=k(i,j);  % update next starting point
            test(i,j)=4;
        else  % mid-range use full dispersion         
            f=@(k)om(i,j).^2-g*k*tanh(k*h);
            k(i,j)=fzero(f,k0); % start with previous k0
            k0=k(i,j);  % update next starting point
            test(i,j)=3;
        end
    end
end