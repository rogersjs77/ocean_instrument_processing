function k = wavenum(om,h);
% function that takes the radian wave frequency and
% a vector of depths and returns the wavenumber at that
% depth by solving the dispersion relationship
% using newtons method

% the initial guess will be the shallow water wavenumber
% frequency runs down (om needs to be column) and cross-shore depth runs
% across (h needs to be a row)
%%
som=size(om);
sh=size(h);

om=repmat(om,[1 sh(2)]);
h=repmat(h,[som(1) 1]);

g = 9.81;

k = om./sqrt(g*h);


f = g*k.*tanh(k.*h) - om.^2;
%%
for a=1:10 
  dfdk = g*k.*h.*(sech(k.*h)).^2 + g*tanh(k.*h);
  k = k - f./dfdk;
  f = g*k.*tanh(k.*h) - om.^2;
end
