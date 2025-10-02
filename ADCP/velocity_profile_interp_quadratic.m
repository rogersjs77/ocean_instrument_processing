function [vel_extrap] = velocity_profile_interp_quadratic(vel,bin_h,bin_h_extrap,beta)
% this function extrapolates data from bottom to surface
% no slip bottom
% assume top is best estimate for top
% rough quadratic log layer u = ah^2+bh+c;
% no slip bottom, 
% beta is free parameter between 1 and 2. 1=linear, 2=no stress at h_bin(1)
 
vbot = zeros(1,size(vel,2)); 

vtop = mean(vel(end,:),1); 

a = (1-beta)*vel(1,:)./bin_h(1)^2;
b = beta*vel(1,:)/bin_h(1);
c = 0;
indx = bin_h_extrap <  bin_h(1); % find lowest bins;

hll = bin_h_extrap(indx);
for j = 1:length(hll)
vll(j,:) = a.*hll(j).^2+b.*hll(j)+c;
end
tempv = [vll; vel; vtop];
temph = [hll bin_h bin_h_extrap(end)];

vel_extrap = interp1(temph,tempv,bin_h_extrap,'linear','extrap');
end