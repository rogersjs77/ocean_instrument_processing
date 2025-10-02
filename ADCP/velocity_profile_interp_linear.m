function [vel_extrap] = velocity_profile_interp_linear(vel,bin_MAB,bin_MAB_extrap)
% this function extrapolates data from bottom to surface
% assume top data point is best estimate for top
% log, law of the wall parameterization for bottom 1/2 of column
% u(z) = ustar/k*ln((z+h)/z0)
% where ustar is friction velocity, k is vonkarman constant 0.41
% z is vertical from surface, h is total depth, z0 is roughmness height

if  min(isfinite(vel))==0 % any vel data is NaN
    vel_extrap = bin_MAB_extrap+NaN;
    return
end


vbot = zeros(1,size(vel,2)); 

vtop = mean(vel(end,:),1); 

% [~,BIN_MAB]=meshgrid(vel(1,:),bin_MAB);
% [~,BIN_MAB_EXTRAP]=meshgrid(vel(1,:),bin_MAB_extrap);

% assume log layer applies to lower 1/2 of column
% [xData, yData] = prepareCurveData( bin_MAB(1:end/2), vel(1:end/2)' );
% % Set up fittype and options.
% ft = fittype( 'real((a/0.41)*log((x)/b))', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( ft );
% opts.Display = 'Off';
% opts.Lower = [-Inf 0.001]; % ustar can be any value, z0 must be +
% opts.StartPoint = [0.01 0.2];
% opts.Upper = [Inf 0.80]; % ustar can be any value, 
% % z0 realistically must be smaller than zref=1m
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% coeffs = coeffvalues(fitresult);
% ustar = coeffs(1);
% z0 = coeffs(2);
% rmse = gof.rmse;
% r2 = gof.rsquare;

% if r2 < 0.8 % if log fit is bad, force linear fit to bottom
% clear fitresult gof xData yData ft opts    
[xData, yData] = prepareCurveData( [0 bin_MAB(1:2)], [0 vel(1:2)'] );
% Set up fittype and options.
ft = fittype( 'a*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf];
opts.StartPoint = [1];
opts.Upper = [Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
r2 = gof.rsquare;    
    
    
% 
% end

indx = bin_MAB_extrap <  bin_MAB(1); % find lowest bins;
hll = bin_MAB_extrap(indx);
vll = feval(fitresult,hll); % evaluated log layer section

% indx = hll <=z0; % set v below z0 to zero 
% vll(indx)=0;
indx = ~isfinite(vll);
vll(indx)=0; % set any inf values to 0

tempv = [vll; vel; vtop];
temph = [hll bin_MAB bin_MAB_extrap(end)];

vel_extrap = interp1(temph,tempv,bin_MAB_extrap,'linear','extrap');
end