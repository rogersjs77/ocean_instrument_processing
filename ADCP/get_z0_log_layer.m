function [roughness] = get_z0_log_layer(vel_x,vel_y,bin_MAB,depth,max_depth_fract,vel_min_corr,zref)
% This function finds roughnes height z0 from filtered velocities,
% assuming a log layer fit to the data
%  u(z,t) = ustar(t)/k*ln((z-d)/z0), Reidenbach, 2006
% limit fit to <U> > vel_min_corr, then performs least squares fit to parameters
% only uses results with good correlations: r2> r2min, 
%   and realistic roughness: 30*z0 + d < max_depth_fract*Depth/2
% mean values computed as <z0> = exp(<ln(z0(t))>)
%                         <d> = exp(<ln(d(t)>)
% 
% Justin Rogers, Stanford EFML, 5/2014
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 0.41; % von karman const
r2min = 0.80; % per reidenbach et al, gives 95% confidence
% zref = 1; % reference height above bottom for Cd calc.
z0min = 0.01/30; % set depending on conditions; z0*30 ~ physical roughness scale
dmin = 0.01; % set depending on conditions;


%% calculations


%% compute over time
for i=1:size(vel_x,2)
disp(['processing ' num2str(i) ' of ' num2str(size(vel_x,2))]);    
 
h = depth(i);
roughness.h(i) = h;

n = bin_MAB < max_depth_fract*depth(i);
ft = fittype( 'real((a/0.41)*(log((x-c)/b) + 2*0.2*(sin(pi*(x-c)/(2*h)).^2)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 1E-4 0.01 h]; % ustar can be any value, z0 must be +, d must be 0
opts.StartPoint = [0.1 0.01 0.01 h];
opts.Upper = [Inf max_depth_fract*depth(i) 0.01 h]; 

% opts.Lower = [-Inf 1E-4 1E-4 h]; % ustar can be any value, z0 must be +, d must be +
% opts.StartPoint = [0.1 0.01 0.01 h];
% ustar can be any value, % z0 must be smaller than depth, % d must be less than depth
% opts.Upper = [Inf max_depth_fract*depth(i) max_depth_fract*depth(i) h]; 


xin =  bin_MAB(n);
% if length(xin)<3
%     xin = bin_MAB;
%     disp('not enough points in z, (<3), exiting')
% %     return
% end

    % check if real data and bin_MAB>3
    if  length(xin)<4 || isnan(h) || min(isfinite(vel_x(n,i)))==0 || min(isfinite(vel_y(n,i)))==0 % any vel data is NaN
        roughness.ustar(i) = nan;
        roughness.z0.ut(i) = nan;
        roughness.d.ut(i) = nan;
        roughness.r2u(i) = nan;
        
        roughness.vstar(i) = nan;
        roughness.z0.vt(i) = nan;
        roughness.d.vt(i) = nan;
        roughness.r2v(i) = nan;
    else

    % Fit model to data.
    
    % major axis
    if abs(mean(vel_x(n,i),1)) < vel_min_corr
    roughness.ustar(i) = nan;
    roughness.z0.ut(i) = nan;
    roughness.d.ut(i) = nan;
    roughness.r2u(i) = nan;
    else
    yin = vel_x(n,i)';
    [xData, yData] = prepareCurveData(xin,yin);
    [fitresult, gof] = fit( xData, yData, ft, opts );
    coeffs = coeffvalues(fitresult);
    roughness.ustar(i) = coeffs(1);
    roughness.z0.ut(i) = coeffs(2);
    roughness.d.ut(i) = coeffs(3);
    roughness.r2u(i) = gof.rsquare;
    end

    % minor axis
    if abs(mean(vel_y(n,i),1)) < vel_min_corr
    roughness.vstar(i) = nan;
    roughness.z0.vt(i) = nan;
    roughness.d.vt(i) = nan;
    roughness.r2v(i) = nan;
    else
    yin = vel_y(n,i)';
    [xData, yData] = prepareCurveData(xin,yin);
    [fitresult, gof] = fit( xData, yData, ft, opts );
    coeffs = coeffvalues(fitresult);
    roughness.vstar(i) = coeffs(1);
    roughness.z0.vt(i) = coeffs(2);
    roughness.d.vt(i) = coeffs(3);
    roughness.r2v(i) = gof.rsquare; 
    end
    end
end

% remove low correlation values
indx = roughness.r2u < r2min;
roughness.z0.ut(indx) = NaN;
roughness.d.ut(indx) = NaN;
indx = roughness.r2v < r2min;
roughness.z0.vt(indx) = NaN;
roughness.d.vt(indx) = NaN;

% remove unrealistic results 30*z0+d > 2*depth
% Depth = nanmean(depth);
indx = 30*roughness.z0.ut+roughness.d.ut > 4*roughness.h; 
roughness.z0.ut(indx) = NaN;
roughness.d.ut(indx) = NaN;
indx = 30*roughness.z0.vt+roughness.d.vt > 4*roughness.h; 
roughness.z0.vt(indx) = NaN;
roughness.d.vt(indx) = NaN;
% remove very low z0 < z0min 
indx = roughness.z0.ut < z0min;
roughness.z0.ut(indx) = NaN;
roughness.d.ut(indx) = NaN;
indx = roughness.z0.vt < z0min;
roughness.z0.vt(indx) = NaN;
roughness.d.vt(indx) = NaN;
% remove very low d < dmin
indx = roughness.d.ut < dmin;
roughness.z0.ut(indx) = NaN;
roughness.d.ut(indx) = NaN;
indx = roughness.d.vt < dmin;
roughness.z0.vt(indx) = NaN;
roughness.d.vt(indx) = NaN;

% compute drag coeff Cd for depth avg flow
PI=0.2;
roughness.Cd.ut = k^2*(log(roughness.h./roughness.z0.ut) + (PI-1)).^-2;
roughness.Cd.vt = k^2*(log(roughness.h./roughness.z0.vt) + (PI-1)).^-2;

% compute drag coeff Cd at z reference
roughness.Cd_zref.ut = real(k^2./(log((zref-roughness.d.ut)./roughness.z0.ut)).^2);
roughness.Cd_zref.vt = real(k^2./(log((zref-roughness.d.vt)./roughness.z0.vt)).^2);

% compute statistics for z0, mean is log mean
[roughness.z0.u] = logmean(roughness.z0.ut);
[roughness.z0.v] = logmean(roughness.z0.vt);

% compute statistics for d, mean is log mean
[roughness.d.u] = logmean(roughness.d.ut);
[roughness.d.v] = logmean(roughness.d.vt);

% compute statistics for Cd depth avg
[roughness.Cd.u] = logmean(roughness.Cd.ut);
[roughness.Cd.v] = logmean(roughness.Cd.vt);

% compute statistics for Cd at z reference
[roughness.Cd_zref.u] = logmean(roughness.Cd_zref.ut);
[roughness.Cd_zref.v] = logmean(roughness.Cd_zref.vt);

end


