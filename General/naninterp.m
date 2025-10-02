function X = naninterp(X,fill)
% Interpolate over NaNs
% See INTERP1 for more info
% interpolate all interior NaNs

if max(max(isnan(X)))==1
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)),'linear',NaN);
end
% interpolate exterior nans with nearest neighbor
if max(max(isnan(X)))==1
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)),'nearest','extrap');
end

% if all else fails, enter fill value
if nargin==2
   X(isnan(X))=fill; 
end

return