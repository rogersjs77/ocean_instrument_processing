function [C] = Correlation(X,Y,corr_method,poly_order)
% takes data x, and y returns srtructure file C with
% correlations, lag, polynomial fit coefficients,
% r2 and rmse errors
%

% limit data to smallest good data
indxx = find(~isnan(X));
indxy = find(~isnan(Y));
if length(indxx)>length(indxy) 
    X = X(indxy);
    Y = Y(indxy);
else
    X = X(indxx);
    Y = Y(indxx);
end

X = naninterp(X);
Y = naninterp(Y);

[C.Corr,C.lag] = xcorr(X,Y,corr_method);
        
C.Poly = polyfit(X,Y,poly_order);

C.fit = polyval(C.Poly,X);
      
[C.r2,C.rmse] = rsquare(Y,C.fit);

C.X=X;
C.Y=Y;

end