function [lm] = logmean(x)
% this function computes the logmean of a vector x
% xbar = exp(<log(x)>)
% and associated error in x, xstd
% returns structure file with
% lm.xbar, log mean
% lm.xstd, log std
% lm.percent_nan, percent of nans in set
% lm.good_samples, number of good samples used for averaging

% take log
f1 = real(log(x));
s1 = abs(nanstd(x)./x);

% take mean & std
f2 = nanmean(f1);
s2 = nanstd(f1);
% s2 = (1/sum(~isnan(f1))).*sqrt(nansum(s1.^2));

% take exponential
lm.xbar = exp(f2);
lm.xstd = abs(lm.xbar*s2);

lm.percent_nan = sum(isnan(x))/length(x)*100;
lm.good_samples = sum(~isnan(x));

end