%------------------------------------------------------------------------------
%       PGM : LPASS                       Ren-Chieh Lien
%
%             Low pass time series using Butterworth filter
%          
%               xlpass = lpass(xoriginal, tinterval, cutofft, npole,[XO])
% 
%             EX : xlpass = lpass(x0, 25, 120, 2)
%                           
%                  tinterval = 25 units (eg. secs, mins, hours, ...etc.)
%                  
%                  cutoff time = 120 units
%
%                  npole = 2 (poles used in Butterworth filter, default 2)
%
%                  [XO] : optional, 'Y' ===> plot original time series
%------------------------------------------------------------------------------
%
%
        function [yl,ff,hl]=lpass(x,delt,cutt,npole,XOPLOT);
        cutt = cutt/2;
%
%             CHECK THE INPUT ARGUMENT
%
if any(isnan(x));
   !echo time series has NaN or Inf. Moving averages will be proceeded.
   yl = movingavg(x,floor(cutt/delt),-Inf,Inf,NaN);
   ff = []; hl = [];
else
        if ~exist('XOPLOT');
           XOPLOT = 'N';
        end
	nptslost = ceil(cutt/delt);
%
%             REMOVE THE MEAN FROM THE TIME SERIES FIRST
%
    orignalx = x;
	bad = find(isnan(x));
	good = find(~isnan(x));
	npts = length(x);
	if (length(bad) == 1 & (bad == 1 | bad == length(x)));
	   x(bad )= []; add = 'y';
        end
	if (length(bad) > 1);
	   if (min(bad) == 1 & sum(diff(diff(bad)).^2) == 0);
	      x(bad) = []; add = 'y';
           end
	   if (max(bad) == npts & sum(diff(diff(bad)).^2) == 0);
	      x(bad) = []; add = 'y';
           end
	end  

        xo = x;
        x0=mean(x);
        x = x - x0;
%
%             SET UP THE SAMPLING AND CUTOFF FREQUENCY
%
        fs = 1/delt;                                   % sampling frequency
        fcut = 1/cutt;                                 % cutoff frequency
        lband = fcut/fs;                               % normalized lpass band
%
%             CONSTRUCT THE BUTTERWORTH FILTER
%
        time = (0:(length(x)-1))*delt;
%    npole
%    lband
        [cl,ch] = butter(npole,lband);            
%       [cl,ch] = cheby1(npole,4,lband);            
%
%             CHECK THE FREQUENCY RESPONSE FUNCTION
%
        nfreq = 128*2;
        ff = fs*(0:nfreq-1)/(2*nfreq);
        hl = freqz(cl,ch,nfreq);
        hl = hl.*conj(hl);
%       subplot(221)
%       semilogx(ff(2:nfreq),hl(2:nfreq))
        i=1:length(x);
%
%             APPLYING THE FILTER FIRST
%
        y = filter(cl,ch,x);
%
%             REVERSING THE TIME SERIES AND APPLYING THE FILTER AGAIN
%             TO REMOVE THE PHASE SHIFT DUE TO THE FILTERING
%
        z0 = y(length(y):-1:1);
        z = filter(cl,ch,z0);
%
%             REVERSE THE TIME SERIES BACK AGAIN AND ADD THE MEAN
% 
        w = z(length(z):-1:1);
	    w1 = filtfilt(cl,ch,x);
        xl =w + x0;
    	yl = w1 + x0;

        npts = length(xl);

        add = 'y';
        if strcmp(add,'y');
     	  x(good) = xl;
          xl = x;
        end
   end
