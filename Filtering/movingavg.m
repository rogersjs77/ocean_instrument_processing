%
%        avgx = movingavg(x,m,lbx,ubx,badx,rmend);
%
%         x : original time series
%         m : number of points of averaging
%
%       ex. avgx = movingavg(x,6,-1.0e+20,1.0e+20,NaN,'y');
%
%        if x is hourly data
%           avgx(1:3) = NaN;
%           avgx(4) = mean(x(1:7)
%           ie, aveaging over 7 hours.
%        if rmend ='y'; the end points will be removed
%
   function avgx = movingavg(x,m,lbx,ubx,badx,rmend);

   x = x(:)';
   if (rem(m,2) == 1);
      nlost = (m-1)/2;
   else
      nlost = m/2;
   end
   begpt = nlost;
   [nrow,ncol]=size(x); 
   endpt = ncol - nlost;
   avgx(:,1:begpt) = NaN*ones(nrow,begpt);
   avgx(:,endpt+1:ncol) = NaN*ones(nrow,nlost);
   for i=begpt+1:endpt;
       bins = (i-nlost):(i+nlost);
       ax = mymean(x(:,bins)',lbx,ubx,badx);
       avgx(:,i) = ax';
   end
   if exist('rmend');
      if (rmend == 'y');
          avgx([1:nlost (ncol-nlost+1:ncol)]) = [];
      end
   end
