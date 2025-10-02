%------------------------------------------------------------------------------
%    mymean                   06-11-92                   Ren-Chieh Lien
%
%    mean after excluding NaN or Inf.
%    function mx = mymean(x,lcut,hcut,badflag);
%
%------------------------------------------------------------------------------
    function mx = mymean(x,lcut,hcut,badflag,ncrt);
%   if size(x,1) == 1; x = x(:); end;
    if (nargin <= 4); ncrt = 1; end                  % set default
    if (nargin <=3); badflag = NaN; end              % set default
    if (nargin <= 2); hcut = Inf; end                % set default
    if (nargin <= 1); lcut = -Inf; end               % set default
     
    [row col]=size(x); 
    mx = NaN*ones(1,col);
%   if row==1; x = x(:); end   % convert row to column vector
    k = 0;
    for i = 1:col;
        gxind = ...
	  find( ~isnan(x(:,i)) & ~isinf(x(:,i)) & ...
	  x(:,i) <= hcut & x(:,i) >= lcut & x(:,i) ~= badflag);
	  if (length(gxind) >= ncrt);
	     mx(i) =mean(x(gxind,i));
%            indx(i) = i;
          else
	     mx(i) = NaN;
	     k = k+1;
%            noindx(k) = i;
          end
       clear gxind
    end
