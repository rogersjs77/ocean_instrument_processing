%
%    function lpu = mylpass(u,dt,lpdt,npole);
%

    function lpu = mylpass(u,dt,lpdt,npole);
    [nj,nt] = size(u);
    lpu = NaN*ones(size(u));
    for j = 1:nj;
        g = find(~isnan(u(j,:)));
        if length(g) >= 10;
           tmpu = u(j,g);
           intpu = interp1(g,tmpu,1:nt);
           g = find(~isnan(intpu));
           tmp = lpass(intpu(g),dt,lpdt,npole,'b');
           lpu(j,g) = tmp(:)';
        end
    end
