  function pubfiglabel(x0,y0);

  figlabel=str2mat('(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)',...
                   '(j)','(k)','(l)','(m)','(n)','(o)','(p)','(q)','(r)',...
                   '(s)','(t)','(u)','(v)','(w)','(x)','(y)','(z)');
  h_allaxes = getallaxes;
  k = 0;
  for i = 1:length(h_allaxes);
      tag = get(h_allaxes(i),'Tag');
      if ~strcmp(tag,'Colorbar') & ~strcmp(tag,'legend');
         k = k+1;
         h_axes(k) = h_allaxes(i);
      end
  end
  h_allaxes = h_axes;
  n = length(h_allaxes); 
  for i = 1:length(h_allaxes);
      if nargin == 0; 
         sgtext(figlabel(i,:)); 
      end
      if nargin == 2;
         axes(h_allaxes(n-i+1))
         text(x0,y0,figlabel(i,:), 'unit','normalized');
      end 
  end
