function [Cg] = get_cg(k,depth);
g=9.81;
%get_cg(k,depth);   
%0.5 * omega.* k .* ( 1 + (2*k*depth)./sinh(2*k*depth) ); 
Cp = sqrt(g./k.*tanh(k.*depth));
Cg = 0.5*Cp.*( 1 + (k*depth).*(1-(tanh(k*depth)).^2)./tanh(k*depth) ); 

end

