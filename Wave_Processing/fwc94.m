function fwc = fwc94( cmu, cukw )
% FWC94 - Wave-current friction factor
% Equations 32 and 33 in Madsen, 1994
fwc = .00999; %meaningless (small) return value
if( cukw <= 0. ),
fprintf(1,'ERROR: cukw too small in fwc94: %9.4f/n',cukw)
return
end
if( cukw < 0.2 ),
fwc = exp( 7.02*0.2^(-0.078) - 8.82 );
fprintf(1,'WARNING: cukw very small in fwc94: %9.4f/n',cukw)
end
if( (cukw >= 0.2) && (cukw <= 100.) ),
fwc = cmu*exp( 7.02*cukw^(-0.078)-8.82 );
elseif( (cukw > 100.) && (cukw <= 10000.) ),
fwc = cmu*exp( 5.61*cukw^(-0.109)-7.30 );
elseif( cukw > 10000.),
fwc = cmu*exp( 5.61*10000.^(-0.109)-7.30 );
else
fprintf(1,'WARNING: cukw very large in fwc94: % 9.4f/n',cukw)
end