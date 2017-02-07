function r = JohnsonRND(type, gamma, delta, mi, sigma, out_size)

% http://www.dtic.mil/dtic/tr/fulltext/u2/762722.pdf

r = random('normal',0,1,out_size);

switch type
	case 'SU' % unbounded      
        r = mi + sigma * sinh((r - gamma)./delta);
 	case 'SL' % semi-bounded - untested
        r = mi + exp((r - gamma)./delta);
	case 'SB' % bounded - untested
        r_tmp = exp((r - gamma)./delta);
        r = mi + sigma * r_tmp ./ (1 + r_tmp);
   otherwise
      error('Unknown distribution type. Possible options: SU, SL, SB');
end
