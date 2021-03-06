function p = JohnsonCDF(x, gamma, delta, xi, sigma, type)

% save Johnson_tmp
% return

%         gamma - shape parameter #1
%         delta - shape parameter #2
%         xi - location parameter
%         sigma - scale parameter

if isscalar(gamma)
    gamma = repmat(gamma,size(x,1),1);
end
if isvector(gamma) && size(gamma,1) ~= size(x,1)
    error('The length of the vector of the Johnson parameter gamma does not match the number of elements in x')
end
if isscalar(delta)
    delta = repmat(delta,size(x,1),1);
end
if isvector(delta) && size(delta,1) ~= size(x,1)
    error('The length of the vector of the Johnson parameter delta does not match the number of elements in x')
end
if isscalar(xi)
    xi = repmat(xi,size(x,1),1);
end
if isvector(xi) && size(xi,1) ~= size(x,1)
    error('The length of the vector of the Johnson parameter mi does not match the number of elements in x')
end
if isscalar(sigma)
    sigma = repmat(sigma,size(x,1),1);
end
if isvector(sigma) && size(sigma,1) ~= size(x,1)
    error('The length of the vector of the Johnson parameter sigma does not match the number of elements in x')
end

switch type
	case 'SU' % unbounded
        p = 0.5.*(1+erf((gamma+delta.*asinh((x-xi)./sigma))./(2^0.5)));
	case 'SL' % semi-bounded
        p = zeros(size(x));
        ind1 = (x > xi) & (x <= xi+sigma);
        if sum(ind1) > 0
            p(ind1) = 0.5.*erfc(-(gamma(ind1)+delta(ind1).*log((x(ind1)-xi(ind1))./sigma(ind1)))./(2^0.5));
        end
        ind2 = (x > xi + sigma) & (x < Inf);
        p(ind2) = 0.5.*(1+erf((gamma(ind2)+delta(ind2).*log((x(ind2)-xi(ind2))./sigma(ind2)))./(2^0.5)));
        p(x == Inf) = 1; 
	case 'SB' % bounded
        p = zeros(size(x));
        ind1 = (x > xi) & (x < xi+sigma/2);
        p(ind1) = 0.5.*erfc(-(gamma(ind1)+delta(ind1).*log((x(ind1)-xi(ind1))./(-x(ind1)+xi(ind1)+sigma(ind1))))./(2^0.5));
        ind2 = (x >= xi + sigma/2) & (x < xi + sigma);
        p(ind2) = 0.5.*(1+erf((gamma(ind2)+delta(ind2).*log((x(ind2)-xi(ind2))./(-x(ind2)+xi(ind2)+sigma(ind2))))./(2^0.5)));
        p(x >=  xi + sigma) = 1;    
%         p(x <= mi) = 0;
   otherwise
      error('Unknown distribution type. Possible options: SU, SL, SB');
end
