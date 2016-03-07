function p = JohnsonCDF(x, gamma, delta, mi, sigma, type)

% save Johnson_tmp
% return

%         gamma - shape parameter #1
%         delta - shape parameter #2
%         mi - location parameter
%         sigma - scale parameter

if isscalar(gamma)
    gamma = repmat(gamma,size(x,1),1);
end
if isvector(gamma) && size(gamma,1) ~= size(x,1)
    error('The lenght of the vector of Johnson parameter gamma does not match the number of elements in x')
end
if isscalar(delta)
    delta = repmat(delta,size(x,1),1);
end
if isvector(delta) && size(delta,1) ~= size(x,1)
    error('The lenght of the vector of Johnson parameter delta does not match the number of elements in x')
end
if isscalar(mi)
    mi = repmat(mi,size(x,1),1);
end
if isvector(mi) && size(mi,1) ~= size(x,1)
    error('The lenght of the vector of Johnson parameter mi does not match the number of elements in x')
end
if isscalar(sigma)
    sigma = repmat(sigma,size(x,1),1);
end
if isvector(sigma) && size(sigma,1) ~= size(x,1)
    error('The lenght of the vector of Johnson parameter sigma does not match the number of elements in x')
end

switch type
	case 'SU' % unbounded
        p = 0.5.*(1+erf((gamma+delta.*asinh((x-mi)./sigma))./(2^0.5)));
	case 'SL' % semi-bounded
        p = zeros(size(x));
        ind1 = (x > mi) & (x <= mi+sigma);
        p(ind1) = 0.5.*erfc(-(gamma(ind1)+delta(ind1).*log((x(ind1)-mi(ind1))./sigma(ind1)))./(2^0.5));
        ind2 = (x > mi + sigma) & (x < Inf);
        p(ind2) = 0.5.*(1+erf((gamma(ind2)+delta(ind2).*log((x(ind2)-mi(ind2))./sigma(ind2)))./(2^0.5)));
        p(x == Inf) = 1; 
	case 'SB' % bounded
        p = zeros(size(x));
        p((x > mi) & (x < mi+sigma/2)) = 0.5.*erfc(-(gamma((x > mi) & (x < mi+sigma/2))+delta((x > mi) & (x < mi+sigma/2)).*log((x((x > mi) & (x < mi+sigma/2))-mi((x > mi) & (x < mi+sigma/2)))./(-x((x > mi) & (x < mi+sigma/2))+mi((x > mi) & (x < mi+sigma/2))+sigma((x > mi) & (x < mi+sigma/2)))))./(2^0.5));
        p((x >= mi + sigma/2) & (x < mi + sigma)) = 0.5.*(1+erf((gamma((x >= mi + sigma/2) & (x < mi + sigma))+delta((x >= mi + sigma/2) & (x < mi + sigma)).*log((x((x >= mi + sigma/2) & (x < mi + sigma))-mi((x >= mi + sigma/2) & (x < mi + sigma)))./(-x((x >= mi + sigma/2) & (x < mi + sigma))+mi((x >= mi + sigma/2) & (x < mi + sigma))+sigma((x >= mi + sigma/2) & (x < mi + sigma)))))./(2^0.5)));
        p(x >=  mi + sigma) = 1;    
%         p(x <= mi) = 0;
   otherwise
      error('Unknown distribution type');
end
