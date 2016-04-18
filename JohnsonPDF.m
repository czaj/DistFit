function p = JohnsonPDF(x, gamma, delta, mi, sigma, type)

%         gamma - shape parameter #1
%         delta - shape parameter #2
%         mi - location parameter
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
if isscalar(mi)
    mi = repmat(mi,size(x,1),1);
end
if isvector(mi) && size(mi,1) ~= size(x,1)
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
        p = (exp(-0.5.*(gamma+delta.*asinh((x-mi)./sigma)).^2).*delta)./((2*pi).^(0.5).*((x-mi).^2+sigma.^2).^(0.5));
	case 'SL' % semi-bounded
        p = zeros(size(x));
        p(x > mi) = (exp(-0.5*(gamma(x > mi)+delta(x > mi).*log((x(x > mi)-mi(x > mi))./sigma(x > mi))).^2).*delta(x > mi)) ./ ((2*pi)^(0.5).*(x(x > mi)-mi(x > mi)));
%         p(x <= mi) = 0;
	case 'SB' % bounded
        p = zeros(size(x));
        p((x > mi) & (x < mi + sigma)) = (exp(-0.5.*(gamma((x > mi) & (x < mi + sigma))+delta((x > mi) & (x < mi + sigma)).*log((x((x > mi) & (x < mi + sigma))-mi((x > mi) & (x < mi + sigma)))./(-x((x > mi) & (x < mi + sigma))+mi((x > mi) & (x < mi + sigma))+sigma((x > mi) & (x < mi + sigma))))).^2).*delta((x > mi) & (x < mi + sigma)).*sigma((x > mi) & (x < mi + sigma)))./((2*pi)^(0.5).*(x((x > mi) & (x < mi + sigma))-mi((x > mi) & (x < mi + sigma))).*(-x((x > mi) & (x < mi + sigma))+mi((x > mi) & (x < mi + sigma))+sigma((x > mi) & (x < mi + sigma))));
%         p((x <= mi) | (x >= mi + sigma)) = 0;  
   otherwise
      error('Unknown distribution type. Possible options: SU, SL, SB');
end