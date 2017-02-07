function p = JohnsonPDF(x, gamma, delta, xi, sigma, type)

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
        p = (exp(-0.5.*(gamma+delta.*asinh((x-xi)./sigma)).^2).*delta)./((2*pi).^(0.5).*((x-xi).^2+sigma.^2).^(0.5));
	case 'SL' % semi-bounded
        p = zeros(size(x));
        p(x > xi) = (exp(-0.5*(gamma(x > xi)+delta(x > xi).*log((x(x > xi)-xi(x > xi))./sigma(x > xi))).^2).*delta(x > xi)) ./ ((2*pi)^(0.5).*(x(x > xi)-xi(x > xi)));
%         p(x <= mi) = 0;
	case 'SB' % bounded
        p = zeros(size(x));
        p((x > xi) & (x < xi + sigma)) = (exp(-0.5.*(gamma((x > xi) & (x < xi + sigma))+delta((x > xi) & (x < xi + sigma)).*log((x((x > xi) & (x < xi + sigma))-xi((x > xi) & (x < xi + sigma)))./(-x((x > xi) & (x < xi + sigma))+xi((x > xi) & (x < xi + sigma))+sigma((x > xi) & (x < xi + sigma))))).^2).*delta((x > xi) & (x < xi + sigma)).*sigma((x > xi) & (x < xi + sigma)))./((2*pi)^(0.5).*(x((x > xi) & (x < xi + sigma))-xi((x > xi) & (x < xi + sigma))).*(-x((x > xi) & (x < xi + sigma))+xi((x > xi) & (x < xi + sigma))+sigma((x > xi) & (x < xi + sigma))));
%         p((x <= mi) | (x >= mi + sigma)) = 0;  
   otherwise
      error('Unknown distribution type. Possible options: SU, SL, SB');
end