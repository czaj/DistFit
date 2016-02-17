function p = Johnson_PDF(x, gamma, delta, mi, sigma, type)

%         gamma - shape parameter #1
%         delta - shape parameter #2
%         mi - location parameter
%         sigma - scale parameter

if isscalar(gamma)
    gamma = repmat(gamma,size(x,1),1)';
end
if isvector(gamma) && size(gamma,1) ~= size(x,1)
    error('The lenght of the vector of Johnson parameter gamma does not match the number of elements in x')
end
if isscalar(delta)
    delta = repmat(delta,size(x,1),1)';
end
if isvector(delta) && size(delta,1) ~= size(x,1)
    error('The lenght of the vector of Johnson parameter delta does not match the number of elements in x')
end
if isscalar(mi)
    mi = repmat(mi,size(x,1),1)';
end
if isvector(mi) && size(mi,1) ~= size(x,1)
    error('The lenght of the vector of Johnson parameter mi does not match the number of elements in x')
end
if isscalar(sigma)
    sigma = repmat(sigma,size(x,1),1)';
end
if isvector(sigma) && size(sigma,1) ~= size(x,1)
    error('The lenght of the vector of Johnson parameter sigma does not match the number of elements in x')
end

switch type
	case 'SU' % unbounded
        p = (exp(-0.5.*(gamma+delta.*arcsinh((x-mi)./sigma))^2).*delta)./((2*pi)^(0.5).*((x-mi)^2+sigma^2)^(0.5));
	case 'SL' % semi-bounded
        if mi == 0
            cprintf(rgb('DarkOrange'), 'Semi-bounded SL at 0')
        end
        p = zeros(size(x));
        p(x > mi) = (exp(-0.5.*(gamma+delta.*arcsinh((x(x > mi)-mi)./sigma))^2)*delta)./((2*pi)^(0.5).*(x(x > mi)-mi));
        p(x <= mi) = 0; % Czy to si� zgadza? 
	case 'SB' % bounded
        if mi == 0
            cprintf(rgb('DarkOrange'), 'Semi-bounded SB at 0')
        end
        p = zeros(size(x));
        p((x > mi) & (x < mi + sigma)) = (exp(-0.5.*(gamma+delta.*log((x((x > mi) & (x < mi + sigma))-mi)./(-x((x > mi) & (x < mi + sigma))+mi+sigma)))^2).*delta.*sigma)./((2*pi)^(0.5).*(x((x > mi) & (x < mi + sigma))-mi).*(-x((x > mi) & (x < mi + sigma))+mi+sigma));
        p((x <= mi) | (x >= mi + sigma)) = 0; % Czy to si� zgadza? 
   otherwise
      error('Unknown distribution type');
end