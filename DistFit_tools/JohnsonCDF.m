function p = JohnsonCDF(x, gamma, delta, mi, sigma, type)

% save Johnson_tmp
% return

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
        p = 0.5.*(1+erf((gamma+delta.*asinh((x-mi)./sigma))./(2^0.5)));
	case 'SL' % semi-bounded
        if mi == 0
            cprintf(rgb('DarkOrange'), 'Semi-bounded SL at 0')
        end
        p = zeros(size(x));
        p((x > mi) & (x <= mi+sigma)) = 0.5.*erfc(-(gamma+delta.*log((x((x > mi) & (x <= mi+sigma))-mi)./sigma))./(2^0.5));
        p(x > mi + sigma) = 0.5.*(1+erf((gamma+delta.*log((x(x > mi + sigma)-mi)./sigma))./(2^0.5)));
        p(x <= mi) = 0; % Co dla x <=  mi? Czy wg Wolframa to 0? Bo to wg mnie jest
        % niejasne, tam jest s³owo "true" zamiast przedzia³u.
	case 'SB' % bounded
        if mi == 0
            cprintf(rgb('DarkOrange'), 'Semi-bounded SB at 0')
        end
        p = zeros(size(x));
        p((x > mi) & (x < mi+sigma/2)) = 0.5.*erfc(-(gamma+delta.*log((x((x > mi) & (x < mi+sigma/2))-mi)./(-x((x > mi) & (x < mi+sigma/2))+mi+sigma)))./(2^0.5));
        p((x >= mi + sigma/2) & (x < mi + sigma)) = 0.5.*(1+erf((gamma+delta.*log((x((x >= mi + sigma/2) & (x < mi + sigma))-mi)./(-x((x >= mi + sigma/2) & (x < mi + sigma))+mi+sigma)))./(2^0.5)));
        p(x >=  mi + sigma) = 1;    
        p(x <= mi) = 0; % Zgadza siê?
   otherwise
      error('Unknown distribution type');
end
