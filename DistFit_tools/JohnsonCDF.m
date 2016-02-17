function p = JohnsonCDF(x,B,type)

% save Johnson_tmp
% return

B = B(:);

switch type
	case 'SU' % unbounded
        if size(B,1)~=4
            error('Incorrect number of parameters - Johnson SU requires 4 parameters')
        end
        gamma = B(1); % shape parameter #1
        delta = B(2); % shape parameter #2
        mi = B(3); % location parameter
        sigma = B(4); % scale parameter
        p = 0.5*(1+erf((gamma+delta*asinh((x-mi)/sigma))/(2^0.5)));
	case 'SL' % semi-bounded
        if size(B,1)~=4
            error('Incorrect number of parameters - Johnson SU requires 4 parameters')
        end
        gamma = B(1); % shape parameter #1
        delta = B(2); % shape parameter #2
        mi = B(3); % location parameter
        sigma = B(4); % scale parameter
        p = zeros(size(x));
        p((x > mi) & (x <= mi+sigma)) = 0.5*erfc(-(gamma+delta*log((x((x > mi) & (x <= mi+sigma))-mi)./sigma))/(2^0.5));
        p(x > mi + sigma) = 0.5*(1+erf((gamma+delta*log((x(x > mi + sigma)-mi)./sigma))/(2^0.5)));
        p(x <= mi) = 0; % Co dla x <=  mi? Czy wg Wolframa to 0? Bo to wg mnie jest
        % niejasne, tam jest s�owo "true" zamiast przedzia�u.
	case 'SL0' % semi-bounded (at 0)
        if size(B,1)~=3
            if size(B,1)==4 && B(3) == 0
                gamma = B(1); % shape parameter #1
                delta = B(2); % shape parameter #2
                mi = B(3); % location parameter
                sigma = B(4); % scale parameter
            else
                error('Incorrect number of parameters - Johnson SL0 requires 3 parameters')
            end
        else
            gamma = B(1); % shape parameter #1
            delta = B(2); % shape parameter #2
            mi = 0; % location parameter
            sigma = B(3); % scale parameter            
        end        
        p = zeros(size(x));
        % Rozumiem, �e mo�na tu zostawi� mi, bo zostanie za to podstawione
        % 0, prawda? To chyba taki sam wz�r jak na SL, tak? Skopiowa�am.
        p((x > mi) & (x <= mi+sigma)) = 0.5*erfc(-(gamma+delta*log((x((x > mi) & (x <= mi+sigma))-mi)./sigma))/(2^0.5));
        p(x > mi + sigma) = 0.5*(1+erf((gamma+delta*log((x(x > mi + sigma)-mi)./sigma))/(2^0.5)));
        p(x <= mi) = 0; % Zgadza si�?
	case 'SB' % bounded
        if size(B,1)~=4
            error('Incorrect number of parameters - Johnson SU requires 4 parameters')
        end
        gamma = B(1); % shape parameter #1
        delta = B(2); % shape parameter #2
        mi = B(3); % location parameter
        sigma = B(4); % scale parameter
        p = zeros(size(x));
        p((x > mi) & (x < mi+sigma/2)) = 0.5*erfc(-(gamma+delta*log((x((x > mi) & (x < mi+sigma/2))-mi)./(-x((x > mi) & (x < mi+sigma/2))+mi+sigma)))/(2^0.5));
        p((x >= mi + sigma/2) & (x < mi + sigma)) = 0.5*(1+erf((gamma+delta*log((x((x >= mi + sigma/2) & (x < mi + sigma))-mi)./(-x((x >= mi + sigma/2) & (x < mi + sigma))+mi+sigma)))/(2^0.5)));
        p(x >=  mi + sigma) = 1;    
        p(x <= mi) = 0; % Zgadza si�?
   case 'SB0' % bounded (at 0)
        if size(B,1)~=3
            if size(B,1)==4 && B(3) == 0
                gamma = B(1); % shape parameter #1
                delta = B(2); % shape parameter #2
                mi = B(3); % location parameter
                sigma = B(4); % scale parameter
            else
                error('Incorrect number of parameters - Johnson SL0 requires 3 parameters')
            end
        else
            gamma = B(1); % shape parameter #1
            delta = B(2); % shape parameter #2
            mi = 0; % location parameter
            sigma = B(3); % scale parameter            
        end        
        % Znowu przeklei�am z SB, bo to chyba to samo...
        p = zeros(size(x));
        p((x > mi) & (x < mi+sigma/2)) = 0.5*erfc(-(gamma+delta*log((x((x > mi) & (x < mi+sigma/2))-mi)./(-x((x > mi) & (x < mi+sigma/2))+mi+sigma)))/(2^0.5));
        p((x >= mi + sigma/2) & (x < mi + sigma)) = 0.5*(1+erf((gamma+delta*log((x((x >= mi + sigma/2) & (x < mi + sigma))-mi)./(-x((x >= mi + sigma/2) & (x < mi + sigma))+mi+sigma)))/(2^0.5)));
        p(x >=  mi + sigma) = 1;    
        p(x <= mi) = 0; % Zgadza si�?
   otherwise
      error('Unknown distribution type');
end