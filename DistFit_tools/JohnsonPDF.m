function p = JohnsonPDF(x,B,type)

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
        p = (exp(-0.5*(gamma+delta*arcsinh((x-mi)/sigma))^2)*delta)./((2*pi)^(0.5)*((x-mi)^2+sigma^2)^(0.5));
	case 'SL' % semi-bounded
        if size(B,1)~=4
            error('Incorrect number of parameters - Johnson SU requires 4 parameters')
        end
        gamma = B(1); % shape parameter #1
        delta = B(2); % shape parameter #2
        mi = B(3); % location parameter
        sigma = B(4); % scale parameter
        p = zeros(size(x));
        p(x > mi) = (exp(-0.5*(gamma+delta*arcsinh((x(x > mi)-mi)/sigma))^2)*delta)./((2*pi)^(0.5)*(x(x > mi)-mi));
        p(x <= mi) = 0; % Czy to siê zgadza? 
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
        % Pozostawi³am mi i skopiowa³am z SL.
        p(x > mi) = (exp(-0.5*(gamma+delta*arcsinh((x(x > mi)-mi)/sigma))^2)*delta)./((2*pi)^(0.5)*(x(x > mi)-mi));
        p(x <= mi) = 0; % Czy to siê zgadza? 
	case 'SB' % bounded
        if size(B,1)~=4
            error('Incorrect number of parameters - Johnson SU requires 4 parameters')
        end
        gamma = B(1); % shape parameter #1
        delta = B(2); % shape parameter #2
        mi = B(3); % location parameter
        sigma = B(4); % scale parameter
        p = zeros(size(x));
        p((x > mi) & (x < mi + sigma)) = (exp(-0.5*(gamma+delta*log((x((x > mi) & (x < mi + sigma))-mi)./(-x((x > mi) & (x < mi + sigma))+mi+sigma)))^2)*delta*sigma)./((2*pi)^(0.5)*(x((x > mi) & (x < mi + sigma))-mi)*(-x((x > mi) & (x < mi + sigma))+mi+sigma));
        p((x <= mi) | (x >= mi + sigma)) = 0; % Czy to siê zgadza? 
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
        % Znowu przeklei³am z SB, bo to chyba to samo...
        p = zeros(size(x));
        p((x > mi) & (x < mi + sigma)) = (exp(-0.5*(gamma+delta*log((x((x > mi) & (x < mi + sigma))-mi)./(-x((x > mi) & (x < mi + sigma))+mi+sigma)))^2)*delta*sigma)./((2*pi)^(0.5)*(x((x > mi) & (x < mi + sigma))-mi)*(-x((x > mi) & (x < mi + sigma))+mi+sigma));
        p((x <= mi) | (x >= mi + sigma)) = 0; % Czy to siê zgadza? 
   otherwise
      error('Unknown distribution type');
end