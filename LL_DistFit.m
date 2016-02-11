function f = LL_DistFit(bounds, dist, b0)

% save CDF_WTP_tmp
% return

switch dist

% unbounded
    case 0 % Normal % mu, sigma
        p0 = cdf('Normal',bounds,b0(1),b0(2));
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('Normal',bounds(p==0,1),b0(1),b0(2));
        f = log(p);
    case 1 % Logistic % mu, sigma
        p0 = cdf('Logistic',bounds,b0(1),b0(2)); 
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('Logistic',bounds(p==0,1),b0(1),b0(2));
        f = log(p);
    case 2 % Extreme Value % mu sigma
        p0 = cdf('Extreme Value',bounds,b0(1),b0(2));        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('Extreme Value',bounds(p==0,1),b0(1),b0(2));
        f = log(p);
    case 3 % Generalized Extreme Value % k, sigma, mu 
        p0 = cdf('Generalized Extreme Value',bounds,b0(1),b0(2),b0(3));        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('Generalized Extreme Value',bounds(p==0,1),b0(1),b0(2),b0(3));
        f = log(p);
    case 4 % tLocationScale % mu, sigma, nu 
        p0 = cdf('tLocationScale',bounds,b0(1),b0(2),b0(3));        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('tLocationScale',bounds(p==0,1),b0(1),b0(2),b0(3));
        f = log(p);
    case 5 % Uniform % min, max
        p0 = cdf('Uniform',bounds,b0(1),b0(2));        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('Uniform',bounds(p==0,1),b0(1),b0(2));
        f = log(p);
    case 6 % Johnson SU % gamma delta xi lambda
        p0 = [f_johnson_cdf(bounds(:,1),b0,'SU'),f_johnson_cdf(bounds(:,2),b0,'SU')];        
        p = p0(:,2) - p0(:,1);
        p(p==0) = f_johnson_pdf(bounds(p==0,1),b0,'SU');
        f = log(p);
    % stable
        
% bounded (0,Inf)
    case 10 % Exponential % mu
        p0 = cdf('Exponential',bounds,b0);
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('Exponential',bounds(p==0,1),b0);
        f = log(p);
    case 11 % Lognormal % mu, sigma
        p0 = cdf('Lognormal',bounds,b0(1),b0(2));        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('Lognormal',bounds(p==0,1),b0(1),b0(2));
        f = log(p);
    case 12 % Loglogistic % mu, sigma
        p0 = cdf('LogLogistic',bounds,b0(1),b0(2));        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('LogLogistic',bounds(p==0,1),b0(1),b0(2));
        f = log(p);
    case 13 % Weibull % A, b0
        p0 = cdf('Weibull',bounds,b0(1),b0(2));        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('Weibull',bounds(p==0,1),b0(1),b0(2));
        f = log(p);
    case 14 % Rayleigh % b0
        p0 = cdf('Rayleigh',bounds,b0);        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('Rayleigh',bounds(p==0,1),b0);
        f = log(p);
    case 15 % Gamma % a, b
        p0 = cdf('Gamma',bounds,b0(1),b0(2));        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('Gamma',bounds(p==0,1),b0(1),b0(2));
        f = log(p);
    case 16 % BirnbaumSaunders % beta, gamma
        p0 = cdf('BirnbaumSaunders',bounds,b0(1),b0(2));        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('BirnbaumSaunders',bounds(p==0,1),b0(1),b0(2));
        f = log(p);
    case 17 % Generalized Pareto % k, sigma, theta
        p0 = cdf('Generalized Pareto',bounds,b0(1),b0(2),b0(3));        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('Generalized Pareto',bounds(p==0,1),b0(1),b0(2),b0(3));
        f = log(p);
    case 18 % InverseGaussian % k, sigma
        p0 = cdf('InverseGaussian',bounds,b0(1),b0(2));        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('InverseGaussian',bounds(p==0,1),b0(1),b0(2));
        f = log(p);
    case 19 % Nakagami % mu, omega
        p0 = cdf('Nakagami',bounds,b0(1),b0(2));        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('Nakagami',bounds(p==0,1),b0(1),b0(2));
        f = log(p);
    case 20 % Rician % s, sigma
        p0 = cdf('Rician',bounds,b0(1),b0(2));        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('Rician',bounds(p==0,1),b0(1),b0(2));
        f = log(p);
%     case 21 % Johnson SB % gamma delta xi lambda 
%         save tmp2
%         p0 = [f_johnson_cdf(bounds(:,1),b0,'SB'),f_johnson_cdf(bounds(:,2),b0,'SB')];        
%         p = p0(:,2) - p0(:,1);
%         p(p==0) = f_johnson_pdf(bounds(p==0,1),b0,'SB');
%         f = log(p);
        %     case x % Burr % Error - Parto distribution fits better
        %           p = cdf('Burr',bounds,b0(1),b0(2),b0(3));        
        %           f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
        %     elseif dist==23 % generalized inverse Gaussian
        % %         p =
        %     elseif dist==24 % sinh-arcsinh
        % %         p =
        
 % discrete
    case 31 % Poisson % lambda
        p0 = cdf('Poisson',bounds,b0);        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('Poisson',bounds(p==0,1),b0);
        f = log(p);
    case 32 % Negative Binomial % R, P
        p0 = cdf('Negative Binomial',bounds,b0(1),b0(2));        
        p = p0(:,2) - p0(:,1);
        p(p==0) = pdf('Negative Binomial',bounds(p==0,1),b0(1),b0(2));
        f = log(p);
end

end

