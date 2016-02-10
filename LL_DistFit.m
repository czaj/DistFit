function f = LL_DistFit(bounds, dist, RealMin, b0)

% save CDF_WTP_tmp
% return

switch dist

% unbounded
    case 0 % Nnormal % mu, sigma
        p = cdf('Normal',bounds,b0(1),b0(2));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 1 % Logistic % mu, sigma
        p = cdf('Logistic',bounds,b0(1),b0(2));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 2 % Extreme Value % mu sigma
        p = cdf('Extreme Value',bounds,b0(1),b0(2));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 3 % Generalized Extreme Value % k, sigma, mu 
        p = cdf('Generalized Extreme Value',bounds,b0(1),b0(2),b0(3));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 4 % tLocationScale % mu, sigma, nu 
        p = cdf('tLocationScale',bounds,b0(1),b0(2),b0(3));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 5 % Uniform % min, max
        p = cdf('Uniform',bounds,b0(1),b0(2));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 6 % Johnson SU % gamma delta xi lambda
        p = [f_johnson_cdf(bounds(:,1),b0, 'SU'),f_johnson_cdf(bounds(:,2),b0, 'SU')];        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    % stable
        
% bounded (0,Inf)
    case 10 % Exponential % mu
        p = cdf('Exponential',bounds,b0);
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 11 % Lognormal % mu, sigma
        p = cdf('Lognormal',bounds,b0(1),b0(2));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 12 % Loglogistic % mu, sigma
        p = cdf('LogLogistic',bounds,b0(1),b0(2));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 13 % Weibull % A, b0
        p = cdf('Weibull',bounds,b0(1),b0(2));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 14 % Rayleigh % b0
        p = cdf('Rayleigh',bounds,b0);        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 15 % Gamma % a, b
        p = cdf('Gamma',bounds,b0(1),b0(2));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 16 % BirnbaumSaunders % beta, gamma
        p = cdf('BirnbaumSaunders',bounds,b0(1),b0(2));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 17 % Generalized Pareto % k, sigma, theta
        p = cdf('Generalized Pareto',bounds,b0(1),b0(2),b0(3));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 18 % InverseGaussian % k, sigma
        p = cdf('InverseGaussian',bounds,b0(1),b0(2));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 19 % Nakagami % mu, omega
        p = cdf('Nakagami',bounds,b0(1),b0(2));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 20 % Rician % s, sigma
        p = cdf('Rician',bounds,b0(1),b0(2));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 21 % Johnson SB % gamma delta xi lambda % ERROR - values out of range!
        p = [f_johnson_cdf(bounds(:,1),b0, 'SB'),f_johnson_cdf(bounds(:,2),b0, 'SB')];        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
        %     case x % Burr % Error - Parto distribution fits better
        %           p = cdf('Burr',bounds,b0(1),b0(2),b0(3));        
        %           f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
        %     elseif dist==23 % generalized inverse Gaussian
        % %         p =
        %     elseif dist==24 % sinh-arcsinh
        % %         p =
        
 % discrete
    case 31 % Poisson % lambda
        p = cdf('Poisson',bounds,b0);        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 32 % Negative Binomial % R, P
        p = cdf('Negative Binomial',bounds,b0(1),b0(2));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
end

end

