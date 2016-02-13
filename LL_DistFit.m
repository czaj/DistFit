function f = LL_DistFit(bounds, dist, Spike, b0)

% save CDF_WTP_tmp
% return

if Spike
    pSpike = normcdf(b0(end));
else
    pSpike = 0;
end

switch dist
    
% unbounded 
    case 0 % Normal % mu, sigma
        p0 = cdf('Normal',bounds,b0(1),b0(2));
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('Normal',bounds(dp==0,1),b0(1),b0(2));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 1 % Logistic % mu, sigma
        p0 = cdf('Logistic',bounds,b0(1),b0(2)); 
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('Logistic',bounds(dp==0,1),b0(1),b0(2));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 2 % Extreme Value % mu sigma
        p0 = cdf('Extreme Value',bounds,b0(1),b0(2));        
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('Extreme Value',bounds(dp==0,1),b0(1),b0(2));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 3 % Generalized Extreme Value % k, sigma, mu 
        p0 = cdf('Generalized Extreme Value',bounds,b0(1),b0(2),b0(3));        
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('Generalized Extreme Value',bounds(dp==0,1),b0(1),b0(2),b0(3));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 4 % tLocationScale % mu, sigma, nu 
        p0 = cdf('tLocationScale',bounds,b0(1),b0(2),b0(3));        
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('tLocationScale',bounds(dp==0,1),b0(1),b0(2),b0(3));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 5 % Uniform % min, max
        p0 = cdf('Uniform',bounds,b0(1),b0(2));        
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('Uniform',bounds(dp==0,1),b0(1),b0(2));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 6 % Johnson SU % gamma delta xi lambda
%         save tmp1
%         return
        p0 = JohnsonCDF(bounds,b0(1:4),'SB');
%         ... 

%         p0 = [f_johnson_cdf(bounds(:,1),b0,'SU'),f_johnson_cdf(bounds(:,2),b0,'SU')];        
%         dp = p0(:,2) - p0(:,1);
%         dp(dp==0) = f_johnson_pdf(bounds(dp==0,1),b0,'SU');
%         p = (1-pSpike)*dp;
%         p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
%         f = log(p);
    % stable
        
% bounded (0,Inf)
    case 10 % Exponential % mu
        p0 = cdf('Exponential',bounds,b0(1));
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('Exponential',bounds(dp==0,1),b0(1));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 11 % Lognormal % mu, sigma
        p0 = cdf('Lognormal',bounds,b0(1),b0(2));        
        dp = p0(:,2) - p0(:,1);        
        dp(dp==0) = pdf('Lognormal',bounds(dp==0,1),b0(1),b0(2));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);        
    case 12 % Loglogistic % mu, sigma
        p0 = cdf('LogLogistic',bounds,b0(1),b0(2));        
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('LogLogistic',bounds(dp==0,1),b0(1),b0(2));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 13 % Weibull % A, b0
        p0 = cdf('Weibull',bounds,b0(1),b0(2));        
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('Weibull',bounds(dp==0,1),b0(1),b0(2));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 14 % Rayleigh % b0
        p0 = cdf('Rayleigh',bounds,b0(1));        
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('Rayleigh',bounds(dp==0,1),b0(1));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 15 % Gamma % a, b
        p0 = cdf('Gamma',bounds,b0(1),b0(2));        
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('Gamma',bounds(dp==0,1),b0(1),b0(2));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 16 % BirnbaumSaunders % beta, gamma
        p0 = cdf('BirnbaumSaunders',bounds,b0(1),b0(2));        
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('BirnbaumSaunders',bounds(dp==0,1),b0(1),b0(2));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 17 % Generalized Pareto % k, sigma, theta
        p0 = cdf('Generalized Pareto',bounds,b0(1),b0(2),b0(3));        
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('Generalized Pareto',bounds(dp==0,1),b0(1),b0(2),b0(3));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 18 % InverseGaussian % k, sigma
        p0 = cdf('InverseGaussian',bounds,b0(1),b0(2));        
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('InverseGaussian',bounds(dp==0,1),b0(1),b0(2));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 19 % Nakagami % mu, omega
        p0 = cdf('Nakagami',bounds,b0(1),b0(2));        
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('Nakagami',bounds(dp==0,1),b0(1),b0(2));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 20 % Rician % s, sigma
        p0 = cdf('Rician',bounds,b0(1),b0(2));        
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('Rician',bounds(dp==0,1),b0(1),b0(2));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
%     case 21 % Johnson SB % gamma delta xi lambda 
%         save tmp2
%         p0 = [f_johnson_cdf(bounds(:,1),b0,'SB'),f_johnson_cdf(bounds(:,2),b0,'SB')];        
%         p = normcdf(b0(5))*(p0(:,2) - p0(:,1))+(1-normcdf(b0(5)));
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
        p0 = cdf('Poisson',bounds,b0(1));        
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('Poisson',bounds(dp==0,1),b0(1));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
    case 32 % Negative Binomial % R, P
        p0 = cdf('Negative Binomial',bounds,b0(1),b0(2));        
        dp = p0(:,2) - p0(:,1);
        dp(dp==0) = pdf('Negative Binomial',bounds(dp==0,1),b0(1),b0(2));
        p = (1-pSpike)*dp;
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike;
        f = log(p);
end

end

