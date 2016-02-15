function f = LL_DistFit(bounds, X, dist, SpikeTrue, b0)

% save CDF_WTP_tmp
% return

b0 = b0(:);

numDistParam = 1*any(dist == [10,14,31]) + 2*any(dist == [0:2,5,11:13,15,16,18:20,32]) + 3*any(dist == [3,4,17]) + 4*any(dist == [6,21]);

XCovSize = size(X,2);
% XCovTrue = size(X,2)>0;

BDist = b0(1:numDistParam); % distribution parameters
if SpikeTrue 
    if XCovSize > 0 % Spike and X
        BpSpike = b0(numDistParam+1);
        BCovDist = b0(numDistParam+1+1:numDistParam*(1+XCovSize)+1);
        BCovSpike = b0(numDistParam*(1+XCovSize)+1+1:(numDistParam+1)*(1+XCovSize));
    else % Spike only
        BpSpike = b0(numDistParam+1);
        BCovDist = zeros(numDistParam*XCovSize,0); % []
        BCovSpike = zeros(XCovSize,0); % []
    end
    pSpike = normcdf(BpSpike+X*BCovSpike);
else 
    if XCovSize > 0 % X only
%         BpSpike = zeros(1,0); % []
        BCovDist = b0(numDistParam+1:numDistParam*(1+XCovSize));
%         BCovSpike = zeros(XCovSize,0); % []
    else % baseline distribution only
%         BpSpike = zeros(1,0); % []
        BCovDist = zeros(numDistParam*XCovSize,0); % []
%         BCovSpike = zeros(XCovSize,0); % []
    end
    pSpike = 0;
end


switch dist
    
% unbounded 
    case 0 % Normal % mu, sigma
%         p0 = cdf('Normal',bounds,BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2));
        dp = cdf('Normal',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2)) - ...
             cdf('Normal',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2));
%        bounds_min = (min(abs(bounds.')) .* (2*(abs(min(bounds.')) < max(bounds.'))-1)).'; 
        [~,I] = min(abs(bounds),[],2); ... % lower of the absolute value of bounds
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Normal',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2)); % replace 0 cdf difference with pdf at lower absolute bound (if extreme bounds or exact x known)
        p = (1-pSpike).*dp; % scale down to allow for the probability of spike
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2)); % add spike probability to observations with 0 in bounds
        f = log(p);
    case 1 % Logistic % mu, sigma
        dp = cdf('Logistic',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2)) - ...
             cdf('Logistic',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Logistic',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);
    case 2 % Extreme Value % mu sigma
        dp = cdf('Extreme Value',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2)) - ...
             cdf('Extreme Value',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Extreme Value',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);
    case 3 % Generalized Extreme Value % k, sigma, mu 
        dp = cdf('Generalized Extreme Value',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2),BDist(3)+X*BCovDist(XCovSize*2+1:XCovSize*3)) - ...
             cdf('Generalized Extreme Value',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2),BDist(3)+X*BCovDist(XCovSize*2+1:XCovSize*3));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Generalized Extreme Value',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2),BDist(3)+X(dp==0,:)*BCovDist(XCovSize*2+1:XCovSize*3)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);
    case 4 % tLocationScale % mu, sigma, nu 
        dp = cdf('tLocationScale',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2),BDist(3)+X*BCovDist(XCovSize*2+1:XCovSize*3)) - ...
             cdf('tLocationScale',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2),BDist(3)+X*BCovDist(XCovSize*2+1:XCovSize*3));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('tLocationScale',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2),BDist(3)+X(dp==0,:)*BCovDist(XCovSize*2+1:XCovSize*3)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);
    case 5 % Uniform % min, max
        dp = cdf('Uniform',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2)) - ...
             cdf('Uniform',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Uniform',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
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
        dp = cdf('Exponential',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize)) - ...
             cdf('Exponential',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Exponential',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);
    case 11 % Lognormal % mu, sigma
        dp = cdf('Lognormal',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2)) - ...
             cdf('Lognormal',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Lognormal',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);        
    case 12 % Loglogistic % mu, sigma
        dp = cdf('Loglogistic',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2)) - ...
             cdf('Loglogistic',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Loglogistic',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);
    case 13 % Weibull % A, b0
        dp = cdf('Weibull',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2)) - ...
             cdf('Weibull',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Weibull',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);
    case 14 % Rayleigh % b0
        dp = cdf('Rayleigh',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize)) - ...
             cdf('Rayleigh',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Rayleigh',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);
    case 15 % Gamma % a, b
        dp = cdf('Gamma',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2)) - ...
             cdf('Gamma',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Gamma',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);
    case 16 % BirnbaumSaunders % beta, gamma
        dp = cdf('BirnbaumSaunders',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2)) - ...
             cdf('BirnbaumSaunders',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('BirnbaumSaunders',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);
    case 17 % Generalized Pareto % k, sigma, theta
        dp = cdf('Generalized Pareto',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2),BDist(3)+X*BCovDist(XCovSize*2+1:XCovSize*3)) - ...
             cdf('Generalized Pareto',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2),BDist(3)+X*BCovDist(XCovSize*2+1:XCovSize*3));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Generalized Pareto',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2),BDist(3)+X(dp==0,:)*BCovDist(XCovSize*2+1:XCovSize*3)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);
    case 18 % InverseGaussian % k, sigma
        dp = cdf('InverseGaussian',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2)) - ...
             cdf('InverseGaussian',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('InverseGaussian',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);
    case 19 % Nakagami % mu, omega
        dp = cdf('Nakagami',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2)) - ...
             cdf('Nakagami',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Nakagami',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);
    case 20 % Rician % s, sigma
        dp = cdf('Rician',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2)) - ...
             cdf('Rician',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Rician',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
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
        dp = cdf('Poisson',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize)) - ...
             cdf('Poisson',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Poisson',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);
    case 32 % Negative Binomial % R, P
        dp = cdf('Negative Binomial',bounds(:,2),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2)) - ...
             cdf('Negative Binomial',bounds(:,1),BDist(1)+X*BCovDist(1:XCovSize),BDist(2)+X*BCovDist(XCovSize+1:XCovSize*2));
        [~,I] = min(abs(bounds),[],2); ... 
        bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); 
        dp(dp==0) = pdf('Negative Binomial',bounds_min(dp==0,1),BDist(1)+X(dp==0,:)*BCovDist(1:XCovSize),BDist(2)+X(dp==0,:)*BCovDist(XCovSize+1:XCovSize*2)); 
        p = (1-pSpike).*dp; 
        p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
        f = log(p);
end

end

