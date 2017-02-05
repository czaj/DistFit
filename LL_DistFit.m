function f = LL_DistFit(bounds, X, weights, dist, Spike, RealMin, b0)

% save CDF_WTP_tmp
% return

b0 = b0(:);

numDistParam = 1*any(dist == [10,14,31]) + 2*any(dist == [0:2,5,11:13,15,16,18:20,32]) + 3*any(dist == [3,4,17]) + 4*any(dist == [6,21:22]);

numX = size(X,2);

if Spike
    bSpike = b0(numDistParam+1);
    if numX > 0 % Spike and X
        bDist0 = b0(1:numDistParam);
        bCovDist = b0(numDistParam+1+1:numDistParam*(1+numX)+1);
        bDist(size(bounds,1),numDistParam) = 0;
        for i = 1:numDistParam
            bDist(:,i) = bDist0(i) + X*bCovDist(numX*(i-1)+1:numX*i);
        end
        bCovSpike = b0(numDistParam*(1+numX)+1+1:(numDistParam+1)*(1+numX));
        pSpike = normcdf(bSpike+X*bCovSpike);
    else % Spike only
        bDist = repmat(b0(1:numDistParam)',[size(bounds,1),1]); % distribution parameters
        pSpike = repmat(normcdf(bSpike),[size(bounds,1),1]);
    end
else
    if numX > 0 % X only
        bDist0 = b0(1:numDistParam);
        bCovDist = b0(numDistParam+1:numDistParam*(1+numX));
        bDist(size(bounds,1),numDistParam) = 0;
        for i = 1:numDistParam
            bDist(:,i) = bDist0(i) + X*bCovDist(numX*(i-1)+1:numX*i);
        end
    else % baseline distribution only
        bDist = repmat(b0(1:numDistParam)',[size(bounds,1),1]); % distribution parameters
    end
    pSpike(size(bounds,1),1) = 0;
end

[~,I] = min(abs(bounds),[],2);
bounds_min = bounds(sub2ind(size(bounds),(1:size(bounds,1)).',I)); % lower of the absolute value of bounds  

switch dist
    
    % unbounded
    case 0 % Normal % mu, sigma>0
        dp = cdf('Normal',bounds(:,2),bDist(:,1),bDist(:,2)) - ...
            cdf('Normal',bounds(:,1),bDist(:,1),bDist(:,2));
        dp(dp==0) = pdf('Normal',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2)); % replace 0 cdf difference with pdf at lower absolute bound (if extreme bounds or exact x known)
    case 1 % Logistic % mu, sigma>=0
        dp = cdf('Logistic',bounds(:,2),bDist(:,1),bDist(:,2)) - ...
            cdf('Logistic',bounds(:,1),bDist(:,1),bDist(:,2));
        dp(dp==0) = pdf('Logistic',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2));
    case 2 % Extreme Value % mu, sigma>0
        dp = cdf('Extreme Value',bounds(:,2),bDist(:,1),bDist(:,2)) - ...
            cdf('Extreme Value',bounds(:,1),bDist(:,1),bDist(:,2));
        dp(dp==0) = pdf('Extreme Value',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2));
    case 3 % Generalized Extreme Value % k, sigma>0, mu
        dp = cdf('Generalized Extreme Value',bounds(:,2),bDist(:,1),bDist(:,2),bDist(:,3)) - ...
            cdf('Generalized Extreme Value',bounds(:,1),bDist(:,1),bDist(:,2),bDist(:,3));
        dp(dp==0) = pdf('Generalized Extreme Value',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2),bDist(dp==0,3));
    case 4 % tLocationScale % mu, sigma>0, nu>0
        dp = cdf('tLocationScale',bounds(:,2),bDist(:,1),bDist(:,2),bDist(:,3)) - ...
            cdf('tLocationScale',bounds(:,1),bDist(:,1),bDist(:,2),bDist(:,3));
        dp(dp==0) = pdf('tLocationScale',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2),bDist(dp==0,3));
    case 5 % Uniform % min, max
        dp = cdf('Uniform',bounds(:,2),bDist(:,1),bDist(:,2)) - ...
            cdf('Uniform',bounds(:,1),bDist(:,1),bDist(:,2));
        dp(dp==0) = pdf('Uniform',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2));
    case 6 % Johnson SU % gamma, delta>0, mi, sigma>0
        %         save tmp1
        dp = JohnsonCDF(bounds(:,2),bDist(:,1),bDist(:,2),bDist(:,3),bDist(:,4),'SU') - ...
            JohnsonCDF(bounds(:,1),bDist(:,1),bDist(:,2),bDist(:,3),bDist(:,4),'SU');
        dp(dp==0) = JohnsonPDF(bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2),bDist(dp==0,3),bDist(dp==0,4),'SU');
        %         save tmp1
        %         return
        
        % stable
        
        % bounded (0,Inf)
    case 10 % Exponential % mu>0
        dp = cdf('Exponential',bounds(:,2),bDist(:,1)) - ...
            cdf('Exponential',bounds(:,1),bDist(:,1));
        dp(dp==0) = pdf('Exponential',bounds_min(dp==0,1),bDist(dp==0,1));
    case 11 % Lognormal % mu, sigma>0
        dp = cdf('Lognormal',bounds(:,2),bDist(:,1),bDist(:,2)) - ...
            cdf('Lognormal',bounds(:,1),bDist(:,1),bDist(:,2));
        dp(dp==0) = pdf('Lognormal',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2));
    case 12 % Loglogistic % mu>0, sigma>0
        dp = cdf('Loglogistic',bounds(:,2),bDist(:,1),bDist(:,2)) - ...
            cdf('Loglogistic',bounds(:,1),bDist(:,1),bDist(:,2));
        dp(dp==0) = pdf('Loglogistic',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2));
    case 13 % Weibull % A>0, b0>0
        dp = cdf('Weibull',bounds(:,2),bDist(:,1),bDist(:,2)) - ...
            cdf('Weibull',bounds(:,1),bDist(:,1),bDist(:,2));
        dp(dp==0) = pdf('Weibull',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2));
    case 14 % Rayleigh % b0>0
        dp = cdf('Rayleigh',bounds(:,2),bDist(:,1)) - ...
            cdf('Rayleigh',bounds(:,1),bDist(:,1));
        dp(dp==0) = pdf('Rayleigh',bounds_min(dp==0,1),bDist(dp==0,1));
    case 15 % Gamma % a>0, b>0
        dp = cdf('Gamma',bounds(:,2),bDist(:,1),bDist(:,2)) - ...
            cdf('Gamma',bounds(:,1),bDist(:,1),bDist(:,2));
        dp(dp==0) = pdf('Gamma',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2));
    case 16 % BirnbaumSaunders % beta>0, gamma>0
        dp = cdf('BirnbaumSaunders',bounds(:,2),bDist(:,1),bDist(:,2)) - ...
            cdf('BirnbaumSaunders',bounds(:,1),bDist(:,1),bDist(:,2));
        dp(dp==0) = pdf('BirnbaumSaunders',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2));
    case 17 % Generalized Pareto % k, sigma>0, theta
        dp = cdf('Generalized Pareto',bounds(:,2),bDist(:,1),bDist(:,2),bDist(:,3)) - ...
            cdf('Generalized Pareto',bounds(:,1),bDist(:,1),bDist(:,2),bDist(:,3));
        dp(dp==0) = pdf('Generalized Pareto',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2),bDist(dp==0,3));
    case 18 % InverseGaussian % k>0, sigma>0
        dp = cdf('InverseGaussian',bounds(:,2),bDist(:,1),bDist(:,2)) - ...
            cdf('InverseGaussian',bounds(:,1),bDist(:,1),bDist(:,2));
        dp(dp==0) = pdf('InverseGaussian',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2));
    case 19 % Nakagami % mu>=0.5, omega>0
        dp = cdf('Nakagami',bounds(:,2),bDist(:,1),bDist(:,2)) - ...
            cdf('Nakagami',bounds(:,1),bDist(:,1),bDist(:,2));
        dp(dp==0) = pdf('Nakagami',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2));
%         p = (1-pSpike).*dp;
%         p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) = p(bounds(:,1) <= 0 & 0 <= bounds(:,2)) + pSpike(bounds(:,1) <= 0 & 0 <= bounds(:,2));
%         f = log(p).*weights;
    case 20 % Rician % s>=0, sigma>0
        dp = cdf('Rician',bounds(:,2),bDist(:,1),bDist(:,2)) - ...
            cdf('Rician',bounds(:,1),bDist(:,1),bDist(:,2));
        dp(dp==0) = pdf('Rician',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2));
    case 21 % Johnson SB % gamma, delta>0, mi, sigma>0
        dp = JohnsonCDF(bounds(:,2),bDist(:,1),bDist(:,2),bDist(:,3),bDist(:,4),'SB') - ...
            JohnsonCDF(bounds(:,1),bDist(:,1),bDist(:,2),bDist(:,3),bDist(:,4),'SB');
        dp(dp==0) = JohnsonPDF(bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2),bDist(dp==0,3),bDist(dp==0,4),'SB');
    case 22 % Johnson SL % gamma, delta>0, mi, sigma>0
        dp = JohnsonCDF(bounds(:,2),bDist(:,1),bDist(:,2),bDist(:,3),bDist(:,4),'SL') - ...
            JohnsonCDF(bounds(:,1),bDist(:,1),bDist(:,2),bDist(:,3),bDist(:,4),'SL');
        dp(dp==0) = JohnsonPDF(bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2),bDist(dp==0,3),bDist(dp==0,4),'SL');
        
        % discrete
    case 31 % Poisson % lambda>0
        dp = cdf('Poisson',bounds(:,2),bDist(:,1)) - ...
            cdf('Poisson',bounds(:,1),bDist(:,1));
        dp(dp==0) = pdf('Poisson',bounds_min(dp==0,1),bDist(dp==0,1));
    case 32 % Negative Binomial % R>0 and R is an integer, 0<P<1
        dp = cdf('Negative Binomial',bounds(:,2),bDist(:,1),bDist(:,2)) - ...
            cdf('Negative Binomial',bounds(:,1),bDist(:,1),bDist(:,2));
        dp(dp==0) = pdf('Negative Binomial',bounds_min(dp==0,1),bDist(dp==0,1),bDist(dp==0,2));       
end

% dp(dp == Inf) = 1 - eps;

p = (1-pSpike).*dp; % scale down to allow for the probability of spike
I0 = ((bounds(:,1) == 0 & bounds(:,2) == 0) | (bounds(:,1) <= 0 & bounds(:,2) > 0));
p(I0) = p(I0) + pSpike(I0); % add spike probability to observations with 0 in bounds

if RealMin == 0
    f = log(p).*weights;
else
f = log(max(p,eps)).*weights;
% f = log(max(p,realmin)).*weights;
end



