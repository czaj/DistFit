function f = LL_DistFit(bounds, dist, RealMin, b0)

% save CDF_WTP_tmp
% return

% czy mo?esz wszystkie poprawi? tak jak rozk?ad lognormalny? 
% to znaczy niech u?ywaj? jednolicie funkcji cdf i do tego tej opcji RealMin

switch dist

% unbounded
    case 0 % normal % mu, sigma
        f = log(normcdf(bounds(:,2),b0(1),b0(2)) - normcdf(bounds(:,1),b0(1),b0(2)));        
    case 1 % logistic
        pd = fitdist(midpoint,'Logistic','Options',OptimOptFit);
        b0 = pd.ParameterValues; % mu, sigma
    case 2 % Extreme Value
        pd = fitdist(midpoint,'ExtremeValue','Options',OptimOptFit);
        b0 = pd.ParameterValues; % mu sigma
    case 3 % Generalized Extreme Value
        pd = fitdist(midpoint,'GeneralizedExtremeValue','Options',OptimOptFit);
        b0 = pd.ParameterValues; % k, sigma, mu
    case 4 % tLocationScale
        pd = fitdist(midpoint,'tLocationScale','Options',OptimOptFit);
        b0 = pd.ParameterValues; % mu, sigma, nu
    case 5 % uniform
        b0 = [min(midpoint),max(midpoint)];
    case 6 % Rician
        pd = fitdist(midpoint,'Rician','Options',OptimOptFit);
        b0 = pd.ParameterValues; % s, sigma
    case 7 % Johnson SU
        pd = f_johnson_fit(midpoint);
        b0 =  pd.coef; % gamma delta xi lambda
    % stable
        
% bounded (0,Inf)
    case 10 % exponential
        pd = fitdist(midpoint,'Exponential','Options',OptimOptFit);
        b0 = pd.ParameterValues; % mu
    case 12 % lognormal
        p = cdf('Lognormal',bounds,b0(1),b0(2));        
        f = RealMin*log(max(eps,p(:,2) - p(:,1))) + ~RealMin*log(p(:,2) - p(:,1));
    case 13 % loglogistic
        pd = fitdist(midpoint,'Loglogistic','Options',OptimOptFit);
        b0 = pd.ParameterValues; % mu, sigma
    case 14 % Weibull
        pd = fitdist(midpoint,'Weibull','Options',OptimOptFit);
        b0 = pd.ParameterValues; % A, b0
    case 11 % Rayleigh
        pd = fitdist(midpoint,'Rayleigh','Options',OptimOptFit);
        b0 = pd.ParameterValues; % b0
    case 15 % Gamma
        pd = fitdist(midpoint,'Gamma','Options',OptimOptFit);
        b0 = pd.ParameterValues; % a, b
    case 16 % b0irnbaumSaunders
        pd = fitdist(midpoint,'b0irnbaumSaunders','Options',OptimOptFit);
        b0 = pd.ParameterValues; % beta, gamma
    case 17 % Generalized Pareto
        pd = fitdist(midpoint,'GeneralizedPareto','Options',OptimOptFit);
        b0 = pd.ParameterValues; % k, sigma, theta
    case 18 % InverseGaussian
        pd = fitdist(midpoint,'InverseGaussian','Options',OptimOptFit);
        b0 = pd.ParameterValues; % k, sigma, theta
    case 19 % Nakagami
        pd = fitdist(midpoint,'Nakagami','Options',OptimOptFit);
        b0 = pd.ParameterValues; % mu, omega
    case 20 % Rician
        pd = fitdist(midpoint,'Rician','Options',OptimOptFit);
        b0 = pd.ParameterValues; % s, sigma
    case 21 % Johnson Sb0
        pd = f_johnson_fit(midpoint);
        b0 = pd.coef;
        %     case x % b0urr
        %         pd = fitdist(midpoint,'b0urr','Options',OptimOptFit); % Error - Parto distribution fits better
        %         b0 = pd.ParameterValues; %
        %     elseif dist==23 % generalized inverse Gaussian
        % %         b0 =
        %     elseif dist==24 % sinh-arcsinh
        % %         b0 =
        
        % discrete
    case 31 % Poisson
        pd = fitdist(midpoint,'Poisson','Options',OptimOptFit);
        b0 = pd.ParameterValues; % lambda
    case 32 % negative binomial
        pd = fitdist(round(midpoint),'Negativeb0inomial','Options',OptimOptFit);
        b0 = pd.ParameterValues; % R, P
end

    
%     
% elseif dist==1
% elseif dist==2
%     f = log(cdf('Exponential',bounds(:,2),b0) - cdf('Exponential',bounds(:,1),b0));
% elseif dist==3
%     f = log(cdf('Logistic',bounds(:,2),b0(1),b0(2)) - cdf('Logistic',bounds(:,1),b0(1),b0(2)));
% elseif dist==4
%     f = log(cdf('LogLogistic',bounds(:,2),b0(1),b0(2)) - cdf('LogLogistic',bounds(:,1),b0(1),b0(2)));
% elseif dist==5 
% %     f = log(cdf('Poisson',bounds(:,2),b0) - cdf('Poisson',bounds(:,1),b0));
%     f = log(max(cdf('Poisson',bounds(:,2),b0) - cdf('Poisson',bounds(:,1),b0),realmin));
% elseif dist==6
%     f = log(cdf('Weibull',bounds(:,2),b0(1),b0(2)) - cdf('Weibull',bounds(:,1),b0(1),b0(2)));
% elseif dist==7
%     f = log(cdf('Gamma',bounds(:,2),b0(1),b0(2)) - cdf('Gamma',bounds(:,1),b0(1),b0(2)));
% elseif dist==8
% %     f = log(cdf('Rayleigh',bounds(:,2),b0) - cdf('Rayleigh',bounds(:,1),b0));
%     f = log(max(cdf('Rayleigh',bounds(:,2),b0) - cdf('Rayleigh',bounds(:,1),b0),realmin));
% elseif dist==9
%     f = log(cdf('b0irnbaumSaunders',bounds(:,2),b0(1),b0(2)) - cdf('b0irnbaumSaunders',bounds(:,1),b0(1),b0(2)));
% elseif dist==10
%     f = log(cdf('ExtremeValue',bounds(:,2),b0(1),b0(2)) - cdf('ExtremeValue',bounds(:,1),b0(1),b0(2)));
% elseif dist==11
%     f = log(cdf('GeneralizedExtremeValue',bounds(:,2),b0(1),b0(2),b0(3)) - cdf('GeneralizedExtremeValue',bounds(:,1),b0(1),b0(2),b0(3)));
% elseif dist==12
%     f = log(cdf('GeneralizedPareto',bounds(:,2),b0(1),b0(2)) - cdf('GeneralizedPareto',bounds(:,1),b0(1),b0(2)));
% elseif dist==13
%     f = log(cdf('InverseGaussian',bounds(:,2),b0(1),b0(2)) - cdf('InverseGaussian',bounds(:,1),b0(1),b0(2)));
% elseif dist==14
%     f = log(cdf('Nakagami',bounds(:,2),b0(1),b0(2)) - cdf('Nakagami',bounds(:,1),b0(1),b0(2)));
% elseif dist==15
%     f = log(cdf('Rician',bounds(:,2),b0(1),b0(2)) - cdf('Rician',bounds(:,1),b0(1),b0(2)));
% elseif dist==16
%     f = log(cdf('tLocationScale',bounds(:,2),b0(1),b0(2),b0(3)) - cdf('tLocationScale',bounds(:,1),b0(1),b0(2),b0(3)));
% elseif dist==17 % Johnson SU
%     f = log(f_johnson_cdf(bounds(:,2),b0, 'SU') - f_johnson_cdf(bounds(:,1),b0, 'SU'));
% elseif dist==18 % Johnson Sb0 - nie dzia�a, gdy any(u <= 0) || any(u>=1.0)
%     f = log(f_johnson_cdf(bounds(:,2),b0, 'Sb0') - f_johnson_cdf(bounds(:,1),b0, 'Sb0'));
% elseif dist==19 % b0urr
%     f = log(cdf('b0urr',bounds(:,2),b0(1),b0(2),b0(3)) - cdf('b0urr',bounds(:,1),b0(1),b0(2),b0(3)));
% elseif dist==20 % F
%     f = log(cdf('F',bounds(:,2),b0(1),b0(2)) - cdf('F',bounds(:,1),b0(1),b0(2)));
% elseif dist==21 % negative binomial
%     f = log(cdf('Negative b0inomial',bounds(:,2),b0(1),b0(2)) - cdf('Negative b0inomial',bounds(:,1),b0(1),b0(2)));
% elseif dist==22 % uniform
%     f = log(cdf('Uniform',bounds(:,2),b0(1),b0(2)) - cdf('Uniform',bounds(:,1),b0(1),b0(2)));
% elseif dist==23 % generalized inverse Gaussian
% %     f = 
% elseif dist==24 % sinh-arcsinh
% %     f = 
end
