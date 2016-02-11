function WTP = DistFit(INPUT, varargin)


%% check no. of inputs

if nargin < 1 
	error('Too few input arguments for WTPfit(INPUT,dist,b0,OptimOpt)')
elseif nargin == 1
    dist = [];
    b0 = [];
    EstimOpt = [];    
elseif nargin == 2
    dist = varargin{1};
    b0 = [];
    EstimOpt = [];
elseif nargin == 3
    dist = varargin{1};
    b0 = varargin{2};
    EstimOpt = [];
elseif nargin == 4
	dist = varargin{1};
    b0 = varargin{2};
    EstimOpt = varargin{3};
end
   
%% check INPUT

if size(INPUT.bounds,2) < 2
    error('Providing lower and upper bounds for WTP required')
elseif any(isnan(INPUT.bounds(:)))
    error('Some of the lower or upper bounds for WTP missing')
elseif any(sum(isfinite(INPUT.bounds),2)==0) % both bounds undefined
    error('Dataset includes observations with infinite bounds')
elseif any(INPUT.bounds(:,1) > INPUT.bounds(:,2))
    error('Some lower bounds are greater than upper bounds')
elseif (dist==[10:21] & INPUT.bounds(:)<=0)
    cprintf(rgb('DarkOrange'), 'WARNING: Adjusting bounds to match the distribution type \n')
    INPUT.bounds(INPUT.bounds(:,1)<=0,1)= eps;
    INPUT.bounds(INPUT.bounds(:,2)<=0,2)= INPUT.bounds(INPUT.bounds(:,2)<=0,1) + eps;
elseif (dist==[31:32] & INPUT.bounds(:)<0)
    cprintf(rgb('DarkOrange'), 'WARNING: Adjusting bounds to match the distribution type \n')
    INPUT.bounds(INPUT.bounds(:,1)<0,1)= 0;
    INPUT.bounds(INPUT.bounds(:,2)<0,2)= INPUT.bounds(INPUT.bounds(:,2)<0,1) + 1;
elseif (dist==[31:32] & any(isinteger(INPUT.bounds(:))==0))
    cprintf(rgb('DarkOrange'), 'WARNING: Adjusting bounds to match the distribution type \n')
    if round(INPUT.bounds(INPUT.bounds(:,1))==0,1)<0
        INPUT.bounds(isinteger(INPUT.bounds(:,1))==0,1) = 0;
    elseif round(INPUT.bounds(INPUT.bounds(:,1))==0,1)>=0
        INPUT.bounds(isinteger(INPUT.bounds(:,1))==0,1) = round(INPUT.bounds(INPUT.bounds(:,1))==0,1);
    elseif round(INPUT.bounds(INPUT.bounds(:,2))==0,2)<0
        INPUT.bounds(isinteger(INPUT.bounds(:,2))==0,2) = 0;
    elseif round(INPUT.bounds(INPUT.bounds(:,2))==0,2)>=0
        INPUT.bounds(isinteger(INPUT.bounds(:,2))==0,2) = round(INPUT.bounds(INPUT.bounds(:,2))==0,2);
    end
end
    
%% distribution type:

if ~exist('dist','var') || isempty(dist)
    disp('Assuming normally distributed WTP')
    dist = 0;
elseif ~any(dist == [0:7,10:21,31:32])
	error('Unsupported distribution type')
end

%% starting values:

if exist('b0','var') && ~isempty(b0)
    b0 = b0(:);
    if size(b0,1) ~= 1 || (size(b0,2) ~= 1 & dist == [10,14,31]) || (size(b0,2) ~= 2 & dist == [0:2,5,11:13,15,16,18:20,32]) || (size(b0,2) ~= 3 & dist == [3,4,17]) || (size(b0,2) ~= 4 & dist == [6,21])
        
        % Ewa - size,2 te? trzeba tu sprawdza?, czy si? zgadza z liczb? potrzebnych parametrów dla danego rozk?adu
        % i w zaleznosci od tego, czy dopuszczamy obciecie danych (czyli spike w 0) czy nie. 
        % Nie rozumiem tej drugiej czêœci - co trzeba wiêcej sprawdziæ?
        
        cprintf(rgb('DarkOrange'), 'WARNING: Incorrect number of starting values - using midpoint-based estimates as starting values \n')
        b0 = [];
    end
end

midpoint = mean(INPUT.bounds,2);
midpoint(~isfinite(midpoint)) = INPUT.bounds(~isfinite(midpoint),2);
midpoint(~isfinite(midpoint)) = INPUT.bounds(~isfinite(midpoint),1);

if ~exist('b0','var') || isempty(b0)
    OptimOptFit = statset('gevfit');
    OptimOptFit.MaxFunEvals = 1e12;
    OptimOptFit.MaxIter = 1e6;
    warning('off','stats:gevfit:ConvergedToBoundary')
    warning('off','stats:mlecov:NonPosDefHessian')
%     OptimOptFit.Display = 'iter';

    switch dist

% unbounded
        case 0 % normal
            pd = fitdist(midpoint,'Normal','Options',OptimOptFit);
            b0 = pd.ParameterValues; % mu, sigma
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
        case 6 % Johnson SU        
            pd = f_johnson_fit(midpoint);
            b0 =  pd.coef; % gamma delta xi lambda
% stable
        
% bounded (0,Inf)        
        case 10 % exponential
            pd = fitdist(midpoint,'Exponential','Options',OptimOptFit);
            b0 = pd.ParameterValues; % mu        
        case 11 % lognormal
            pd = fitdist(midpoint,'Lognormal','Options',OptimOptFit);
            b0 = pd.ParameterValues; % mu, sigma
        case 12 % loglogistic
            pd = fitdist(midpoint,'Loglogistic','Options',OptimOptFit);
            b0 = pd.ParameterValues; % mu, sigma
        case 13 % Weibull
            pd = fitdist(midpoint,'Weibull','Options',OptimOptFit);
            b0 = pd.ParameterValues; % A, B
        case 14 % Rayleigh 
            pd = fitdist(midpoint,'Rayleigh','Options',OptimOptFit);
            b0 = pd.ParameterValues; % B
        case 15 % Gamma
            pd = fitdist(midpoint,'Gamma','Options',OptimOptFit);
            b0 = pd.ParameterValues; % a, b
        case 16 % BirnbaumSaunders
            pd = fitdist(midpoint,'BirnbaumSaunders','Options',OptimOptFit);
            b0 = pd.ParameterValues; % beta, gamma 
        case 17 % Generalized Pareto
            pd = fitdist(midpoint,'GeneralizedPareto','Options',OptimOptFit);
            b0 = pd.ParameterValues; % k, sigma, theta
        case 18 % InverseGaussian
            pd = fitdist(midpoint,'InverseGaussian','Options',OptimOptFit);
            b0 = pd.ParameterValues; % k, sigma
        case 19 % Nakagami
            pd = fitdist(midpoint,'Nakagami','Options',OptimOptFit);
            b0 = pd.ParameterValues; % mu, omega
        case 20 % Rician
            pd = fitdist(midpoint,'Rician','Options',OptimOptFit);
            b0 = pd.ParameterValues; % s, sigma
        case 21 % Johnson SB
            pd = f_johnson_fit(midpoint);
            b0 = pd.coef; % gamma delta xi lambda
    %     case x % Burr
    %         pd = fitdist(midpoint,'Burr','Options',OptimOptFit); % Error - Parto distribution fits better
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
            pd = fitdist(round(midpoint),'NegativeBinomial','Options',OptimOptFit);
            b0 = pd.ParameterValues; % R, P    
    end
end

% check RealMin
% if ~exist('EstimOpt','var') || isempty(EstimOpt) || ~isfield(EstimOpt,'RealMin') || isempty(EstimOpt.RealMin)
%     EstimOpt.RealMin = true;
% end
% if ~islogical(EstimOpt.RealMin)
%     if any(EstimOpt.RealMin == [0,1])
%         EstimOpt.RealMin = (EstimOpt.RealMin==1);
%     else
%         error('EstimOpt.RealMin value not logical')
%     end
% end    
    
% check optimizer options:
if ~isfield(EstimOpt,'OptimOpt') || isempty(EstimOpt.OptimOpt)
    EstimOpt.OptimOpt = optimoptions('fminunc');
    EstimOpt.OptimOpt.Algorithm = 'quasi-newton';
    EstimOpt.OptimOpt.MaxFunEvals = 1e6; % Maximum number of function evaluations allowed (1000)
    EstimOpt.OptimOpt.MaxIter = 1e3; % Maximum number of iterations allowed (500)    
    EstimOpt.OptimOpt.GradObj = 'off';
    EstimOpt.OptimOpt.FinDiffType = 'central'; % ('forward')
    EstimOpt.OptimOpt.TolFun = 1e-12;
    EstimOpt.OptimOpt.TolX = 1e-12;
    EstimOpt.OptimOpt.OutputFcn = @outputf;
%     OptimOpt.Display = 'iter-detailed';
end


[WTP.beta, WTP.fval, WTP.flag, WTP.out, WTP.grad, WTP.hess] = fminunc(@(b) -sum(LL_DistFit(INPUT.bounds,dist,b)), b0, EstimOpt.OptimOpt);
