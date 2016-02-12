function WTP = DistFit(INPUT, Spike, varargin)


%% check no. of inputs

if nargin < 2 
	error('Too few input arguments for DistFit(INPUT,Spike,dist,b0,EstimOpt)')
elseif nargin == 2
    dist = [];
    b0 = [];
    EstimOpt = [];    
elseif nargin == 3
    dist = varargin{1};
    b0 = [];
    EstimOpt = [];
elseif nargin == 4
    dist = varargin{1};
    b0 = varargin{2};
    EstimOpt = [];
elseif nargin == 5
	dist = varargin{1};
    b0 = varargin{2};
    EstimOpt = varargin{3};
end

% save DistFit_tmp
% return

%% distribution type:

if ~exist('dist','var') || isempty(dist)
    disp('Assuming normally distributed WTP')
    dist = 0;
elseif ~any(dist == [0:7,10:21,31:32])
	error('Unsupported distribution type')
end

%% check INPUT

if size(INPUT.bounds,2) < 2
    error('Providing lower and upper bounds for WTP required')
elseif any(isnan(INPUT.bounds(:)))
    error('Some of the lower or upper bounds for WTP missing (NaN)')
elseif any(sum(isfinite(INPUT.bounds),2)==0) % both bounds infinite
    error('Dataset includes observations with infinite bounds')
elseif any(INPUT.bounds(:,1) > INPUT.bounds(:,2))
    error('Some lower bounds are greater than upper bounds')
elseif any(dist == [10:21,31:32]) && any(INPUT.bounds(:,2) < 0)
    error('Negative upper bounds not consistent with the distribution type')
elseif any(dist == [10:21,31:32]) && any(INPUT.bounds(:,1) < 0)
    cprintf(rgb('DarkOrange'), 'WARNING: Negative lower bounds not consistent with the distribution type - censoring to 0 \n')
    INPUT.bounds(INPUT.bounds(:,1)<0,1) = 0;
    INPUT.bounds(:,2) = max(INPUT.bounds,[],2);
elseif dist == 5 && any(INPUT.bounds(:,1) == -Inf)
    cprintf(rgb('DarkOrange'), 'WARNING: -Inf lower bounds not consistent with uniform distribution - censoring to minimum bound \n')
    INPUT.bounds(INPUT.bounds(:,1)==-Inf,1) = min([INPUT.bounds(isfinite(INPUT.bounds));0]);
    INPUT.bounds(:,2) = max(INPUT.bounds,[],2);
end
if any(dist == 11:20) && any(INPUT.bounds(:,2) == 0)
    cprintf(rgb('DarkOrange'), 'WARNING: 0 upper bounds not consistent with lognormal distribution - censoring to eps \n')
    INPUT.bounds(INPUT.bounds(:,2)==0,2) = eps;
end

% perhaps we will need to check if bounds == 0, but let's leave if for later, when the spike option is added. 
if any(dist == 31:32) && isinteger(INPUT.bounds(isfinite(INPUT.bounds)))
    cprintf(rgb('DarkOrange'), 'WARNING: Rounding up finite bounds to integer values to match the distribution type \n')
    INPUT.bounds(isfinite(INPUT.bounds)) = round(INPUT.bounds(isfinite(INPUT.bounds)));
end

% this is temporary:
if dist == 21 && any(INPUT.bounds(:) == 0)
    cprintf(rgb('DarkOrange'), 'WARNING: Bounds = 0 not consistent with the distribution type - censoring to eps \n')
    INPUT.bounds(INPUT.bounds==0) = eps;    
end
if dist == 21 && any(isinf(INPUT.bounds(:)))
    cprintf(rgb('DarkOrange'), 'WARNING: Bounds = Inf not consistent with the distribution type - censoring to max bound * 2 \n')
    INPUT.bounds(isinf(INPUT.bounds)) = max(INPUT.bounds(~isinf(INPUT.bounds(:,2)),2)) * 2;
end



%% starting values:

if exist('b0','var') && ~isempty(b0)
    b0 = b0(:);
    k = 2*any(dist == [10,14,31]) + 3*any(dist == [0:2,5,11:13,15,16,18:20,32]) + 4*any(dist == [3,4,17]) + 5*any(dist == [6,21]);    
    if size(b0,2) ~= k
       cprintf(rgb('DarkOrange'), 'WARNING: Incorrect number of starting values - using midpoint-based estimates as starting values \n')
       b0 = [];
    end

        % Ewa - size,2 te? trzeba tu sprawdza?, czy si? zgadza z liczb? potrzebnych parametrów dla danego rozk?adu
        % i w zaleznosci od tego, czy dopuszczamy obciecie danych (czyli spike w 0) czy nie. 
        
        % Nie rozumiem tej drugiej czêœci - co trzeba wiêcej sprawdziæ?
        
        % Jak b?dzie model ze spike i podanymi warto?ciami startowymi, to trzeba b?dzie sprawdza? czy ich liczba jest ok. 
        % modele ze spike b?d? mia?y wi?cej parametrów. 
        % Potem trzeba te? b?dzie sprawdza?, czy modele ze zmiennymi obja?niaj?cymi maj? wystarczaj?c? liczb? parametrów. 
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
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % mu, sigma
        case 1 % logistic
            pd = fitdist(midpoint,'Logistic','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % mu, sigma
        case 2 % Extreme Value
            pd = fitdist(midpoint,'ExtremeValue','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % mu sigma
        case 3 % Generalized Extreme Value
            pd = fitdist(midpoint,'GeneralizedExtremeValue','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % k, sigma, mu
        case 4 % tLocationScale
            pd = fitdist(midpoint,'tLocationScale','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % mu, sigma, nu
        case 5 % uniform
            b0 = [min([INPUT.bounds(isfinite(INPUT.bounds));0]) max(INPUT.bounds(isfinite(INPUT.bounds))) Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10];            
        case 6 % Johnson SU        
            pd = f_johnson_fit(midpoint);
            b0 =  [pd.coef Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % gamma delta xi lambda
% stable
        
% bounded (0,Inf)        
        case 10 % exponential
            pd = fitdist(midpoint,'Exponential','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % mu        
        case 11 % lognormal
            midpoint(midpoint <= 0) = eps;
            pd = fitdist(midpoint,'Lognormal','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % mu, sigma
        case 12 % loglogistic
            pd = fitdist(midpoint,'Loglogistic','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % mu, sigma
        case 13 % Weibull
            pd = fitdist(midpoint,'Weibull','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % A, B
        case 14 % Rayleigh 
            pd = fitdist(midpoint,'Rayleigh','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % B
        case 15 % Gamma
            pd = fitdist(midpoint,'Gamma','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % a, b
        case 16 % BirnbaumSaunders
            pd = fitdist(midpoint,'BirnbaumSaunders','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % beta, gamma 
        case 17 % Generalized Pareto
            pd = fitdist(midpoint,'GeneralizedPareto','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % k, sigma, theta
        case 18 % InverseGaussian
            pd = fitdist(midpoint,'InverseGaussian','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % k, sigma
        case 19 % Nakagami
            pd = fitdist(midpoint,'Nakagami','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % mu, omega
        case 20 % Rician
            pd = fitdist(midpoint,'Rician','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % s, sigma
%         case 21 % Johnson SB
%             pd = f_johnson_fit(midpoint);
%             b0 = [pd.coef Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % gamma delta xi lambda
%             save tmp1
%             b0(3) = min(b0(3),min(INPUT.bounds(:))-eps);
%             if any((INPUT.bounds(:)-b0(3))./b0(4) >= 1)
%                 b0(4) = 1 ./ (max((INPUT.bounds(:)-b0(3))) + eps);
%             end

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
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % lambda
        case 32 % negative binomial
            pd = fitdist(round(midpoint),'NegativeBinomial','Options',OptimOptFit);
            b0 = [pd.ParameterValues Spike*sum(INPUT.bounds(INPUT.bounds(:,1)==INPUT.bounds(:,2),1)==0)/size(INPUT.bounds,1)+(1-Spike)*10]; % R, P    
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
