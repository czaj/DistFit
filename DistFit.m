function Results = DistFit(INPUT,varargin)


%% check no. of inputs
cprintf('\n\n')
if nargin < 1
    error('Too few input arguments for DistFit(INPUT,dist,b0,EstimOpt)')
elseif nargin == 1
    dist = [];
    b0 = [];
elseif nargin == 2
    dist = varargin{1};
    b0 = [];
elseif nargin == 3
    dist = varargin{1};
    b0 = varargin{2};
end

% save DistFit_tmp
% return

%% distribution type:

if ~exist('dist','var') || isempty(dist)
    disp('Assuming normally distributed Results')
    dist = 0;
elseif ~any(dist == [0:6,10:22,31:32])
    error('Unsupported distribution type')
end

if ~isfield(INPUT,'Spike') || isempty(INPUT.Spike)
    INPUT.Spike = false;
end
if ~islogical(INPUT.Spike)
    if any(INPUT.Spike == [0,1])
        INPUT.Spike = (INPUT.Spike==1);
    else
        error('INPUT.Spike value not logical')
    end
end
if ~isfield(INPUT,'SimStats') || isempty(INPUT.SimStats)
    INPUT.SimStats = false;
end
if ~islogical(INPUT.SimStats)
    if any(INPUT.SimStats == [0,1])
        INPUT.SimStats = (INPUT.SimStats==1);
    else
        error('INPUT.SimStats value not logical')
    end
end
if ~isfield(INPUT,'HessEstFix') || isempty(INPUT.HessEstFix)
    INPUT.HessEstFix = 2;
elseif all(INPUT.HessEstFix ~= 0:3)
    error('Incorrect INPUT.HessEstFix setting - available options are: 0 - retained from optimization, 1 - BHHH based, 2 - numerical high-precision Jacobian based, 3 - numerical Hessian based')
end

%% check INPUT

if size(INPUT.bounds,2) < 2
    error('Providing lower and upper bounds for Results required')
elseif any(isnan(INPUT.bounds(:)))
    error('Some of the lower or upper bounds for Results missing (NaN)')
elseif any(sum(isfinite(INPUT.bounds),2)==0) % both bounds infinite
    error('Dataset includes observations with infinite bounds')
elseif any(INPUT.bounds(:,1) > INPUT.bounds(:,2))
    error('Some lower bounds are greater than upper bounds')
    % elseif any(dist == [10:21,31:32]) && any(INPUT.bounds(:,2) < 0)
    %     error('Negative upper bounds not consistent with the distribution type')
    % elseif any(dist == [10:20,31:32]) && any(INPUT.bounds(:,1) < 0)
    %     cprintf(rgb('DarkOrange'), 'WARNING: Negative lower bounds not consistent with the distribution type - censoring to 0 \n')
    %     INPUT.bounds(INPUT.bounds(:,1)<0,1) = 0;
    %     INPUT.bounds(:,2) = max(INPUT.bounds,[],2);
    % elseif dist == 5 && any(INPUT.bounds(:,1) == -Inf)
    %     cprintf(rgb('DarkOrange'), 'WARNING: -Inf lower bounds not consistent with uniform distribution - censoring to minimum finite bound \n')
    %     INPUT.bounds(INPUT.bounds(:,1)==-Inf,1) = min([INPUT.bounds(isfinite(INPUT.bounds));0]);
    %     INPUT.bounds(:,2) = max(INPUT.bounds,[],2);
end


% Czy mog?aby? sprawdzi?, czy osobno trzeba rozpatrywa? przypadek ujemnych warto?ci i niedodatnich warto?ci - a je?li tak to rozdzieli? te 2 przypadki?
% Spike troch? to zmienia, bo jak jest spike w zerze to nawet rozk?ady, które wymagaj? dodatnich warto?ci a która? jest równa 0 to i tak jest ok.

if any(dist == [10:21,31:32]) && (any(INPUT.bounds(:,2) <= 0) || (INPUT.Spike && any(INPUT.bounds(:,2) < 0)))
    cprintf(rgb('DarkOrange'), 'WARNING: Data contains observations not consistent with the distribution type (negative) \n')
end

% if any(dist == [31:32]) && ~isinteger(INPUT.bounds(isfinite(INPUT.bounds))) % also 32?
% %     cprintf(rgb('DarkOrange'), 'WARNING: Rounding up finite bounds to integer values to match the distribution type \n')
% %     INPUT.bounds(isfinite(INPUT.bounds)) = round(INPUT.bounds(isfinite(INPUT.bounds)));
%     cprintf(rgb('DarkOrange'), 'WARNING: Data contains observations not consistent with the distribution type (not integer)\n')
% end

% if any(dist == 11:20) && any(INPUT.bounds(:,2) == 0)
%     cprintf(rgb('DarkOrange'), 'WARNING: 0 upper bounds not consistent with lognormal distribution - censoring to eps \n')
%     INPUT.bounds(INPUT.bounds(:,2)==0,2) = eps;
% end
% if any(dist == [21:22]) && any(INPUT.bounds(:,1) == -Inf)
%     cprintf(rgb('DarkOrange'), 'WARNING: -Inf lower bounds not consistent with Johnson SB / SL distribution - censoring to minimum finite bound or 0 \n')
%     INPUT.bounds(INPUT.bounds(:,1)==-Inf,1) = min([INPUT.bounds(isfinite(INPUT.bounds));0]);
%     INPUT.bounds(:,2) = max(INPUT.bounds,[],2);
% end
if dist == 21 && any(INPUT.bounds(:,2) == Inf)
    cprintf(rgb('DarkOrange'), 'WARNING: Inf upper bounds not consistent with Johnson SB distribution - censoring to maximum finite bound *10\n')
    INPUT.bounds(INPUT.bounds(:,2)==Inf,2) = max(INPUT.bounds(isfinite(INPUT.bounds)))*10;
    INPUT.bounds(:,1) = min(INPUT.bounds,[],2);
end


% this is temporary:
% if dist == 21 && any(INPUT.bounds(:) == 0)
% cprintf(rgb('DarkOrange'), 'WARNING: Bounds = 0 not consistent with the distribution type - censoring to eps \n')
% INPUT.bounds(INPUT.bounds==0) = eps;
% end
% if dist == 21 && any(isinf(INPUT.bounds(:)))
% cprintf(rgb('DarkOrange'), 'WARNING: Bounds = Inf not consistent with the distribution type - censoring to max bound\n')
% INPUT.bounds(isinf(INPUT.bounds)) = max(INPUT.bounds(~isinf(INPUT.bounds(:,2)),2));
% end

if ~isfield(INPUT,'WT') || isempty(INPUT.WT)
    INPUT.WT = ones(size(INPUT.bounds,1),1);
elseif (size(INPUT.WT,1) ~= size(INPUT.bounds,1))
    error('Number of weights not consistent with the number of observations')
elseif size(INPUT.WT,2) ~= 1
    error('Matrix of weights is not a single column vector')
else
    disp('Observations are weighted')
end

if ~isfield(INPUT,'X')
    INPUT.X = zeros(size(INPUT.bounds,1),0);
elseif size(INPUT.X,1) ~= size(INPUT.bounds,1)
    error('Inconsistent size of bounds and explanatory variables (X)')
else
    NaN_index = sum(isnan(INPUT.X),2) > 0;
    if sum(NaN_index)>0
        cprintf(rgb('DarkOrange'), 'WARNING: Skipping ')
        cprintf(['*',num2str(rgb('DarkOrange'))], [num2str(sum(NaN_index)),' '])
        cprintf(rgb('DarkOrange'), 'out of ')
        cprintf(['*',num2str(rgb('DarkOrange'))], [num2str(size(NaN_index,1)),' '])
        cprintf(rgb('DarkOrange'), 'observatons with missing values of the explanatory variables \n')
        INPUT.bounds = INPUT.bounds(~NaN_index,:);
        INPUT.X = INPUT.X(~NaN_index,:);
        INPUT.WT = INPUT.WT(~NaN_index,:);
    end
    if isfield(INPUT,'NamesX') == 0 || isempty(INPUT.NamesX) || length(INPUT.NamesX) ~= size(INPUT.X,2)
        INPUT.NamesX = (1:size(INPUT.X,2))';
        INPUT.NamesX = cellstr(num2str(INPUT.NamesX));
        % INPUT.NamesX = num2cell(INPUT.NamesX);
    else
        INPUT.NamesX = INPUT.NamesX(:);
    end
end

if sum(INPUT.WT) ~= size(INPUT.bounds,1)
    cprintf(rgb('DarkOrange'), 'WARNING: Scaling weights for unit mean \n')
    INPUT.WT = INPUT.WT - mean(INPUT.WT) + 1;
end

numDistParam = 1*any(dist == [10,14,31]) + 2*any(dist == [0:2,5,11:13,15,16,18:21,32]) + 3*any(dist == [3,4,17,22]) + 4*any(dist == [6]);
numX = size(INPUT.X,2);
numB = (numDistParam + INPUT.Spike) * (1 + numX);


%% starting values:


if exist('b0','var') && ~isempty(b0)
    b0 = b0(:);
    if size(b0,1) ~= numB
        cprintf(rgb('DarkOrange'), 'WARNING: Incorrect number of starting values - using midpoint-based estimates as starting values \n')
        b0 = [];
    end
end

if INPUT.Spike
    bounds_tmp = INPUT.bounds(INPUT.bounds(:,1)~=INPUT.bounds(:,2)~=0,:);
    pSpike = norminv(max(1 - size(bounds_tmp,1)./size(INPUT.bounds,1),0.1),0,1);
else
    bounds_tmp = INPUT.bounds;
    pSpike = [];
end
midpoint = mean(bounds_tmp,2);
midpoint(~isfinite(midpoint)) = bounds_tmp(~isfinite(midpoint),2);
midpoint(~isfinite(midpoint)) = bounds_tmp(~isfinite(midpoint),1);

if ~exist('b0','var') || isempty(b0)
    OptimOptFit = statset('gevfit');
    OptimOptFit.MaxFunEvals = 1e12;
    OptimOptFit.MaxIter = 1e6;
    warning('off','stats:gevfit:ConvergedToBoundary')
    warning('off','stats:mlecov:NonPosDefHessian')
    % OptimOptFit.Display = 'iter';
    
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
            b0 = [min([bounds_tmp(isfinite(bounds_tmp));0]) max(bounds_tmp(isfinite(bounds_tmp)))];
        case 6 % Johnson SU
            % pd = f_johnson_fit(midpoint);
            % b0 = pd.coef;
            b0 = [0,1,mean(midpoint),std(midpoint)/1.78]; % gamma, delta, mi, sigma
            % stable
            
            % bounded (0,Inf)
        case 10 % exponential
            pd = fitdist(midpoint(midpoint>0),'Exponential','Options',OptimOptFit);
            b0 = pd.ParameterValues; % mu
        case 11 % lognormal
            pd = fitdist(midpoint(midpoint>0),'Lognormal','Options',OptimOptFit);
            b0 = pd.ParameterValues; % mu, sigma
        case 12 % loglogistic
            pd = fitdist(midpoint(midpoint>0),'Loglogistic','Options',OptimOptFit);
            b0 = pd.ParameterValues; % mu, sigma
        case 13 % Weibull
            pd = fitdist(midpoint(midpoint>0),'Weibull','Options',OptimOptFit);
            b0 = pd.ParameterValues; % A, B
        case 14 % Rayleigh
            pd = fitdist(midpoint(midpoint>0),'Rayleigh','Options',OptimOptFit);
            b0 = pd.ParameterValues; % B
        case 15 % Gamma
            pd = fitdist(midpoint(midpoint>0),'Gamma','Options',OptimOptFit);
            b0 = pd.ParameterValues; % a, b
        case 16 % BirnbaumSaunders
            pd = fitdist(midpoint(midpoint>0),'BirnbaumSaunders','Options',OptimOptFit);
            b0 = pd.ParameterValues; % beta, gamma
        case 17 % Generalized Pareto
            pd = fitdist(midpoint(midpoint>0),'GeneralizedPareto','Options',OptimOptFit);
            b0 = pd.ParameterValues; % k, sigma, theta
        case 18 % InverseGaussian
            pd = fitdist(midpoint(midpoint>0),'InverseGaussian','Options',OptimOptFit);
            b0 = pd.ParameterValues; % k, sigma
        case 19 % Nakagami
            pd = fitdist(midpoint(midpoint>0),'Nakagami','Options',OptimOptFit);
            b0 = pd.ParameterValues; % mu, omega
        case 20 % Rician
            pd = fitdist(midpoint(midpoint>0),'Rician','Options',OptimOptFit);
            b0 = pd.ParameterValues; % s, sigma
        case 21 % Johnson SB % gamma, delta, mi, sigma
            % b0 = [skewness(midpoint,0),kurtosis(midpoint,0),mean(midpoint),std(midpoint)];
            % pd = f_johnson_fit(midpoint);
            % b0 = pd.coef;
            b0 = [0,1,min(INPUT.bounds(:))-eps,max(INPUT.bounds(:)) - min(INPUT.bounds(:))+eps]; % gamma, delta, mi, sigma
        case 22 % Johnson SL
            b0 = [0,1,min(INPUT.bounds(:))-eps,max(INPUT.bounds(:)) - min(INPUT.bounds(:))+eps]; % gamma, delta, mi, sigma
            
            % case x % Burr
            % pd = fitdist(midpoint,'Burr','Options',OptimOptFit); % Error - Parto distribution fits better
            % b0 = pd.ParameterValues; %
            % elseif dist==23 % generalized inverse Gaussian
            % % b0 =
            % elseif dist==24 % sinh-arcsinh
            % % b0 =
            
            % discrete
        case 31 % Poisson
            pd = fitdist(round(midpoint),'Poisson','Options',OptimOptFit);
            b0 = pd.ParameterValues; % lambda
        case 32 % negative binomial
            pd = fitdist(round(midpoint),'NegativeBinomial','Options',OptimOptFit);
            b0 = pd.ParameterValues; % R, P
    end
    
    b0 = [b0 pSpike zeros(1,numX*(numDistParam + INPUT.Spike))];
    
end


%% check optimizer options:
% if ~isfield(INPUT,'OptimOpt') || isempty(INPUT.OptimOpt)
%     INPUT.OptimOpt = optimoptions('fminunc');
%     INPUT.OptimOpt.Algorithm = 'quasi-newton';
%     INPUT.OptimOpt.MaxFunEvals = 1e6; % Maximum number of function evaluations allowed (1000)
%     INPUT.OptimOpt.MaxIter = 1e3; % Maximum number of iterations allowed (500)
%     INPUT.OptimOpt.GradObj = 'off';
%     INPUT.OptimOpt.FinDiffType = 'central'; % ('forward')
%     INPUT.OptimOpt.TolFun = 1e-12;
%     INPUT.OptimOpt.TolX = 1e-12;
%     INPUT.OptimOpt.OutputFcn = @outputf;
%     % OptimOpt.Display = 'iter-detailed';
% end

if ~isfield(INPUT,'OptimOpt') || isempty(INPUT.OptimOpt)
    INPUT.OptimOpt = optimoptions('fmincon');
    INPUT.OptimOpt.Algorithm = 'interior-point'; %'sqp'; 'active-set'; 'trust-region-reflective';
    INPUT.OptimOpt.MaxFunEvals = 1e6; % Maximum number of function evaluations allowed (1000)
    INPUT.OptimOpt.MaxIter = 1e3; % Maximum number of iterations allowed (500)
    INPUT.OptimOpt.GradObj = 'off';
    INPUT.OptimOpt.FinDiffType = 'central'; % ('forward')
    INPUT.OptimOpt.TolFun = 1e-12;
    INPUT.OptimOpt.TolX = 1e-12;
    INPUT.OptimOpt.OutputFcn = @outputf;
    % OptimOpt.Display = 'iter-detailed';
end


%% Print model specification


Distributions = {...
    0 'Normal'; ...
    1 'Logistic'; ...
    2 'Extreme_Value'; ...
    3 'Generalized_Extreme_Value'; ...
    4 'tLocationScale'; ...
    5 'Uniform'; ...
    6 'Johnson_SU'; ...
    
    10 'Expotential'; ...
    11 'Lognormal'; ...
    12 'Loglogistic'; ...
    13 'Weibull'; ...
    14 'Rayleigh'; ...
    15 'Gamma'; ...
    16 'BirnbaumSaunders'; ...
    17 'Generalized_Pareto'; ...
    18 'Inverse_Gaussian'; ...
    19 'Nakagami'; ...
    20 'Rician'; ...
    21 'Johnson_SB'; ...
    22 'Johnson_SL'; ...
    
    31 'Poisson'; ...
    32 'Negative_Binomial'};

cprintf('Fitting ')
cprintf('*blue',[Distributions{[Distributions{:,1}]==dist,2},' '])
cprintf('distribution')
if INPUT.Spike
    if numX>0
        cprintf(' with ')
        cprintf('*blue','spike ')
        cprintf('density at 0 and ')
        cprintf('*blue',[num2str(numX),' '])
        cprintf('covariates of distribution parameters.\n\n')
    else
        cprintf(' with ')
        cprintf('*blue','spike ')
        cprintf('density at 0.\n\n')
    end
else
    if numX>0
        cprintf(' with ')
        cprintf('*blue',[num2str(numX),' '])
        cprintf('covariates of distribution parameters.\n\n')
    else
        cprintf('.\n\n')
    end
end


%% Specify constraints


switch dist
    
    % unbounded
    case {0,1,2,11} % Normal: mu, sigma>=0; Logistic: mu, sigma>=0; % Extreme Value: mu, sigma>=0; Lognormal: mu, sigma>0
        A = -[zeros(size(INPUT.bounds,1),1), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX), INPUT.X, zeros(size(INPUT.bounds,1),numX*INPUT.Spike)];
        C = zeros(size(INPUT.bounds,1),1);
    case 3 % Generalized Extreme Value % k, sigma>0, mu
        A = -[zeros(size(INPUT.bounds,1),1), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX), INPUT.X, zeros(size(INPUT.bounds,1),numX), zeros(size(INPUT.bounds,1),numX*INPUT.Spike)];
        C = zeros(size(INPUT.bounds,1),1);
    case 4 % tLocationScale % mu, sigma>0, nu>0
        A = [zeros(size(INPUT.bounds,1),1), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX), INPUT.X, zeros(size(INPUT.bounds,1),numX), zeros(size(INPUT.bounds,1),numX*INPUT.Spike)];
        A = -[A; [zeros(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),1), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX), zeros(size(INPUT.bounds,1),numX), INPUT.X, zeros(size(INPUT.bounds,1),numX*INPUT.Spike)]];
        C = zeros(2*size(INPUT.bounds,1),1);
    case 5 % Uniform % min, max
        A = [ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            INPUT.X, zeros(size(INPUT.bounds,1),numX), zeros(size(INPUT.bounds,1),numX*INPUT.Spike)];
        A = [A; -[zeros(size(INPUT.bounds,1),1), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX), INPUT.X, zeros(size(INPUT.bounds,1),numX*INPUT.Spike)]];
        C = [INPUT.bounds(:,1); ...
            INPUT.bounds(:,2)];
    case 6 % Johnson SU: gamma, delta>0, mi, sigma>0; 
        A = [zeros(size(INPUT.bounds,1),1), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),2), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX), INPUT.X, zeros(size(INPUT.bounds,1),numX*2), zeros(size(INPUT.bounds,1),numX*INPUT.Spike)];
        A = -[A; [zeros(size(INPUT.bounds,1),3), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX*3), INPUT.X, zeros(size(INPUT.bounds,1),numX*INPUT.Spike)]];
        C = zeros(2*size(INPUT.bounds,1),1);
    case {10,14,31} % Exponential: mu>0; Rayleigh: b0>0; Poisson: lambda>0
        A = -[ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            INPUT.X, zeros(size(INPUT.bounds,1),numX*INPUT.Spike)];
        C = zeros(size(INPUT.bounds,1),1);
    case {12,13,15,16,18,20} % Loglogistic: mu>0, sigma>0; Weibull: A>0, b0>0; Gamma: a>0, b>0; BirnbaumSaunders: beta>0, gamma>0; InverseGaussian: k>0, sigma>0; Rician: s>=0, sigma>0
        A = [ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            INPUT.X, zeros(size(INPUT.bounds,1),numX), zeros(size(INPUT.bounds,1),numX*INPUT.Spike)];
        A = -[A; [zeros(size(INPUT.bounds,1),1), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX), INPUT.X, zeros(size(INPUT.bounds,1),numX*INPUT.Spike)]];
        C = zeros(2*size(INPUT.bounds,1),1);
    case 17 % Generalized Pareto: k, sigma>0, theta
        A = -[zeros(size(INPUT.bounds,1),1), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX), INPUT.X, zeros(size(INPUT.bounds,1),numX), zeros(size(INPUT.bounds,1),numX*INPUT.Spike)];
        C = zeros(size(INPUT.bounds,1),1);
    case 19 % Nakagami: mu>=0.5, omega>0
        A = [ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            INPUT.X, zeros(size(INPUT.bounds,1),numX), zeros(size(INPUT.bounds,1),numX*INPUT.Spike)];
        A = -[A; [zeros(size(INPUT.bounds,1),1), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX), INPUT.X, zeros(size(INPUT.bounds,1),numX*INPUT.Spike)]];
        C = [-0.5*ones(size(INPUT.bounds,1),1); zeros(size(INPUT.bounds,1),1)];
    case 21 % Johnson SB: gamma, delta>0, mi<lbound, x-mi<sigma
        A = -[zeros(size(INPUT.bounds,1),1), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),2), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX), INPUT.X, zeros(size(INPUT.bounds,1),numX*2), zeros(size(INPUT.bounds,1),numX*INPUT.Spike)];
        A = [A; [zeros(size(INPUT.bounds,1),2), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX*2), INPUT.X, zeros(size(INPUT.bounds,1),numX), zeros(size(INPUT.bounds,1),numX*INPUT.Spike)]];
        A = [A; [zeros(size(INPUT.bounds,1),3), -ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX*3), -INPUT.X, zeros(size(INPUT.bounds,1),numX*INPUT.Spike)]];
        C = [zeros(size(INPUT.bounds,1),1);...
            INPUT.bounds(:,1);...
            -INPUT.bounds(:,2)+INPUT.bounds(:,1)];
    case 22 % Johnson SL: gamma, delta>0, mi<lbound, sigma>0
        A = -[zeros(size(INPUT.bounds,1),1), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),2), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX), INPUT.X, zeros(size(INPUT.bounds,1),numX*2), zeros(size(INPUT.bounds,1),numX*INPUT.Spike)];
        A = [A; [zeros(size(INPUT.bounds,1),2), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX*2), INPUT.X, zeros(size(INPUT.bounds,1),numX), zeros(size(INPUT.bounds,1),numX*INPUT.Spike)]];
        A = [A; [zeros(size(INPUT.bounds,1),3), -ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX*3), -INPUT.X, zeros(size(INPUT.bounds,1),numX*INPUT.Spike)]];
        C = [zeros(size(INPUT.bounds,1),1);...
            INPUT.bounds(:,1);...
            zeros(size(INPUT.bounds,1),1)];
    case 32 % Negative Binomial % R>0 and R is an integer, 0<P<1
        A = [ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            INPUT.X, zeros(size(INPUT.bounds,1),numX), zeros(size(INPUT.bounds,1),numX*INPUT.Spike)];
        A = -[A; [zeros(size(INPUT.bounds,1),1), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX), INPUT.X, zeros(size(INPUT.bounds,1),numX*INPUT.Spike)]];
        A = [A; [zeros(size(INPUT.bounds,1),1), ones(size(INPUT.bounds,1),1), zeros(size(INPUT.bounds,1),INPUT.Spike), ...
            zeros(size(INPUT.bounds,1),numX), INPUT.X, zeros(size(INPUT.bounds,1),numX*INPUT.Spike)]];
        C = [zeros(2*size(INPUT.bounds,1),1);ones(size(INPUT.bounds,1),1)];
end



%% Optimization


% [Results.beta, Results.fval, Results.flag, Results.out, Results.grad, Results.hess] = fminunc(@(b) -sum(LL_DistFit(INPUT.bounds,INPUT.X,INPUT.WT, dist,INPUT.Spike,b)), b0, INPUT.OptimOpt);

[Results.beta, Results.fval, Results.flag, Results.out, Results.lambda, Results.grad, Results.hess] = fmincon(@(b) -sum(LL_DistFit(INPUT.bounds,INPUT.X,INPUT.WT, dist,INPUT.Spike,b)), b0,A,C,[],[],[],[],[],INPUT.OptimOpt);


%% generate output


% save tmpout1

Results.beta = Results.beta(:);
Results.fval = -Results.fval;
if INPUT.HessEstFix == 0
    Results.ihess = inv(Results.hess);
elseif INPUT.HessEstFix == 1
    Results.f = LL_DistFit(INPUT.bounds,INPUT.X,INPUT.WT,dist,INPUT.Spike,Results.beta);
    Results.jacobian1 = numdiff(@(B) LL_DistFit(INPUT.bounds,INPUT.X,INPUT.WT,dist,INPUT.Spike,B),Results.f,Results.beta,isequal(INPUT.OptimOpt.FinDiffType,'central'));
    Results.hess1 = Results.jacobian1'*Results.jacobian1;
    Results.ihess = inv(Results.hess1);
elseif INPUT.HessEstFix == 2
    Results.jacobian2 = jacobianest(@(B) LL_DistFit(INPUT.bounds,INPUT.X,INPUT.WT,dist,INPUT.Spike,B),Results.beta);
    Results.hess2 = Results.jacobian2'*Results.jacobian2;
    Results.ihess = inv(Results.hess2);
elseif INPUT.HessEstFix == 3
    Results.hess3 = hessian(@(B) -sum(LL_DistFit(INPUT.bounds,INPUT.X,INPUT.WT,dist,INPUT.Spike,B)),Results.beta);
    Results.ihess = inv(Results.hess3);
end

[~,err] = cholcov(Results.ihess);
if err ~= 0
    cprintf(rgb('DarkOrange'), 'WARNING: AVC is not positive semi-definite - try a different Hessian approximation method using HessEstFix \n')
end

Results.std = sqrt(diag(Results.ihess));
Results.std(logical(imag(Results.std))) = NaN;
Results.pv = pv(Results.beta,Results.std);
Results.stars = cell(size(Results.beta));
Results.stars(Results.pv < 0.01) = {'***'};
Results.stars((Results.pv >= 0.01) & (Results.pv < 0.05)) = {'** '};
Results.stars((Results.pv >= 0.05) & (Results.pv < 0.1)) = {'* '};
Results.stars(Results.pv >= 0.1) = {' '};
Results.stars(isnan(Results.pv)) = {' '};

betaX = Results.beta(numDistParam+INPUT.Spike+1:end);
betaX = num2cell(reshape(betaX,numX,numDistParam+INPUT.Spike));
stdX = Results.std(numDistParam+INPUT.Spike+1:end);
stdX = num2cell(reshape(stdX,numX,numDistParam+INPUT.Spike));
pvX = Results.pv(numDistParam+INPUT.Spike+1:end);
pvX = num2cell(reshape(pvX,numX,numDistParam+INPUT.Spike));
starsX = Results.stars(numDistParam+INPUT.Spike+1:end);
starsX = reshape(starsX,numX,numDistParam+INPUT.Spike);

R_out(1,1) = {['Fitted ',Distributions{[Distributions{:,1}]==dist,2},' distribution parameters']};
R_out(4,1) = {'dist. parameters:'};
head = {'var.','coef.',[],'s.e.','p-value'};
switch dist
    case 0 % normal
        R_out(2,2) = {'mu'};
        R_out(2,6) = {'sigma'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
    case 1 % logistic % mu, sigma
        R_out(2,2) = {'mu'};
        R_out(2,6) = {'sigma'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
    case 2 % Extreme Value
        R_out(2,2) = {'mu'};
        R_out(2,6) = {'sigma'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
    case 3 % Generalized Extreme Value % k, sigma, mu
        R_out(2,2) = {'k'};
        R_out(2,6) = {'sigma'};
        R_out(2,10) = {'mu'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
        R_out(3,10:13) = head(2:5);
    case 4 % tLocationScale % mu, sigma>0, nu>0
        R_out(2,2) = {'mu'};
        R_out(2,6) = {'sigma'};
        R_out(2,10) = {'nu'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
        R_out(3,10:13) = head(2:5);
    case 5 % uniform
        R_out(2,2) = {'min'};
        R_out(2,6) = {'max'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
    case 6 % Johnson SU % gamma delta xi lambda
        R_out(2,2) = {'gamma'};
        R_out(2,6) = {'delta'};
        R_out(2,10) = {'xi'};
        R_out(2,14) = {'lambda'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
        R_out(3,10:13) = head(2:5);
        R_out(3,14:17) = head(2:5);
    case 10 % exponential % mu
        R_out(2,2) = {'mu'};
        R_out(3,1:5) = head;
    case 11 % lognormal % mu, sigma
        R_out(2,2) = {'mu'};
        R_out(2,6) = {'sigma'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
    case 12 % loglogistic % mu, sigma
        R_out(2,2) = {'mu'};
        R_out(2,6) = {'sigma'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
    case 13 % Weibull % A, B
        R_out(2,2) = {'A'};
        R_out(2,6) = {'B'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
    case 14 % Rayleigh % B
        R_out(2,2) = {'B'};
        R_out(3,1:5) = head;
    case 15 % Gamma % a, b
        R_out(2,2) = {'a'};
        R_out(2,6) = {'b'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
    case 16 % BirnbaumSaunders % beta, gamma
        R_out(2,2) = {'beta'};
        R_out(2,6) = {'gamma'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
    case 17 % Generalized Pareto % k, sigma, theta
        R_out(2,2) = {'k'};
        R_out(2,6) = {'sigma'};
        R_out(2,10) = {'theta'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
        R_out(3,10:13) = head(2:5);
    case 18 % InverseGaussian % k, sigma
        R_out(2,2) = {'k'};
        R_out(2,6) = {'sigma'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
    case 19 % Nakagami % mu, omega
        R_out(2,2) = {'mu'};
        R_out(2,6) = {'omega'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
    case 20 % Rician % s, sigma
        R_out(2,2) = {'s'};
        R_out(2,6) = {'sigma'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
    case 21 % Johnson SB % gamma delta xi lambda
        R_out(2,2) = {'gamma'};
        R_out(2,6) = {'delta'};
        R_out(2,10) = {'xi'};
        R_out(2,14) = {'lambda'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
        R_out(3,10:13) = head(2:5);
        R_out(3,14:17) = head(2:5);
    case 22 % Johnson SL % gamma delta xi lambda
        R_out(2,2) = {'gamma'};
        R_out(2,6) = {'delta'};
        R_out(2,10) = {'xi'};
        R_out(2,14) = {'lambda'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
        R_out(3,10:13) = head(2:5);
        R_out(3,14:17) = head(2:5);
        % % case x % Burr
        % % pd = fitdist(midpoint,'Burr','Options',OptimOptFit); % Error - Parto distribution fits better
        % % b0 = pd.ParameterValues; %
        % % elseif dist==23 % generalized inverse Gaussian
        % % % b0 =
        % % elseif dist==24 % sinh-arcsinh
        % % % b0 =
        
    case 31 % Poisson % lambda
        R_out(2,2) = {'lambda'};
        R_out(3,1:5) = head;
    case 32 % negative binomial % R, P
        R_out(2,2) = {'R'};
        R_out(2,6) = {'P'};
        R_out(3,1:5) = head;
        R_out(3,6:9) = head(2:5);
end

for i = 1:numDistParam
    R_out(4,2+(i-1)*4:2+(i-1)*4+3) = [Results.beta(i),Results.stars(i),Results.std(i),Results.pv(i)];
end

if INPUT.Spike
    R_out(2,1+4*numDistParam+1) = {'spike'};
    R_out(3,1+4*numDistParam+1:1+4*numDistParam+4) = head(2:5);
    R_out(4,1+4*numDistParam+1:1+4*numDistParam+4) = [Results.beta(numDistParam+1),Results.stars(numDistParam+1),Results.std(numDistParam+1),Results.pv(numDistParam+1)];
end

if numX > 0
    R_out(5:numX+4,1) = INPUT.NamesX;
    for i = 1:numDistParam
        R_out(5:numX+4,2+(i-1)*4:2+(i-1)*4+3) = [betaX(:,i),starsX(:,i),stdX(:,i),pvX(:,i)];
    end
    if INPUT.Spike
        R_out(5:numX+4,1+4*numDistParam+1:1+4*numDistParam+4) = [betaX(:,numDistParam+1),starsX(:,numDistParam+1),stdX(:,numDistParam+1),pvX(:,numDistParam+1)];
    end
end

% tu powinny byæ jeszcze wysymulowane:
% œrednia rozk³adu, s.d. rozk³adu, 95% przedzia³ ufnoœci i prawdopodobieñstwo spike
% (wszystko razem z s.e. i gwiazdkami)

R_out(numX+6,1) = {'Model characteristics:'};
R_out(numX+7:numX+10,1) = {'LL';'AIC/n';'n';'k'};
R_out(numX+7:numX+10,2) = num2cell([Results.fval; ((2*numB-2*Results.fval) + 2*numB*(numB+1)/(size(INPUT.bounds,1)-numB-1))/size(INPUT.bounds,1); size(INPUT.bounds,1); numB]);

Results.R_out = R_out;


%% display the results

cprintf('\n\nFitted '); cprintf('*blue',[Distributions{[Distributions{:,1}]==dist,2},' ']); cprintf('distribution parameters:\n')

[~,mCW1] = CellColumnWidth(R_out(4:4+numX,:)); % width and max width of each column
spacing = 2;
spacing2 = 2;
precision = 4;

if size(R_out,2) == 5
    fprintf('%*s%-*s\n',sum(mCW1(1:3))+spacing*2+precision+1,'',sum(mCW1(4:5))+(spacing+precision+1)*2, R_out{2,2})
    fprintf('%-*s%*s%3s%*s%*s\n',mCW1(1)+spacing,R_out{3,1}, mCW1(2)+spacing+precision+1,R_out{3,2}, R_out{3,3}, mCW1(4)+spacing+precision,R_out{3,4}, mCW1(5)+spacing+precision+2,R_out{3,5})
    for i=4:4+numX
        fprintf('%-*s% *.*f%-3s% *.*f% *.*f\n',mCW1(1)+spacing,R_out{i,1}, mCW1(2)+spacing+precision+1,precision,R_out{i,2}, R_out{i,3}, mCW1(4)+spacing+precision+1,precision,R_out{i,4}, mCW1(5)+spacing+precision+1,precision,R_out{i,5})
    end
elseif size(R_out,2) == 9
    fprintf('%*s%-*s%*s%*s%-*s\n',sum(mCW1(1:3))+spacing*2+precision+1,'',sum(mCW1(4:5))+(spacing+precision+1)*2, R_out{2,2}, ...
        spacing2+spacing, '', ...
        sum(mCW1(6:7))+precision+1,'', sum(mCW1(8:9))+(spacing+precision+1)*2, R_out{2,6})
    fprintf('%-*s%*s%3s%*s%*s%s%*s%3s%*s%*s\n',mCW1(1)+spacing,R_out{3,1}, mCW1(2)+spacing+precision+1,R_out{3,2}, R_out{3,3}, mCW1(4)+spacing+precision,R_out{3,4}, mCW1(5)+spacing+precision+2,R_out{3,5},...
        blanks(spacing2),...
        mCW1(6)+spacing+precision+1,R_out{3,6}, R_out{3,7}, mCW1(8)+spacing+precision,R_out{3,8}, mCW1(9)+spacing+precision+2,R_out{3,9})
    for i=4:4+numX
        fprintf('%-*s% *.*f%-3s% *.*f% *.*f%s% *.*f%-3s% *.*f% *.*f\n',mCW1(1)+spacing,R_out{i,1}, mCW1(2)+spacing+precision+1,precision,R_out{i,2}, R_out{i,3}, mCW1(4)+spacing+precision+1,precision,R_out{i,4}, mCW1(5)+spacing+precision+1,precision,R_out{i,5},...
            blanks(spacing2),...
            mCW1(6)+spacing+precision+1,precision,R_out{i,6}, R_out{i,7}, mCW1(8)+spacing+precision+1,precision,R_out{i,8}, mCW1(9)+spacing+precision+1,precision,R_out{i,9})
    end
    
elseif size(R_out,2) == 13
    fprintf('%*s%-*s%*s%*s%-*s%*s%*s%-*s\n',sum(mCW1(1:3))+spacing*2+precision+1,'',sum(mCW1(4:5))+(spacing+precision+1)*2, R_out{2,2}, ...
        spacing2+spacing, '', ...
        sum(mCW1(6:7))+precision+1,'', sum(mCW1(8:9))+(spacing+precision+1)*2, R_out{2,6},...
        spacing2+spacing, '', ...
        sum(mCW1(10:11))+precision+1,'', sum(mCW1(12:13))+(spacing+precision+1)*2, R_out{2,10})
    fprintf('%-*s%*s%3s%*s%*s%s%*s%3s%*s%*s%s%*s%3s%*s%*s\n',mCW1(1)+spacing,R_out{3,1}, mCW1(2)+spacing+precision+1,R_out{3,2}, R_out{3,3}, mCW1(4)+spacing+precision,R_out{3,4}, mCW1(5)+spacing+precision+2,R_out{3,5},...
        blanks(spacing2),...
        mCW1(6)+spacing+precision+1,R_out{3,6}, R_out{3,7}, mCW1(8)+spacing+precision,R_out{3,8}, mCW1(9)+spacing+precision+2,R_out{3,9},...
        blanks(spacing2),...
        mCW1(10)+spacing+precision+1,R_out{3,10}, R_out{3,11}, mCW1(12)+spacing+precision,R_out{3,12}, mCW1(13)+spacing+precision+2,R_out{3,13})
    for i=4:4+numX
        fprintf('%-*s% *.*f%-3s% *.*f% *.*f%s% *.*f%-3s% *.*f% *.*f%s% *.*f%-3s% *.*f% *.*f\n',mCW1(1)+spacing,R_out{i,1}, mCW1(2)+spacing+precision+1,precision,R_out{i,2}, R_out{i,3}, mCW1(4)+spacing+precision+1,precision,R_out{i,4}, mCW1(5)+spacing+precision+1,precision,R_out{i,5},...
            blanks(spacing2),...
            mCW1(6)+spacing+precision+1,precision,R_out{i,6}, R_out{i,7}, mCW1(8)+spacing+precision+1,precision,R_out{i,8}, mCW1(9)+spacing+precision+1,precision,R_out{i,9},...
            blanks(spacing2),...
            mCW1(10)+spacing+precision+1,precision,R_out{i,10}, R_out{i,11}, mCW1(12)+spacing+precision+1,precision,R_out{i,12}, mCW1(13)+spacing+precision+1,precision,R_out{i,13})
    end
elseif size(R_out,2) == 17
    fprintf('%*s%-*s%*s%*s%-*s%*s%*s%-*s%*s%*s%-*s\n',sum(mCW1(1:3))+spacing*2+precision+1,'',sum(mCW1(4:5))+(spacing+precision+1)*2, R_out{2,2}, ...
        spacing2+spacing, '', ...
        sum(mCW1(6:7))+precision+1,'', sum(mCW1(8:9))+(spacing+precision+1)*2, R_out{2,6}, ...
        spacing2+spacing, '', ...
        sum(mCW1(10:11))+precision+1,'', sum(mCW1(12:13))+(spacing+precision+1)*2, R_out{2,10}, ...
        spacing2+spacing, '', ...
        sum(mCW1(14:15))+precision+1,'', sum(mCW1(16:17))+(spacing+precision+1)*2, R_out{2,14})
    fprintf('%-*s%*s%3s%*s%*s%s%*s%3s%*s%*s%s%*s%3s%*s%*s%s%*s%3s%*s%*s\n',mCW1(1)+spacing,R_out{3,1}, mCW1(2)+spacing+precision+1,R_out{3,2}, R_out{3,3}, mCW1(4)+spacing+precision,R_out{3,4}, mCW1(5)+spacing+precision+2,R_out{3,5},...
        blanks(spacing2),...
        mCW1(6)+spacing+precision+1,R_out{3,6}, R_out{3,7}, mCW1(8)+spacing+precision,R_out{3,8}, mCW1(9)+spacing+precision+2,R_out{3,9},...
        blanks(spacing2),...
        mCW1(10)+spacing+precision+1,R_out{3,10}, R_out{3,11}, mCW1(12)+spacing+precision,R_out{3,12}, mCW1(13)+spacing+precision+2,R_out{3,13},...
        blanks(spacing2),...
        mCW1(14)+spacing+precision+1,R_out{3,14}, R_out{3,15}, mCW1(16)+spacing+precision,R_out{3,16}, mCW1(17)+spacing+precision+2,R_out{3,17})
    for i=4:4+numX
        fprintf('%-*s% *.*f%-3s% *.*f% *.*f%s% *.*f%-3s% *.*f% *.*f%s% *.*f%-3s% *.*f% *.*f%s% *.*f%-3s% *.*f% *.*f\n',mCW1(1)+spacing,R_out{i,1}, mCW1(2)+spacing+precision+1,precision,R_out{i,2}, R_out{i,3}, mCW1(4)+spacing+precision+1,precision,R_out{i,4}, mCW1(5)+spacing+precision+1,precision,R_out{i,5},...
            blanks(spacing2),...
            mCW1(6)+spacing+precision+1,precision,R_out{i,6}, R_out{i,7}, mCW1(8)+spacing+precision+1,precision,R_out{i,8}, mCW1(9)+spacing+precision+1,precision,R_out{i,9},...
            blanks(spacing2),...
            mCW1(10)+spacing+precision+1,precision,R_out{i,10}, R_out{i,11}, mCW1(12)+spacing+precision+1,precision,R_out{i,12}, mCW1(13)+spacing+precision+1,precision,R_out{i,13},...
            blanks(spacing2),...
            mCW1(14)+spacing+precision+1,precision,R_out{i,14}, R_out{i,15}, mCW1(16)+spacing+precision+1,precision,R_out{i,16}, mCW1(17)+spacing+precision+1,precision,R_out{i,17})
    end
elseif size(R_out,2) == 21
    fprintf('%*s%-*s%*s%*s%-*s%*s%*s%-*s%*s%*s%-*s%*s%*s%-*s\n',sum(mCW1(1:3))+spacing*2+precision+1,'',sum(mCW1(4:5))+(spacing+precision+1)*2, R_out{2,2}, ...
        spacing2+spacing, '', ...
        sum(mCW1(6:7))+precision+1,'', sum(mCW1(8:9))+(spacing+precision+1)*2, R_out{2,6}, ...
        spacing2+spacing, '', ...
        sum(mCW1(10:11))+precision+1,'', sum(mCW1(12:13))+(spacing+precision+1)*2, R_out{2,10}, ...
        spacing2+spacing, '', ...
        sum(mCW1(14:15))+precision+1,'', sum(mCW1(16:17))+(spacing+precision+1)*2, R_out{2,14}, ...
        spacing2+spacing, '', ...
        sum(mCW1(18:19))+precision+1,'', sum(mCW1(20:21))+(spacing+precision+1)*2, R_out{2,18})
    fprintf('%-*s%*s%3s%*s%*s%s%*s%3s%*s%*s%s%*s%3s%*s%*s%s%*s%3s%*s%*s%s%*s%3s%*s%*s\n',mCW1(1)+spacing,R_out{3,1}, mCW1(2)+spacing+precision+1,R_out{3,2}, R_out{3,3}, mCW1(4)+spacing+precision,R_out{3,4}, mCW1(5)+spacing+precision+2,R_out{3,5},...
        blanks(spacing2),...
        mCW1(6)+spacing+precision+1,R_out{3,6}, R_out{3,7}, mCW1(8)+spacing+precision,R_out{3,8}, mCW1(9)+spacing+precision+2,R_out{3,9},...
        blanks(spacing2),...
        mCW1(10)+spacing+precision+1,R_out{3,10}, R_out{3,11}, mCW1(12)+spacing+precision,R_out{3,12}, mCW1(13)+spacing+precision+2,R_out{3,13},...
        blanks(spacing2),...
        mCW1(14)+spacing+precision+1,R_out{3,14}, R_out{3,15}, mCW1(16)+spacing+precision,R_out{3,16}, mCW1(17)+spacing+precision+2,R_out{3,17},...
        blanks(spacing2),...
        mCW1(18)+spacing+precision+1,R_out{3,18}, R_out{3,19}, mCW1(20)+spacing+precision,R_out{3,20}, mCW1(21)+spacing+precision+2,R_out{3,21})
    for i=4:4+numX
        fprintf('%-*s% *.*f%-3s% *.*f% *.*f%s% *.*f%-3s% *.*f% *.*f%s% *.*f%-3s% *.*f% *.*f%s% *.*f%-3s% *.*f% *.*f%s% *.*f%-3s% *.*f% *.*f\n',mCW1(1)+spacing,R_out{i,1}, mCW1(2)+spacing+precision+1,precision,R_out{i,2}, R_out{i,3}, mCW1(4)+spacing+precision+1,precision,R_out{i,4}, mCW1(5)+spacing+precision+1,precision,R_out{i,5},...
            blanks(spacing2),...
            mCW1(6)+spacing+precision+1,precision,R_out{i,6}, R_out{i,7}, mCW1(8)+spacing+precision+1,precision,R_out{i,8}, mCW1(9)+spacing+precision+1,precision,R_out{i,9},...
            blanks(spacing2),...
            mCW1(10)+spacing+precision+1,precision,R_out{i,10}, R_out{i,11}, mCW1(12)+spacing+precision+1,precision,R_out{i,12}, mCW1(13)+spacing+precision+1,precision,R_out{i,13},...
            blanks(spacing2),...
            mCW1(14)+spacing+precision+1,precision,R_out{i,14}, R_out{i,15}, mCW1(16)+spacing+precision+1,precision,R_out{i,16}, mCW1(17)+spacing+precision+1,precision,R_out{i,17},...
            blanks(spacing2),...
            mCW1(18)+spacing+precision+1,precision,R_out{i,18}, R_out{i,19}, mCW1(20)+spacing+precision+1,precision,R_out{i,20}, mCW1(21)+spacing+precision+1,precision,R_out{i,21})
    end
end

fprintf('\n%s\n',R_out{6+numX,1})
[~,mCW2] = CellColumnWidth(R_out(7+numX:10+numX,1:2)); % width and max width of each column
fprintf('%-*s% *.*f\n',mCW2(1)+spacing+spacing2,R_out{7+numX,1}, mCW2(2)+precision+1,precision,R_out{7+numX,2})
fprintf('%-*s% *.*f\n',mCW2(1)+spacing+spacing2,R_out{8+numX,1}, mCW2(2)+precision+1,precision,R_out{8+numX,2})
fprintf('%-*s% *d\n',mCW2(1)+spacing+spacing2,R_out{9+numX,1}, mCW2(2),R_out{9+numX,2})
fprintf('%-*s% *d\n',mCW2(1)+spacing+spacing2,R_out{10+numX,1}, mCW2(2),R_out{10+numX,2})


% save out2



if INPUT.SimStats
    
    if err ~= 0
        cprintf(rgb('DarkOrange'), 'WARNING: AVC is not positive semi-definite - aborting simulation of descriptive statistics \n')
    else
        
        sim1 = 1e3;
        sim2 = 1e3;
        
        Bi = mvnrnd(Results.beta,Results.ihess,sim1)'; % draw parameter distributions taking uncertainty (standard errors) into account
        
        bDist = Bi(1:numDistParam,:); % numDistParam x sim1
        if INPUT.Spike
            if numX > 0 % Spike and X
                meanX = mean(bsxfun(@times,INPUT.X,INPUT.WT));
                for i = 1:numDistParam
                    %   bDist = bDist + repmat(meanX,[1,numDistParam])*Bi(numDistParam+1+1:numDistParam+1 + numDistParam*numX,:);
                    bDist(i,:) = bDist(i,:) + meanX*Bi(numDistParam+1+1+numX*(i-1):numDistParam+1 + numX*i,:);
                end
                bSpike = Bi(numDistParam+1,:) + meanX*Bi(numDistParam*(1+numX)+1+1:(numDistParam+1)*(1+numX),:);
                pSpike = normcdf(bSpike);
                ru01 = random('unif',0,1,[1,sim1]);
                tSpike = ru01 < pSpike;
            else % Spike only
                bSpike = Bi(numDistParam+1,:);
                pSpike = normcdf(bSpike);
                ru01 = random('unif',0,1,[1,sim1]);
                tSpike = ru01 < pSpike;
            end
        else
            if numX > 0 % X only
                meanX = mean(bsxfun(@times,INPUT.X,INPUT.WT));
                for i = 1:numDistParam
                    bDist(i,:) = bDist(i,:) + meanX*Bi(numDistParam+1+numX*(i-1):numDistParam+numX*i,:);
                end
                % bDist = bDist + repmat(meanX,[1,numDistParam])*Bi(numDistParam+1:numDistParam*(1+numX),:);
            else % baseline distribution only
                % bDist only
            end
            tSpike = false(1,sim1);
        end
        
        Bmtx(sim1,sim2) = 0;
        switch dist
            case 0 % normal
                for i = 1:sim1
                    Bmtx(i,:) = random('normal',bDist(1,i),bDist(2,i),[1,sim2]);
                end
            case 1 % logistic % mu, sigma
                for i = 1:sim1
                    Bmtx(i,:) = random('logistic',bDist(1,i),bDist(2,i),[1,sim2]);
                end
            case 2 % Extreme Value % mu sigma
                for i = 1:sim1
                    Bmtx(i,:) = random('Extreme Value',bDist(1,i),bDist(2,i),[1,sim2]);
                end
            case 3 % Generalized Extreme Value % k, sigma, mu
                for i = 1:sim1
                    Bmtx(i,:) = random('Generalized Extreme Value',bDist(1,i),bDist(2,i),bDist(3,i),[1,sim2]);
                end
            case 4 % tLocationScale % mu, sigma, nu
                for i = 1:sim1
                    Bmtx(i,:) = random('tlocationscale',bDist(1,i),bDist(2,i),bDist(3,i),[1,sim2]);
                end
            case 5 % uniform
                for i = 1:sim1
                    Bmtx(i,:) = random('unif',bDist(1,i),bDist(2,i),[1,sim2]);
                end
            case 6 % Johnson SU % gamma, delta, mi, sigma
                % na to nawet nie mamy funkcji...
            case 10 % exponential % mu
                for i = 1:sim1
                    Bmtx(i,:) = random('exp',bDist(1,i),[1,sim2]);
                end
            case 11 % lognormal % mu, sigma
                for i = 1:sim1
                    Bmtx(i,:) = random('logn',bDist(1,i),bDist(2,i),[1,sim2]);
                end
            case 12 % loglogistic % mu, sigma
                for i = 1:sim1
                    Bmtx(i,:) = random('loglogistic',bDist(1,i),bDist(2,i),[1,sim2]);
                end
            case 13 % Weibull % A, B
                for i = 1:sim1
                    Bmtx(i,:) = random('wbl',bDist(1,i),bDist(2,i),[1,sim2]);
                end
            case 14 % Rayleigh % B
                for i = 1:sim1
                    Bmtx(i,:) = random('rayl',bDist(1,i),[1,sim2]);
                end
            case 15 % Gamma % a, b
                for i = 1:sim1
                    Bmtx(i,:) = random('gam',bDist(1,i),bDist(2,i),[1,sim2]);
                end
            case 16 % BirnbaumSaunders % beta, gamma
                for i = 1:sim1
                    Bmtx(i,:) = random('birnbaumsaunders',bDist(1,i),bDist(2,i),[1,sim2]);
                end
            case 17 % Generalized Pareto % k, sigma, theta
                for i = 1:sim1
                    Bmtx(i,:) = random('gp',bDist(1,i),bDist(2,i),bDist(3,i),[1,sim2]);
                end
            case 18 % InverseGaussian % k, sigma
                for i = 1:sim1
                    Bmtx(i,:) = random('inversegaussian',bDist(1,i),bDist(2,i),[1,sim2]);
                end
            case 19 % Nakagami % mu, omega
                for i = 1:sim1
                    Bmtx(i,:) = random('nakagami',bDist(1,i),bDist(2,i),[1,sim2]);
                end
            case 20 % Rician % s, sigma
                for i = 1:sim1
                    Bmtx(i,:) = random('rician',bDist(1,i),bDist(2,i),[1,sim2]);
                end
            case 21 % Johnson SB % gamma, delta, mi, sigma
                % Brak funkcji
            case 22 % Johnson SL % gamma, delta, mi, sigma
                % Brak funkcji
            case 31 % Poisson % lambda
                for i = 1:sim1
                    Bmtx(i,:) = random('Poisson',bDist(1,i),[1,sim2]);
                end
            case 32 % negative binomial % R, P
                for i = 1:sim1
                    Bmtx(i,:) = random('nbin',bDist(1,i),bDist(2,i),[1,sim2]);
                end
        end
        
        %         save out3;
        Bmtx(repmat(tSpike',[1,sim2])) = 0;
        
        stats1 = [mean(Bmtx(:)) std(Bmtx(:)) median(Bmtx(:)) quantile((Bmtx(:)),0.025) quantile((Bmtx(:)),0.975)];
        stats2 = std([mean(Bmtx,1); std(Bmtx,[],1); median(Bmtx,1); quantile(Bmtx,0.025,1); quantile(Bmtx,0.975,1)],[],2)';
        stats31 = quantile([mean(Bmtx,1); std(Bmtx,[],1); median(Bmtx,1); quantile(Bmtx,0.025,1); quantile(Bmtx,0.975,1)],0.025,2)';
        stats32 = quantile([mean(Bmtx,1); std(Bmtx,[],1); median(Bmtx,1); quantile(Bmtx,0.025,1); quantile(Bmtx,0.975,1)],0.975,2)';
        
        stats = [stats1; stats2; stats31; stats32];
        
        stats_out(1,2:6) = {'mean','s.d.','median','q0.025','q0.975'};
        
        
        if INPUT.Spike
            stats = [stats, [mean(pSpike); std(pSpike); quantile(pSpike,0.025); quantile(pSpike,0.0975)]];
            stats_out(1,7) = {'p.spike'};
        end
        
        stats_out(2:5,1) = {'value','s.e.','lower 95% c.i.','upper 95% c.i.'};
        stats_out(2:5,2:6+INPUT.Spike) = num2cell(stats);
        
        Results.stats = stats;
        Results.stats_out = stats_out;
        
        fprintf('\n\n%s\n','Fitted distribution descriptive statistics:')
        
        [~,mCW3] = CellColumnWidth(stats_out); % width and max width of each column
        if INPUT.Spike
            fprintf('%*s%*s%*s%*s%*s%*s%*s\n',mCW3(1)+spacing2,stats_out{1,1}, mCW3(2)+spacing2+precision,stats_out{1,2}, mCW3(3)+spacing2+precision,stats_out{1,3}, mCW3(4)+spacing2+precision,stats_out{1,4}, mCW3(5)+spacing2+precision,stats_out{1,5}, mCW3(6)+spacing2+precision,stats_out{1,6}, mCW3(7)+spacing2+precision,stats_out{1,7})
            for j = 2:size(stats_out,1)
                fprintf('%-*s% *.*f% *.*f% *.*f% *.*f% *.*f% *.*f\n',mCW3(1)+spacing2,stats_out{j,1}, mCW3(2)+spacing2+precision,precision,stats_out{j,2}, mCW3(3)+spacing2+precision,precision,stats_out{j,3}, mCW3(4)+spacing2+precision,precision,stats_out{j,4}, mCW3(5)+spacing2+precision,precision,stats_out{j,5}, mCW3(6)+spacing2+precision,precision,stats_out{j,6}, mCW3(7)+spacing2+precision,precision,stats_out{j,7})
            end
        else
            fprintf('%*s%*s%*s%*s%*s%*s\n',mCW3(1)+spacing2,stats_out{1,1}, mCW3(2)+spacing2+precision,stats_out{1,2}, mCW3(3)+spacing2+precision,stats_out{1,3}, mCW3(4)+spacing2+precision,stats_out{1,4}, mCW3(5)+spacing2+precision,stats_out{1,5}, mCW3(6)+spacing2+precision,stats_out{1,6})
            for j = 2:size(stats_out,1)
                fprintf('%-*s% *.*f% *.*f% *.*f% *.*f% *.*f\n',mCW3(1)+spacing2,stats_out{j,1}, mCW3(2)+spacing2+precision,precision,stats_out{j,2}, mCW3(3)+spacing2+precision,precision,stats_out{j,3}, mCW3(4)+spacing2+precision,precision,stats_out{j,4}, mCW3(5)+spacing2+precision,precision,stats_out{j,5}, mCW3(6)+spacing2+precision,precision,stats_out{j,6})
            end
        end
    end
end


