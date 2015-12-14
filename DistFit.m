function WTP = DistFit(INPUT, varargin)

% if nargin < 1 % check no. of inputs
% 	error('Too few input arguments for WTPfit(INPUT,dist,b0,OptimOpt)')
% end

if nargin == 1
    dist = [];
    b0 = [];
    OptimOpt = [];    
elseif nargin == 2
    dist = varargin{1};
    b0 = [];
    OptimOpt = [];
elseif nargin == 3
    dist = varargin{1};
    b0 = varargin{2};
    OptimOpt = [];
elseif nargin == 4
	dist = varargin{1};
    b0 = varargin{2};
    OptimOpt = varargin{3};
end
   
% save tmp1

% INPUT
if size(INPUT.bounds,2) < 2
    error('Providing lower and upper bounds for WTP required')
elseif any(isnan(INPUT.bounds(:)))
    error('Some of the lower or upper bounds for WTP missing')
elseif any(INPUT.bounds(:,1) == INPUT.bounds(:,2))
 	cprintf(rgb('DarkOrange'), 'WARNING: Some of the upper and lower bounds equal - shifting the offending upper bounds by eps \n')
    INPUT.bounds(INPUT.bounds(:,1) == INPUT.bounds(:,2),2) = INPUT.bounds(INPUT.bounds(:,1) == INPUT.bounds(:,2),2) + eps;
end
    
% distribution type:
if ~exist('dist','var') || isempty(dist)
    disp('Assuming normally distributed WTP')
    dist = 0;
elseif ~any(dist == [0:8])
	error('Unsupported distribution type')
end

% starting values:
if exist('b0','var') && ~isempty(b0)
    b0 = b0(:);
    if size(b0,1) ~= 2
        cprintf(rgb('DarkOrange'), 'WARNING: Incorrect number of starting values - using midpoint-based estimates as starting values \n')
        b0 = [];
    end
end
if ~exist('b0','var') || isempty(b0)
%     WTPm = INPUT.bounds(isfinite(INPUT.bounds(:,1)),1); 
    if dist == 0 % normal
        [b0(1), b0(2)] = normfit(INPUT.bounds(isfinite(INPUT.bounds(:,1)),1));
    elseif  dist == 1 % lognormal
        Idontknow = INPUT.bounds(:,1);
        Idontknow(Idontknow==0) = eps;
        b0 = lognfit(Idontknow(isfinite(Idontknow)));
    elseif dist == 2 % exponential
        b0 = [mean(INPUT.bounds(:,1),'omitnan')];
    elseif dist == 3 % logistic
        b0 = [mean(INPUT.bounds(:,1),'omitnan'), std(INPUT.bounds(:,1),'omitnan')];
    elseif dist == 4 % loglogistic
        b0 = [log(mean(INPUT.bounds(:,1),'omitnan')), log(std(INPUT.bounds(:,1),'omitnan'))];
    elseif dist == 5 % Poisson
        b0 = poissfit(INPUT.bounds(isfinite(INPUT.bounds(:,1)),1));
    elseif dist == 6 % Weibull
        Idontknow = INPUT.bounds(:,1);
        Idontknow(Idontknow==0) = eps;
        b0 = wblfit(Idontknow(isfinite(Idontknow)));
    elseif dist == 7 % Gamma
        b0 = gamfit(INPUT.bounds(isfinite(INPUT.bounds(:,1)),1));
    elseif dist == 8 % Rayleigh 
        b0 = raylfit(INPUT.bounds(isfinite(INPUT.bounds(:,1)),1));
    end
end

% optimizer options:
if ~exist('OptimOpt','var') || isempty(OptimOpt)
    OptimOpt = optimoptions('fminunc');
    OptimOpt.Algorithm = 'quasi-newton';
    OptimOpt.MaxFunEvals = 1e6; % Maximum number of function evaluations allowed (1000)
    OptimOpt.MaxIter = 1e3; % Maximum number of iterations allowed (500)    
    OptimOpt.GradObj = 'off';
    OptimOpt.FinDiffType = 'central'; % ('forward')
    OptimOpt.TolFun = 1e-12;
    OptimOpt.TolX = 1e-12;
    OptimOpt.OutputFcn = @outputf;
%     OptimOpt.Display = 'iter-detailed';
end


[WTP.beta, WTP.fval, WTP.flag, WTP.out, WTP.grad, WTP.hess] = fminunc(@(b) -sum(LL_DistFit(INPUT.bounds,dist,b)), b0, OptimOpt);
