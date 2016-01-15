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
elseif ~any(dist == [0:18])
	error('Unsupported distribution type')
end

% starting values: % Wczeœniej by³o tu wiêcej opcji, ale zawsze b0 ma
% rozmiar 1, wiêc to zmieni³am
if exist('b0','var') && ~isempty(b0)
    b0 = b0(:);
    if size(b0,1) ~= 1
        cprintf(rgb('DarkOrange'), 'WARNING: Incorrect number of starting values - using midpoint-based estimates as starting values \n')
        b0 = [];
    end
end
if ~exist('b0','var') || isempty(b0)
%     WTPm = INPUT.bounds(isfinite(INPUT.bounds(:,1)),1); 
    if dist == 0 % normal
        [b0(1), b0(2)] = normfit(midpoint);
    elseif  dist == 1 % lognormal
        b0 = lognfit(midpoint);
    elseif dist == 2 % exponential
        b0 = mean(midpoint);
    elseif dist == 3 % logistic
        b0 = [mean(midpoint), std(midpoint)]; 
    elseif dist == 4 % loglogistic
        b0 = [log(mean(midpoint)), log(std(midpoint))];
    elseif dist == 5 % Poisson
        b0 = poissfit(midpoint);
    elseif dist == 6 % Weibull
        b0 = wblfit(midpoint);
    elseif dist == 7 % Gamma
        b0 = gamfit(midpoint);
    elseif dist == 8 % Rayleigh 
        b0 = raylfit(midpoint);
    elseif dist == 9 % BirnbaumSaunders
        result = fitdist(midpoint,'BirnbaumSaunders');
        b0 = [result.beta, result.gamma];
    elseif dist == 10 % Extreme Value
        b0 = evfit(midpoint);
    elseif dist == 11 % Generalized Extreme Value
        b0 = gevfit(midpoint);
    elseif dist == 12 % Generalized Pareto
        b0 = gpfit(midpoint);
    elseif dist == 13 % InverseGaussian
        result = fitdist(midpoint,'InverseGaussian');
        b0 = [result.mu, result.lambda];
    elseif dist == 14 % Nakagami
        result = fitdist(midpoint,'Nakagami');
        b0 = [result.mu, result.omega];
    elseif dist==15 % Rician
        result = fitdist(midpoint,'Rician');
        b0 = [result.s, result.sigma];
    elseif dist==16 % tLocationScale
        result = fitdist(midpoint,'tLocationScale');
        b0 = [result.mu, result.sigma, result.nu];
    elseif dist==17 % Johnson SU
        result = f_johnson_fit(midpoint); 
        b0 =  result.coef;
    elseif dist==18 % Johnson SB
        result = f_johnson_fit(midpoint);
        b0 =  result.coef;  
    elseif dist==19 % Burr
        b0 = fitdist(midpoint,'Burr'); % Error - Parto distribution fits better
        % Podaje oszacowania parametrów dla Pareto, ale jak je wykorzystaæ,
        % inaczej ni¿ kopiuj¹c? Nie s¹ to te same oszacowania, co z funkcji
        % gpfit.
    elseif dist==20 % F
        b0 = 
    elseif dist==21 % negative binomial
        b0 = nbinfit(round(midpoint));
    elseif dist==22 % uniform
        b0 = [min(midpoint) max(midpoint)];
    elseif dist==23 % generalized inverse Gaussian
        b0 = 
    elseif dist==24 % sinh-arcsinh
        b0 = 
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
