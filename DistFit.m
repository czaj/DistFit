function WTP = DistFit(INPUT,varargin)


%% check no. of inputs

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
    disp('Assuming normally distributed WTP')
    dist = 0;
elseif ~any(dist == [0:6,10:21,31:32])
	error('Unsupported distribution type')
end

if ~isfield(INPUT,'SpikeTrue') || isempty(INPUT.SpikeTrue)
    INPUT.SpikeTrue = false;
end
if ~islogical(INPUT.SpikeTrue)
    if any(INPUT.SpikeTrue == [0,1])
        INPUT.SpikeTrue = (INPUT.SpikeTrue==1);
    else
        error('INPUT.SpikeTrue value not logical')
    end
end
if ~isfield(INPUT,'HessEstFix') || isempty(INPUT.HessEstFix)
    INPUT.HessEstFix = 3;
elseif any(INPUT.HessEstFix == 0:3)
    error('Incorrect INPUT.HessEstFix setting - available options are: 0 - retained from optimization, 1 - BHHH based, 2 - numerical high-precision Jacobian based, 3 - numerical Hessian based')
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

if ~isfield(INPUT,'WT') | isempty(INPUT.WT)
     INPUT.WT = ones([size(INPUT.X,1),1]);
end

if ~isempty(INPUT.WT) && (size(INPUT.WT,1) ~= size(INPUT.X,1))
    error('Number of weights not consistent with the number of observations')
end
if ~isempty(INPUT.WT) && (size(INPUT.WT,2) ~= 1)
    error('Matrix of weights is not a single column vector')
end

if ~isfield(INPUT,'X')
    INPUT.X = zeros(size(INPUT.bounds,1),0);
elseif size(INPUT.X,1) ~= size(INPUT.bounds,1)
    error('Inconsistent size of bounds and explanatory variables (X)')
else
    NaN_index = sum(isnan(INPUT.X),2) > 0;
    if sum(NaN_index)>0
        cprintf(rgb('DarkOrange'), 'WARNING: Skipping %d out of %d observatons with missing values of the explanatory variables \n',sum(NaN_index),size(NaN_index,1)) % Ta czêœæ po ró¿owym komentarzu jest potrzebna?
        INPUT.bounds = INPUT.bounds(~NaN_index,:);
        INPUT.X = INPUT.X(~NaN_index,:);
        INPUT.WT = INPUT.WT(~NaN_index,:);
    end
    if isfield(INPUT,'NamesX') == 0 || isempty(INPUT.NamesX) || length(INPUT.NamesX) ~= size(INPUT.X,2)
        INPUT.NamesX = (1:size(INPUT.X,2))';
        INPUT.NamesX = cellstr(num2str(INPUT.NamesX));
%         INPUT.NamesX = num2cell(INPUT.NamesX);
    else
        INPUT.NamesX = INPUT.NamesX(:);
    end
end

numDistParam = 1*any(dist == [10,14,31]) + 2*any(dist == [0:2,5,11:13,15,16,18:20,32]) + 3*any(dist == [3,4,17]) + 4*any(dist == [6,21]);
numX = size(INPUT.X,2);
numB = (numDistParam + INPUT.SpikeTrue) * (1 + numX);


%% starting values:

if exist('b0','var') && ~isempty(b0)
    b0 = b0(:)';
    if size(b0,1) ~= numB
       cprintf(rgb('DarkOrange'), 'WARNING: Incorrect number of starting values - using midpoint-based estimates as starting values \n')
       b0 = [];
    end
end

if INPUT.SpikeTrue
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
            b0 = [min([bounds_tmp(isfinite(bounds_tmp));0]) max(bounds_tmp(isfinite(bounds_tmp)))];
        case 6 % Johnson SU        
            pd = f_johnson_fit(midpoint);
            b0 = pd.coef; % gamma, delta, mi, sigma
% stable
        
% bounded (0,Inf)        
        case 10 % exponential
            pd = fitdist(midpoint,'Exponential','Options',OptimOptFit);
            b0 = pd.ParameterValues; % mu
        case 11 % lognormal
            midpoint(midpoint <= 0) = eps;
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
%        case 21 % Johnson SB

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
    
    b0 = [b0 pSpike zeros(1,numX*(numDistParam + INPUT.SpikeTrue))];
    
end
    
% check optimizer options:
if ~isfield(INPUT,'OptimOpt') || isempty(INPUT.OptimOpt)
    INPUT.OptimOpt = optimoptions('fminunc');
    INPUT.OptimOpt.Algorithm = 'quasi-newton';
    INPUT.OptimOpt.MaxFunEvals = 1e6; % Maximum number of function evaluations allowed (1000)
    INPUT.OptimOpt.MaxIter = 1e3; % Maximum number of iterations allowed (500)    
    INPUT.OptimOpt.GradObj = 'off';
    INPUT.OptimOpt.FinDiffType = 'central'; % ('forward')
    INPUT.OptimOpt.TolFun = 1e-12;
    INPUT.OptimOpt.TolX = 1e-12;
    INPUT.OptimOpt.OutputFcn = @outputf;
%     OptimOpt.Display = 'iter-detailed';
end

Distributions = {...
    0  'Normal'; ...
    1  'Logistic'; ...
    2  'Extreme_Value'; ...
    3  'Generalized_Extreme_Value'; ...
    4  'tLocationScale'; ...
    5  'Uniform'; ...
%    6  'Johnson_SU'; ...

    10  'Expotential'; ...
    11  'Lognormal'; ...
    12  'Loglogistic'; ...
    13  'Weibull'; ...
    14  'Rayleigh'; ...
    15  'Gamma'; ...
    16  'BirnbaumSaunders'; ...
    17  'Generalized_Pareto'; ...
    18  'Inverse_Gaussian'; ...
    19  'Nakagami'; ...
    20  'Rician'; ...
%     21  'Johnson_SB'; ...

    31  'Poisson'; ...
    32  'Negative_Binomial'};

cprintf('Fitting ')
cprintf('*blue',[Distributions{[Distributions{:,1}]==dist,2},' '])
cprintf('distribution')
if INPUT.SpikeTrue
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


[WTP.beta, WTP.fval, WTP.flag, WTP.out, WTP.grad, WTP.hess] = fminunc(@(b) -sum(LL_DistFit(INPUT.bounds,INPUT.X,INPUT.WT, dist,INPUT.SpikeTrue,b)), b0, INPUT.OptimOpt);

% save tmp1
% return

WTP.beta = WTP.beta';
WTP.fval = -WTP.fval;
if INPUT.HessEstFix == 0
    WTP.ihess = inv(WTP.hess);
elseif INPUT.HessEstFix == 1
    WTP.f = LL_DistFit(INPUT.bounds,INPUT.X,INPUT.WT,dist,INPUT.SpikeTrue,WTP.beta);
    WTP.jacobian1 = numdiff(@(B) LL_DistFit(INPUT.bounds,INPUT.X,INPUT.WT,dist,INPUT.SpikeTrue,B),WTP.f,WTP.beta,isequal(INPUT.OptimOpt.FinDiffType,'central'));
    WTP.hess1 = WTP.jacobian1'*WTP.jacobian1;
    WTP.ihess = inv(WTP.hess1);
elseif INPUT.HessEstFix == 2
	WTP.jacobian2 = jacobianest(@(B) LL_DistFit(INPUT.bounds,INPUT.X,INPUT.WT,dist,INPUT.SpikeTrue,B),WTP.beta);
    WTP.hess2 = WTP.jacobian2'*WTP.jacobian2;
    WTP.ihess = inv(WTP.hess2);
elseif INPUT.HessEstFix == 3
	WTP.hess3 = hessian(@(B) -sum(LL_DistFit(INPUT.bounds,INPUT.X,INPUT.WT,dist,INPUT.SpikeTrue,B)),WTP.beta);
    WTP.ihess = inv(WTP.hess3);
end

WTP.std = sqrt(diag(WTP.ihess));
WTP.std(logical(imag(WTP.std)))= NaN;
WTP.pv = pv(WTP.beta,WTP.std);
WTP.stars = cell(size(WTP.beta));
WTP.stars(WTP.pv < 0.01) = {'***'};
WTP.stars((WTP.pv >= 0.01) & (WTP.pv < 0.05)) = {'**'};
WTP.stars((WTP.pv >= 0.05) & (WTP.pv < 0.1)) = {'*'};

betaX = WTP.beta(numDistParam+INPUT.SpikeTrue+1:end);
betaX = num2cell(reshape(betaX,numX,numDistParam+INPUT.SpikeTrue));
stdX = WTP.std(numDistParam+INPUT.SpikeTrue+1:end);
stdX = num2cell(reshape(stdX,numX,numDistParam+INPUT.SpikeTrue));
pvX = WTP.pv(numDistParam+INPUT.SpikeTrue+1:end);
pvX = num2cell(reshape(pvX,numX,numDistParam+INPUT.SpikeTrue));
starsX = WTP.stars(numDistParam+INPUT.SpikeTrue+1:end);
starsX = reshape(starsX,numX,numDistParam+INPUT.SpikeTrue);

R_out(1,1) = {[Distributions{[Distributions{:,1}]==dist,2}]};
R_out(4,1) = {'dist. parameters:'};
head = {'var.','coef.','sig.','s.e.','p-value'};
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
%         %     case x % Burr
%         %         pd = fitdist(midpoint,'Burr','Options',OptimOptFit); % Error - Parto distribution fits better
%         %         b0 = pd.ParameterValues; %
%         %     elseif dist==23 % generalized inverse Gaussian
%         % %         b0 =
%         %     elseif dist==24 % sinh-arcsinh
%         % %         b0 =
      
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
    R_out(4,2+(i-1)*4:2+(i-1)*4+3) = [WTP.beta(i),WTP.stars(i),WTP.std(i),WTP.pv(i)];
end

if INPUT.SpikeTrue
    R_out(2,1+4*numDistParam+1) = {'spike'};
    R_out(3,1+4*numDistParam+1:1+4*numDistParam+4) = head(2:5);
end

if numX > 0
    R_out(5:numX+4,1) = INPUT.NamesX;
    for i = 1:numDistParam
        R_out(5:numX+4,2+(i-1)*4:2+(i-1)*4+3) = [betaX(:,i),starsX(:,i),stdX(:,i),pvX(:,i)];
    end
    if INPUT.SpikeTrue
        R_out(5:numX+4,1+4*numDistParam+1:1+4*numDistParam+4) = [betaX(:,numDistParam+1),starsX(:,numDistParam+1),stdX(:,numDistParam+1),pvX(:,numDistParam+1)];
    end
end

% tu powinny byæ jeszcze wysymulowane: 
% œrednia rozk³adu, s.d. rozk³adu, 95% przedzia³ ufnoœci i prawdopodobieñstwo spike 
% (wszystko razem z s.e. i gwiazdkami)

R_out(numX+6,1) = {'Model characteristics'};
R_out(numX+7:numX+(numX>0)+9,1) = {'LL';'AIC/n';'n';'k'};
R_out(numX+7:numX+(numX>0)+9,2) = num2cell([WTP.fval; ((2*numB-2*WTP.fval) + 2*numB*(numB+1)/(size(INPUT.bounds,1)-numB-1))/size(INPUT.bounds,1); size(INPUT.bounds,1); numB]);

% R_out can be used for display

WTP.R_out = R_out;



% Poczatek do wyswietlania wynikow...

if ~isfield(OptimOptFit, 'NVARA') || isempty(OptimOptFit.NVARA)
    OptimOptFit.NVARA = size(INPUT.X,2);
end
if ~isfield(OptimOptFit, 'spacing') || isempty(OptimOptFit.spacing)
    OptimOptFit.spacing = 4;
end
if ~isfield(OptimOptFit, 'precision') || isempty(OptimOptFit.precision)
    OptimOptFit.precision = 4;
end

s = OptimOptFit.spacing;
prec = OptimOptFit.precision;

c1 = [R_out(4:OptimOptFit.NVARA+4,1)]; % Variable names
s1 = max(cellfun(@(x) numel(x), c1)); 

c2 = [R_out(4:OptimOptFit.NVARA+4,2)]; % BDist(1)
c2 = c2(cellfun(@(x) ~isempty(x),c2));
c2 = c2(cellfun(@(x) isfinite(x),c2)); 
s2 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c2))+1,1])+prec+1;
c3 = [R_out(4:OptimOptFit.NVARA+4,3)]; % stars
s3 = max(cellfun(@(x) numel(x), c3));
c4 = [R_out(4:OptimOptFit.NVARA+4,4)]; % s.e.
c4 = c4(cellfun(@(x) ~isempty(x),c4));
c4 = c4(cellfun(@(x) isfinite(x),c4)); 
s4 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c4))+1,1])+prec+1;
c5 = [R_out(4:OptimOptFit.NVARA+4,5)]; % p-value
c5 = c5(cellfun(@(x) ~isempty(x),c5));
c5 = c5(cellfun(@(x) isfinite(x),c5)); 
s5 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c5))+1,1])+prec+1;

if numDistParam > 1
c6 = [R_out(4:OptimOptFit.NVARA+4,6)]; % BDist(2)
c6 = c6(cellfun(@(x) ~isempty(x),c6));
c6 = c6(cellfun(@(x) isfinite(x),c6)); 
s6 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c6))+1,1])+prec+1;
c7 = [R_out(4:OptimOptFit.NVARA+4,7)]; % stars
s7 = max(cellfun(@(x) numel(x), c7));
c8 = [R_out(4:OptimOptFit.NVARA+4,8)]; % s.e.
c8 = c8(cellfun(@(x) ~isempty(x),c8));
c8 = c8(cellfun(@(x) isfinite(x),c8)); 
s8 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c8))+1,1])+prec+1;
c9 = [R_out(4:OptimOptFit.NVARA+4,9)]; % p-value
c9 = c9(cellfun(@(x) ~isempty(x),c9));
c9 = c9(cellfun(@(x) isfinite(x),c9)); 
s9 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c9))+1,1])+prec+1;
end

if numDistParam > 2
c10 = [R_out(4:OptimOptFit.NVARA+4,10)]; % BDist(3)
c10 = c10(cellfun(@(x) ~isempty(x),c10));
c10 = c10(cellfun(@(x) isfinite(x),c10)); 
s10 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c10))+1,1])+prec+1;
c11 = [R_out(4:OptimOptFit.NVARA+4,11)]; % stars
s11 = max(cellfun(@(x) numel(x), c11));
c12 = [R_out(4:OptimOptFit.NVARA+4,12)]; % s.e.
c12 = c12(cellfun(@(x) ~isempty(x),c12));
c12 = c12(cellfun(@(x) isfinite(x),c12)); 
s12 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c12))+1,1])+prec+1;
c13 = [R_out(4:OptimOptFit.NVARA+4,13)]; % p-value
c13 = c13(cellfun(@(x) ~isempty(x),c13));
c13 = c13(cellfun(@(x) isfinite(x),c13)); 
s13 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c13))+1,1])+prec+1;
end

if numDistParam > 3
c14 = [R_out(4:OptimOptFit.NVARA+4,14)]; % BDist(4)
c14 = c14(cellfun(@(x) ~isempty(x),c14));
c14 = c14(cellfun(@(x) isfinite(x),c14)); 
s14 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c14))+1,1])+prec+1;
c15 = [R_out(4:OptimOptFit.NVARA+4,15)]; % stars
s15 = max(cellfun(@(x) numel(x), c15));
c16 = [R_out(4:OptimOptFit.NVARA+4,16)]; % s.e.
c16 = c16(cellfun(@(x) ~isempty(x),c16));
c16 = c16(cellfun(@(x) isfinite(x),c16)); 
s16 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c16))+1,1])+prec+1;
c17 = [R_out(4:OptimOptFit.NVARA+4,17)]; % p-value
c17 = c17(cellfun(@(x) ~isempty(x),c17));
c17 = c17(cellfun(@(x) isfinite(x),c17)); 
s17 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), c17))+1,1])+prec+1;
end

if INPUT.SpikeTrue
cS1 = [R_out(4:OptimOptFit.NVARA+4,1+4*numDistParam+1)]; % Spike
cS1 = cS1(cellfun(@(x) ~isempty(x),cS1));
cS1 = cS1(cellfun(@(x) isfinite(x),cS1)); 
sS1 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), cS1))+1,1])+prec+1;
cS2 = [R_out(4:OptimOptFit.NVARA+4,1+4*numDistParam+2)]; % stars
sS2 = max(cellfun(@(x) numel(x), cS2));
cS3 = [R_out(4:OptimOptFit.NVARA+4,1+4*numDistParam+3)]; % s.e.
cS3 = cS3(cellfun(@(x) ~isempty(x),cS3));
cS3 = cS3(cellfun(@(x) isfinite(x),cS3)); 
sS3 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), cS3))+1,1])+prec+1;
cS4 = [R_out(4:OptimOptFit.NVARA+4,1+4*numDistParam+4)]; % p-value
cS4 = cS4(cellfun(@(x) ~isempty(x),cS4));
cS4 = cS4(cellfun(@(x) isfinite(x),cS4)); 
sS4 = max([max(cellfun(@(x) floor(log10(max(abs(x)))), cS4))+1,1])+prec+1;
end



% Jesli to, co powyzej ok, to tu sie pojawia dalsze komendy zwiazane z
% fprintf, ktora na razie wyglada dosc magicznie...

% fprintf('\n%*s\n\n', (s1+s2+3+s4+prec+2+s*5 + s6+1+s8+3+s10+prec+2+s*4 + size(R_out{1,1},2))/2,R_out{1,1})
% fprintf('%-*s%-s\n', s1+s2+3+s4+prec+2+s*5, R_out{2,2}, R_out{2,6})
% fprintf('%-*s%*s%*s%*s%s%-*s%*s%*s%*s%*s\n', s1,R_out{3,1}, s2+s,R_out{3,2}, s4+3+s,R_out{3,4}, prec+2+s,R_out{3,5},blanks(2*s), s6,R_out{3,6}, s7,R_out{3,7}, s8+s,R_out{3,8}, s10+3+s,R_out{3,10}, s+prec+2,R_out{3,11})
% fprintf('%-*s%-*s%*.0f%*.*f%-3s%*.*f%*.*f\n', s1+s2+3+s4+prec+2+s*5,R_out{4,1}, s6,R_out{4,6}, s7,R_out{4,7}, s8+s,prec,R_out{4,8}, R_out{4,9}, s10+s,prec,R_out{4,10}, 2+prec+s,prec,R_out{4,11})
% for i = 1:OptimOptFit.NVARA
%     fprintf('%-*s%*.*f%-3s%*.*f%*.*f%s%-*s%*.0f%*.*f%-3s%*.*f%*.*f\n', s1,R_out{4+i,1}, s2+s,prec, R_out{4+i,2}, R_out{4+i,3}, s4+s,prec,R_out{4+i,4}, 2+prec+s,prec,R_out{4+i,5}, blanks(2*s),s6,R_out{4+i,6}, s7,R_out{4+i,7}, s8+s,prec,R_out{4+i,8}, R_out{4+i,9}, s10+s,prec,R_out{4+i,10}, 2+prec+s,prec,R_out{4+i,11})
% end
% fprintf('\n%-s\n', R_out{OptimOptFit.NVARA+5,2})
% fprintf('%-*s%*s%*s%*s\n', s1,R_out{OptimOptFit.NVARA+6,1}, s2+s,R_out{OptimOptFit.NVARA+6,2}, s4+3+s,R_out{OptimOptFit.NVARA+6,4}, prec+2+s,R_out{OptimOptFit.NVARA+6,5})
% fprintf('%-*s%*.*f%-3s%*.*f%*.*f\n', s1,R_out{OptimOptFit.NVARA+7,1}, s2+s,prec, R_out{OptimOptFit.NVARA+7,2}, R_out{OptimOptFit.NVARA+7,3}, s4+s,prec,R_out{OptimOptFit.NVARA+7,4}, 2+prec+s,prec,R_out{OptimOptFit.NVARA+7,5})

