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
%         INPUT.NamesX = num2cell(INPUT.NamesX);
    else
        INPUT.NamesX = INPUT.NamesX(:);
    end
end

if sum(INPUT.WT) ~= size(INPUT.bounds,1)
	cprintf(rgb('DarkOrange'), 'WARNING: Scaling weights for unit mean \n')
    INPUT.WT = INPUT.WT - mean(INPUT.WT) + 1;
end

numDistParam = 1*any(dist == [10,14,31]) + 2*any(dist == [0:2,5,11:13,15,16,18:20,32]) + 3*any(dist == [3,4,17]) + 4*any(dist == [6,21,22]);
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
%             pd = f_johnson_fit(midpoint);
%             b0 = pd.coef; % gamma, delta, mi, sigma
            b0 = [1,1,1,1];
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
        case 21 % Johnson SB
            b0 = [1,1,1,1];
        case 22 % Johnson SL
            b0 = [1,1,1,1];

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
    6  'Johnson_SU'; ...

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
    21  'Johnson_SB'; ...
    22  'Johnson_SL'; ...

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


%% generate output


% save tmpout1

WTP.beta = WTP.beta(:);
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
WTP.std(logical(imag(WTP.std))) = NaN;
WTP.pv = pv(WTP.beta,WTP.std);
WTP.stars = cell(size(WTP.beta));
WTP.stars(WTP.pv < 0.01) = {'***'};
WTP.stars((WTP.pv >= 0.01) & (WTP.pv < 0.05)) = {'** '};
WTP.stars((WTP.pv >= 0.05) & (WTP.pv < 0.1)) = {'*  '};
WTP.stars(WTP.pv >= 0.1) = {'   '};
WTP.stars(isnan(WTP.pv)) = {'   '};

betaX = WTP.beta(numDistParam+INPUT.SpikeTrue+1:end);
betaX = num2cell(reshape(betaX,numX,numDistParam+INPUT.SpikeTrue));
stdX = WTP.std(numDistParam+INPUT.SpikeTrue+1:end);
stdX = num2cell(reshape(stdX,numX,numDistParam+INPUT.SpikeTrue));
pvX = WTP.pv(numDistParam+INPUT.SpikeTrue+1:end);
pvX = num2cell(reshape(pvX,numX,numDistParam+INPUT.SpikeTrue));
starsX = WTP.stars(numDistParam+INPUT.SpikeTrue+1:end);
starsX = reshape(starsX,numX,numDistParam+INPUT.SpikeTrue);

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
    R_out(4,1+4*numDistParam+1:1+4*numDistParam+4) = [WTP.beta(numDistParam+1),WTP.stars(numDistParam+1),WTP.std(numDistParam+1),WTP.pv(numDistParam+1)];
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

% tu powinny by� jeszcze wysymulowane: 
% �rednia rozk�adu, s.d. rozk�adu, 95% przedzia� ufno�ci i prawdopodobie�stwo spike 
% (wszystko razem z s.e. i gwiazdkami)

R_out(numX+6,1) = {'Model characteristics:'};
R_out(numX+7:numX+10,1) = {'LL';'AIC/n';'n';'k'};
R_out(numX+7:numX+10,2) = num2cell([WTP.fval; ((2*numB-2*WTP.fval) + 2*numB*(numB+1)/(size(INPUT.bounds,1)-numB-1))/size(INPUT.bounds,1); size(INPUT.bounds,1); numB]);

WTP.R_out = R_out;


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


