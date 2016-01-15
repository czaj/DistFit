clear all; clc;
DATA = load('dane.mat');

DATA.WTPfpU(DATA.WTPfpU == 401) = Inf;
DATA.WTPhpU(DATA.WTPhpU == 401) = Inf;

% dist = ...
%     0; % normal
%     1; % lognormal
%     2; % exponential
%     3; % logistic
%     4; % loglogistic
%     5; % Poisson
%     6; % Weibull
%     7; % Gamma
%     8; % Rayleigh
%     9; % BirnbaumSaunders
%     10; % Extreme Value
%     11; % Generalized Extreme Value
%     12; % Generalized Pareto
%     13; % Inverse Gaussian
%     14; % Nakagami
%     15; % Rician
%     16; % tLocationScale
%     17; % Johnson SU
%     18; % Johnson SB
%     19; % Burr
%     20; % F 
%     21; % negative binomial
%     22; % uniform
%     23; % generalized inverse Gaussian
%     24; % sinh-arcsinh
    
INPUT.bounds = [DATA.WTPfpL DATA.WTPfpU];
midpoint = zeros(size(INPUT.bounds,1),1);

for ii = 1:length(midpoint)
    if isfinite(DATA.WTPfpU(ii))==1
       midpoint(ii) = (DATA.WTPfpL(ii) + DATA.WTPfpU(ii)) / 2;
    else
       midpoint(ii) = DATA.WTPfpL(ii);
    end
end

% WTP = WTPfit(INPUT);

names = {'normal' 'lognormal' 'exponential' 'logistic' 'loglogistic' 'Poisson' 'Weibull' 'Gamma' 'Rayleigh' 'BirnbaumSaunders' 'ExtremeValue' 'GeneralizedExtremeValue' 'GeneralizedPareto' 'InverseGaussian' 'Nakagami' 'Rician' 'tLocationScale' 'JohnsonSU' 'JohnsonSB' 'Burr' 'F' 'NegativeBinomial' 'uniform' 'GeneralizedInverseGaussian' 'sinh-arcsinh'};

for i = 1:size(names,2);
%     if i == 6 || i == 9; continue; end
    WTP.(names{i}) = DistFit(INPUT,i-1);
end

for i = 1:size(names,2);
%     if i == 6 || i == 9; continue; end    
    LL(i,1:2) = {names{i}, WTP.(names{i}).fval};
end


% for i = 1:8
%     eval(sprintf('WTP.dist%d = WTPfit(INPUT,i)', i));
% end
% 
% for i = 6:9
%     LL(i,1) = eval(sprintf('WTP.dist%d .fval', i))
% end

