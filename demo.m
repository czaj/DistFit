clear all; clc;
DATA = load('data.mat');

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
    
INPUT.bounds = [DATA.WTPfpL DATA.WTPfpU];

% WTP = WTPfit(INPUT);

names = {'normal' 'lognormal' 'exponential' 'logistic' 'loglogistic' 'Poisson' 'Weibull' 'Gamma' 'Rayleigh'};

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

