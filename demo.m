clear all; clc;
DATA = load('data.mat');

DATA.WTPfpU(DATA.WTPfpU == 401) = Inf;
DATA.WTPhpU(DATA.WTPhpU == 401) = Inf;

% dist = ...
%     0; % Normal
%     1; % Logistic
%     2; % Extreme Value
%     3; % Generalized Extreme Value
%     4; % tLocationScale
%     5; % Uniform
%     6; % Johnson SU

%     10; % Expotential
%     11; % Lognormal
%     12; % Loglogistic
%     13; % Weibull
%     14; % Rayleigh
%     15; % Gamma
%     16; % BirnbaumSaunders
%     17; % Generalized Pareto
%     18; % Inverse Gaussian
%     19; % Nakagami
%     20; % Rician
%     21; % Johnson SB

%     31; % Poisson
%     32; % Negative Binomial
    
INPUT.bounds = [DATA.WTPfpL DATA.WTPfpU];


% WTP = WTPfit(INPUT);

% Te nazwy poni¿ej do zmiany, bo zmieni³a siê kolejnoœæ. Ale nie jestem
% pewna, jak poprawiæ, bo numery dist nie s¹ kolejne, a to chyba wp³ywa na
% wykorzystanie wektora nazw do pêtli poni¿ej.
names = {'normal' 'lognormal' 'exponential' 'logistic' 'loglogistic' 'Poisson' 'Weibull' 'Gamma' 'Rayleigh' 'BirnbaumSaunders' 'ExtremeValue' 'GeneralizedExtremeValue' 'GeneralizedPareto' 'InverseGaussian' 'Nakagami' 'Rician' 'tLocationScale' 'JohnsonSU' 'JohnsonSB' 'Burr' 'F' 'NegativeBinomial' 'uniform' 'GeneralizedInverseGaussian' 'sinh-arcsinh'};

% Jeœli chodzi³o o cell matrix, to czy coœ na kszta³t tego:
% A = [0 1 2 3 4 5 6]
% B = {'Normal' 'Logistic' 'ExtremeValue' 'GeneralizedExtremeValue' 'tLocationScale' 'Uniform' 'Johnson SU'}
% V = {A B}
% Wtedy jak to zastosowaæ to do pêtli poni¿ej?

% for i = 1:cellfun('length',V(1));
%     WTP.(A{i}) = DistFit(INPUT,i);
% end

% for i = 1:cellfun('length',V(1));
%     LL(i,1:2) = {A{i}, WTP.(A{i}).fval};
% end

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

