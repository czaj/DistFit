% More examples for working with Johnson Curves
% 
% by David L. Jones, Apr-2014
%
% This file is part of the 'JOHNSON CURVE TOOLBOX FOR MATLAB'
% and is released under the BSD 2-clause license.
% 
% The following examples follow those from the 'JohnsonDistribution package for
% R'.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Find Quantiles of an SL Distribution:                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear the workspace:
clc; clear all; close all;                  

% Get parameters for an SL family of distributions:
jsn = f_johnson_M(1,1,3,-1)
% jsn = 
%     coef: [0.076578 1.3975 -0.2229 1]
%     type: 'SL'

% Define percentage points:
p = [0.01 0.05 0.1 0.25 0.5 0.75 0.9 0.95 0.99]';

% Get corresponding Normal variates:
Z = norminv(p);

% Transform Normal to Johnson variates:
Y = f_johnson_z2y(Z,jsn.coef,jsn.type);

% Create a QQ plot:
figure; set(gcf,'color','w');
plot(Z,Y,'ko');
title('QQ plot');
xlabel('Normal quantiles');
ylabel('Johnson SL quantiles');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Simulate an SL Distribution:                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear the workspace:
clc; clear all; close all;                  

% Get parameters for an SL family of distributions:
jsn = f_johnson_M(1,1,3,-1)
% jsn = 
%     coef: [0.076578 1.3975 -0.2229 1]
%     type: 'SL'

% Generate random SL variates:
Y = f_johnson_rnd(jsn.coef,'SL',1000,1);

% Use kernel smoothing to estimate the probability density estimates of Y:
[f,Yi] = ksdensity(Y); 

% Create PDF plot:
figure; set(gcf,'color','w');
plot(Yi,f,'k-');
title({'Estimated PDF of Johnson SL Curve','(mean = 1, sd = 1, skew = 3)'});
xlabel('Y');
ylabel('estimated PDF(Y)');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Find the percentage points for an SL Distribution:              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear the workspace:
clc; clear all; close all;                  

% Find coefficients for an 'SL' distribution with mean = 1, sd = 1, and skewness
% = 3 (set kurtosis = -1 to force 'SL'):
jsn = f_johnson_M(1,1,3,-1)
% jsn = 
%     coef: [0.076578 1.3975 -0.2229 1]
%     type: 'SL'

% What are the percentage points of this distribution for the observed values of
% 1 through 5?
Y = (1:5)';
Z = f_johnson_y2z(Y,jsn.coef,jsn.type);
normcdf(Z)
% ans =
%       0.63975
%       0.88355
%       0.95656
%       0.98168
%        0.9915

% Alternative method:
f_johnson_cdf(Y,jsn.coef,jsn.type)
% ans =
%       0.63975
%       0.88355
%       0.95656
%       0.98168
%        0.9915
