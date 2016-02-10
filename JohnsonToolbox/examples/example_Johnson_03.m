% Examples of using the Johnson System of Distributions (Johnson Curves)
% 
% by David L. Jones, Apr-2014
%
% This file is part of the 'JOHNSON CURVE TOOLBOX FOR MATLAB'
% and is released under the BSD 2-clause license.


% The following examples follow the simulation work presented in Tables 1-3 of
% George & Ramachandran (2011).
% 
% George, F. and K. M. Ramachandran. 2011. Estimation of parameters of Johnson's
%  system of distributions. Journal of Modern Applied Statistical Methods 10(2):
%  494-504.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           SB Family (Bounded distribution via logistic transform):           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;                  % clear the workspace
rng('default');                             % reset for repeatability
n    = 1000*100;                            % sample size
coef = [1 1 10 10];                         % SB distribution parameters
Y    = f_johnson_rnd(coef,'SB',n);          % generate random SB variates
SB_m = f_johnson_fit(Y,'M',2);              % fit using moments
SB_q = f_johnson_fit(Y,'Q',0);              % fit using quantiles
[coef' SB_m.coef' SB_q.coef']               % compare actual vs. estimates
% ans =
%             1       1.0052        1.001
%             1       1.0071       1.0078
%            10       9.9899       9.9836
%            10       10.025       10.029


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       SU Family (Unbounded distribution via hyperbolic sine transform):      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;                  % clear the workspace
rng('default');                             % reset for repeatability
n    = 1000*100;                            % sample size
coef = [0 2 10 10];                         % SU distribution parameters
Y    = f_johnson_rnd(coef,'SU',n);          % generate random SU variates
SU_m = f_johnson_fit(Y,'M',1);              % fit using moments
SU_q = f_johnson_fit(Y,'Q',0);              % fit using quantiles
[coef' SU_m.coef' SU_q.coef']               % compare actual vs. estimates
% ans =
% 
%             0    -0.027979     0.021172
%             2       2.0369       1.9784
%            10       9.8389       10.119
%            10       10.186       9.8625


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         SL Family (Lognormal distribution via exponential transform):        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;                   % clear the workspace
rng('default');                              % reset for repeatability
n     = 1000*100;                            % sample size
coef  = [1 3 0 1];                           % SL distribution parameters
Y     = f_johnson_rnd(coef,'SL',n);          % generate random SL variates
SL_m  = f_johnson_fit(Y,'M',2);              % fit using moments
SL_q  = f_johnson_fit(Y,'Q',0);              % fit using quantiles
SL_sl = f_johnson_fit(Y,'M_SL',1);           % fit using moments, force SL
[coef' SL_m.coef' SL_q.coef' SL_sl.coef']    % compare actual vs. estimates
% ans =
%             1      -5.3457      -7.0955       1.0264
%             3       2.8709       2.9917       2.9795
%             0     0.068467     0.016351    0.0073852
%             1      0.20633      0.13195            1
% 
% Note: both 'M' & 'Q' methods fit an 'SU' family distribution; to fit an
% 'SL' family we had to specify 'M_SL' as the method. Also, we use a slightly
% different notation to specify the coefficients than those listed in Table 3 of
% George & Ramachandran (2011). They specify coefs = [1 3 0 10] in their third
% example, but the last parameter (i.e, LAMBDA) must equal 1 for 'SL' family
% distributions in our implementation.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            SN Family (Normal distribution via identity transform):           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;                  % clear the workspace
rng('default');                             % reset for repeatability
n    = 1000*100;                            % sample size
coef = [0 1 0 1];                           % SN distribution parameters
Y    = f_johnson_rnd(coef,'SN',n);          % generate random SN variates
SN_m = f_johnson_fit(Y,'M',1);              % fit using moments
SN_q = f_johnson_fit(Y,'Q',0);              % fit using quantiles
[coef' SN_m.coef' SN_q.coef']               % compare actual vs. estimates
% ans =
%             0   0.00079963      0.87006
%             1       1.0031       13.434
%             0            0      0.86896
%             1            1       13.367
% 
% Note 'M' fitted an 'SN' distribution after failing to converge on an 'SB'
% solution, while 'Q' produced an 'SU' distribution.
% 
% 
[mean(Y) std(Y) skewness(Y) kurtosis(Y)]
% ans =
% 
%   -0.00079713      0.99687    0.0042474       2.9891
% -> these should be [0 1 0 3] for a Normal distribution

% See if method 'M' recovers the correct 'identity' transform for a standard
% normal distirubtion:
f_johnson_M(0,1,0,3)
% ans = 
%     coef: [0 1 0 1]
%     type: 'SN'
