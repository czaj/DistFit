function [AIC,AICc,BIC] = f_johnson_aic(X,coef,type)
% - calculate AIC, AICc, and BIC for a Johnson distribution
%
% USAGE: [AIC,AICc,BIC] = f_johnson_aic(X,coef,type);
%
% X    = column vector of Johnson variates
% coef = parameters of the Johnson distribution as: [gamma delta xi lambda]
% type = type of Johnson distribution as: 'SL', 'SU', 'SB', or 'SN'
%
% AIC  = Akaike's Information Criterion
% AICc = AIC corrected to avoid bias between # observations and # parameters
% BIC  = Schwarz's Bayesian Information Criterion
%
% SEE ALSO: f_johnson_lik

% -----Notes:-----
% This function returns Akaike's Information Criterion (AIC), a corrected form
% of AIC for small sample sizes (AICc), and Schwarz's Bayesian Information
% Criterion (BIC) for a Johnson distribution with parameters COEF and family
% TYPE, evaluated for the values in X.

% -----References:-----
% Anderson, D. R., K. P. Burnham, and W. L. Thompson. 2000. Null hypothesis
%  testing: problems, prevalence, and an alternative. J. Wildl. Manage.
%  64(4): 912-923.
% Burnham, K. P. and D. R. Anderson. 2001. Kullback-Leibler information as
%  a basis for strong inference in ecological studies. Wildlife Research 28:
%  111-119.

% -----Author:-----
% by David L. Jones, Apr-2014
%
% This file is part of the 'JOHNSON CURVE TOOLBOX FOR MATLAB'
% and is released under the BSD 2-clause license.

% -----Set defaults & check input:-----
% Check size of input:
if (size(X,2)>1)
   error('X must be a column vector!');
end

% Check for missing values:
if any(isnan(X))
   error('X contains NaN''s!');
end
% -------------------------------------

% Get # of parameters for the distribution:
switch type
   case 'SN' % xi & lambda specify mean and SD
      k = 2;
   case 'SL' % xi specifies the bound, lambda whether it's pos/neg skewed
      k = 3;
   case 'SU'
      k = 4;
   case 'SB' % xi & xi+lambda specify lower & upper bounds
      k = 4;
   case 'ST' % gamma is set to 0
      k = 3;
   otherwise
      error('Unknown TYPE!');
end

n    = size(X,1)    ;                   % # observations
LL   = -1 * f_johnson_lik(X,coef,type); % log-likelihood
AIC  = (-2*LL) + (2*k);                 % Anderson et al. (2000)
BIC  = (-2*LL) + (log(n)*k);
AICc = AIC + ((2*k*(k+1)) / (n-k-1));   % Anderson et al. (2000)
