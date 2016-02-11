function nLL = f_johnson_lik(X,coef,type)
% - negative log-likelihood for a Johnson distribution
%
% USAGE: nLL = f_johnson_lik(X,coef,type);
%
% X    = column vector of Johnson variates
% coef = parameters of the Johnson distribution as: [gamma delta xi lambda]
% type = type of Johnson distribution as: 'SL', 'SU', 'SB', or 'SN'
% 
% nLL = negative log-likelihood
%
% SEE ALSO: f_johnson_aic, f_johnson_pdf, normlike

% -----Notes:-----
% This function returns the negative log-likelihood of the Johnson distribution
% with parameters COEF and family TYPE, evaluated at the values in X.

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

% Compute log-likelihoods:
LL = log(f_johnson_pdf(X,coef,type));

% Negative sum of log-likelihoods (excluding NaN's):
nLL = -sum(LL(~isnan(LL)));
