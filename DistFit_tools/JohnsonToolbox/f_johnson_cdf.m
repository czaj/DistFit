function Y = f_johnson_cdf(X,coef,type)

save tmpFJ

% edited to
% - accept arrays as inputs
% - allow for -Inf, +Inf values
% - allow for 0 in SB

% - cumulative probability density function for a Johnson distribution
%
% USAGE: Y = f_johnson_cdf(X,coef,type);
%
% X    = column vector of Johnson variates
% coef = parameters of the Johnson distribution as: [gamma delta xi lambda]
% type = type of Johnson distribution as: 'SL', 'SU', 'SB', or 'SN'
% 
% Y    = cumulative probability densities
%
% SEE ALSO: f_johnson_inv, f_johnson_pdf, normcdf

% -----Notes:-----
% This function returns the CDF of the Johnson distribution with parameters COEF
% and family TYPE, evaluated at the values in X.

% -----References:-----
% Ported to Matlab from Robert E. Wheeler's 'dists.cc' C code in his 'SuppDists
% package for R'.

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

% Extract coefficients:
gamma  = coef(1);
delta  = coef(2);
xi     = coef(3);
lambda = coef(4);

u = (X-xi)/lambda;

switch type
   case 'SN'
      % do nothing
   case 'SL'
      u = log(u);
   case 'SU'
      u = u + sqrt(1+u.*u);
      u = log(u);
   case 'SB'
      if (any(u <= 0) || any(u>=1.0))
         error('SB values out of range!');
      else
         u = u ./ (1-u);
         u = log(u);
      end
   otherwise
      error('Unknown TYPE!');
end

Z = gamma+delta*u;
Y = normcdf(Z,0,1);
