function X = f_johnson_inv(P,coef,type)
% - inverse of the CDF for a Johnson distribution
%
% USAGE: X = f_johnson_inv(P,coef,type);
%
% P    = column vector of probabilities
% coef = parameters of the Johnson distribution as: [gamma delta xi lambda]
% type = type of Johnson distribution as: 'SL', 'SU', 'SB', or 'SN'
% 
% X    = Johnson variates
%
% SEE ALSO: f_johnson_cdf, norminv

% -----Notes:-----
% This function returns the inverse of the CDF of the Johnson distribution with
% parameters COEF and family TYPE, evaluated at the values in P.

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
if (size(P,2)>1)
   error('X must be a column vector!');
end

% Check for missing values:
if any(isnan(P))
   error('P contains NaN''s!');
end
% -------------------------------------

% Extract coefficients:
gamma  = coef(1);
delta  = coef(2);
xi     = coef(3);
lambda = coef(4);

Z = norminv(P,0,1);
u = (Z-gamma)/delta;

switch type
   case 'SN'
      % do nothing
   case 'SL'
      u = exp(u);
   case 'SU'
      u = exp(u);
      u = ((u.*u)-1)./(2*u);
   case 'SB'
      u = exp(u);
      u = u./(1+u);
   otherwise
      error('Unknown TYPE!');
end

X = xi+lambda*u;
