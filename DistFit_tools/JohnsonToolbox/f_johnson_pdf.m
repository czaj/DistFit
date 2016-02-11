function Y = f_johnson_pdf(X,coef,type)
% - probability density function for a Johnson distribution
%
% USAGE: Y = f_johnson_pdf(X,coef,type);
%
% X    = column vector of Johnson variates
% coef = parameters of the Johnson distribution as: [gamma delta xi lambda]
% type = type of Johnson distribution as: 'SL', 'SU', 'SB', or 'SN'
%
% Y    = probability densities
%
% SEE ALSO: f_johnson_cdf, f_johnson_inv, normpdf

% -----Notes:-----
% This function returns the PDF of the Johnson distribution with parameters COEF
% and family TYPE, evaluated at the values in X.

% -----References:-----
% Based on a ported to Matlab from Robert E. Wheeler's 'dists.cc' C code in his
% 'SuppDists package for R'.
%
% Hill, I. D., R. Hill, and R. L. Holder, 1976. Algorithm AS 99: Fitting Johnson
%  curves by moments. Journal of the Royal Statistical Society. Series C
%  (Applied Statistics) 25(2): 180-189.

% -----Author:-----
% by David L. Jones, Apr-2014
%
% This file is part of the 'JOHNSON CURVE TOOLBOX FOR MATLAB'
% and is released under the BSD 2-clause license.

% Apr-2014: now makes sure values of X don't exceed boundaries of bounded
%           distributions

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

u      = (X-xi)/lambda;
ratio  = delta/lambda;

switch type
   case 'SN'
      fu = u;
      D  = ratio;
   case 'SL'
      D  = ratio./u;
      fu = log(u);
   case 'SU'
      fu = u+sqrt(1+u.*u);
      D  = ratio./sqrt(1.0+u.*u);
      fu = log(fu);
   case 'SB'
      fu = u./(1-u);
      D  = ratio./(u.*(1.0-u));
      fu = log(fu);
   otherwise
      error('Unknown TYPE!');
end

Z = gamma+delta.*fu;

% Return NaN for values exceeding bounds of fitted distribution:
switch type
   case 'SB'
      Z(X<=xi | X>=(xi+lambda)) = NaN;
   case 'SL'
      switch lambda
         case  1 % bounded but positively skewed
            Z(X<xi) = NaN;
         case -1 % bounded but negatively skewed
            Z(X>xi) = NaN;
         otherwise
            error('\nLAMBDA must be +1 or -1!')
      end
end

% Get probability densities:
Y = normpdf(Z,0,1).*D;
