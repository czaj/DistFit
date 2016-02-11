function Y = f_johnson_z2y(Z,coef,type)
% - transform standard normal variates to Johnson variates
% 
% USAGE: Y = f_johnson_z2y(Z,coef,type)
% 
% Z    = column vector of standard normal variates
% coef = parameters of the Johnson distribution as: [gamma delta xi lambda]
% type = type of Johnson distribution as: 'SL', 'SU', 'SB', or 'SN'
% 
% Y = column vector of Johnson variates
% 
% SEE ALSO: f_johnson_fit, f_johnson_y2z

% -----Notes:-----
% This function is used to convert variates from a standard normal distribution
% to those from the corresponding Johnson distribution. It was ported to Matlab
% from Hill et al.'s (1976) original AS-100 FORTRAN source code, which was
% obtained from: http://lib.stat.cmu.edu/apstat/100.
% 
% One of the main purposes of this function is to generate a random sample drawn
% from a specific Johnson distribution. If Z is a random sample taken from the
% standard normal distribution, e.g., Z = randn(n,1), then Y will be a random
% sample taken from a Johnson distribution defined by the gamma, delta, xi, and
% lambda parameters specified in the 'coef' variable.

% -----References:-----
% George, F. and K. M. Ramachandran. 2011. Estimation of parameters of Johnson's
%  system of distributions. Journal of Modern Applied Statistical Methods 10(2):
%  494-504.
% Hill, I. D. 1976. Algorithm AS 100: Normal-Johnson and Johnson-Normal
%  Transformations. Journal of the Royal Statistical Society. Series C (Applied
%  Statistics) 25(2): 190-192.
% Hill, I. D., R. Hill, and R. L. Holder, 1976. Algorithm AS 99: Fitting Johnson
%  Curves by Moments. Journal of the Royal Statistical Society. Series C
%  (Applied Statistics) Vol. 25, No. 2, 180-189.

% -----Author:-----
% by David L. Jones, Mar-2014
%
% This file is part of the 'JOHNSON CURVE TOOLBOX FOR MATLAB'
% and is released under the BSD 2-clause license.

% -----Set defaults & check input:-----
% Check size of input:
if (size(Z,2)>1)
   error('Z must be a column vector!');
end

% Check for missing values:
if any(isnan(Z))
   error('Z contains NaN''s!');
end

% Ensure input is uppercase:
type = upper(type);

% Check coefficients (George & Ramachandran, 2011):
if any([coef(2) coef(4)]<=0)
   error('DELTA & LAMBDA must be > 0!')
end
% ----------------------

nr = size(Z,1); % # rows of input

% Rename input variables to follow original algorithm:
SNV = Z;

% Extract coefficients:
GAMMA = coef(1); % gamma;
DELTA = coef(2); % delta;
XI    = coef(3); % xi;
XLAM  = coef(4); % lambda;

% Define constants:
HALF = repmat(0.5,nr,1);
ONE  = ones(nr,1);

switch type
   case 'SL' % Type SL: Exponential transformation (lognormal distribution):
      AJV = XLAM * exp((XLAM * SNV - GAMMA) / DELTA) + XI;

   case 'SU' % Type SU: Hyperbolic sine transformation (unbounded)
      W   = exp((SNV - GAMMA) / DELTA);
      W   = HALF .* (W - ONE ./ W);
      AJV = XLAM * W + XI;

   case 'SB' % Type SB: Logistic transformation (bounded)
      W   = (SNV - GAMMA) / DELTA;
      V   = exp(-abs(W));
      V   = (ONE - V) ./ (ONE + V);
      AJV = (HALF * XLAM) .* (sub_sign(V, W) + ONE) + XI;
            
   case 'SN' % Type SN: Identity transformation (normal distribution)
      AJV = (SNV - GAMMA) / DELTA;

   otherwise
      error('Unknown distribution specified in TYPE!');
end

% Rename output variables:
Y = AJV;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SUBFUNCTION:                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = sub_sign(A,B)
% - port of FORTRAN 'SIGN' function
% 
% If B\ge 0 then the result is ABS(A), else it is -ABS(A).
A      = abs(A);
A(B<0) = A(B<0) * -1;
