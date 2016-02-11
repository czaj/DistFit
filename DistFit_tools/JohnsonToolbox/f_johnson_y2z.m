function [Z,ifault] = f_johnson_y2z(Y,coef,type)
% - transform Johnson variates to standard normal variates
%
% USAGE: [Z,ifault] = f_johnson_y2z(Y,coef,type);
%
% Y    = column vector of Johnson variates
% coef = parameters of the Johnson distribution as: [gamma delta xi lambda]
% type = type of Johnson distribution as: 'SL', 'SU', 'SB', or 'SN'
%
% Z      = column vector of standard normal variates
% ifault = index to rows not transformed
%
% SEE ALSO: f_johnson_fit, f_johnson_z2y

% -----Notes:-----
% This function is used to convert variates from a Johnson distribution to those
% from the corresponding standard normal distribution. It was ported to Matlab
% from Hill et al.'s (1976) original AS-100 FORTRAN source code, which was
% obtained from: http://lib.stat.cmu.edu/apstat/100.
%
% One of the main purposes of this function is to convert observations of a
% test statistic that follows (as a null hypothesis) a distribution with known
% mean, standard deviation, skewness, and kurtosis to observations that follow a
% standard normal distribution in order to assess statistical significance.

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

% Apr-2014: now properly handles cases where ifault=2

% -----Set defaults & check input:-----
% Check size of input:
if (size(Y,2)>1)
   error('Y must be a column vector!');
end

% Check for missing values:
if any(isnan(Y))
   error('Y contains NaN''s!');
end

% Ensure input is uppercase:
type = upper(type);

% Check coefficients (George & Ramachandran, 2011):
if any([coef(2) coef(4)]<=0)
   error('DELTA & LAMBDA must be > 0!')
end
% ----------------------

nr = size(Y,1); % # rows of input

% Rename input variables to follow original algorithm:
AJV = Y;

% Extract coefficients:
GAMMA = coef(1); % gamma
DELTA = coef(2); % delta
XI    = coef(3); % xi
XLAM  = coef(4); % lambda

% Define constants:
ZERO = 0.0;
HALF = repmat(0.5,nr,1);
ONE  = 1.0;
C    = -63.0;

% Initialize
SNV    = zeros(nr,1);
ifault = zeros(nr,1);

switch type
   case 'SL' % Type SL: Exponential transformation (lognormal distribution):
      W = XLAM * (AJV - XI);
      
      % SNV = 0 if (W<=0):
      idx         = ~(W>ZERO);
      ifault(idx) = 2;
      
      % (W>0)
      idx      = (W>ZERO);
      SNV(idx) = XLAM * (log(W(idx)) * DELTA + GAMMA);
      
   case 'SU' % Type SU: Hyperbolic sine transformation (unbounded)
      W = (AJV - XI) / XLAM;
      
      % Greater than C:
      idx    = (W>C);
      W(idx) = sqrt(W(idx) .* W(idx) + ONE) + W(idx);
      
      % Less than or equal to C:
      idx    = (W<=C);
      W(idx) = -HALF(idx) ./ W(idx);
      
      SNV = log(W) * DELTA + GAMMA;
      
   case 'SB' % Type SB: Logistic transformation (bounded)
      W = AJV - XI;
      V = XLAM - W;
     
      % SNV = 0 if (W <= ZERO or V <= ZERO):
      idx         = ~(W > ZERO) & (V > ZERO);
      ifault(idx) = 2;
      
      % Both W & V > 0:
      idx      = (W > ZERO) & (V > ZERO);
      SNV(idx) = log(W(idx) ./ V(idx)) * DELTA + GAMMA;
   
   case 'SN' % Type SN: Identity transformation (normal distribution)
      SNV = DELTA * AJV + GAMMA;
      
   otherwise
      error('Unknown distribution specified in TYPE!');
end

% Rename output variables:
Z = SNV;
