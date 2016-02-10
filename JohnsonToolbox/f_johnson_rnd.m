function R = f_johnson_rnd(coef,type,n,m)
% - generate random numbers from a Johnson distribution
% 
% USAGE: R = f_johnson_rnd(Y,coef,type);
%
% coef = parameters of the Johnson distribution as: [gamma delta xi lambda]
% type = type of Johnson distribution as: 'SL', 'SU', 'SB', or 'SN'
% n    = number of rows of random numbers required
% m    = number of columns of random numbers required              (default = 1)
% 
% R    = an N X M matrix of Johnson distributed random numbers
%
% SEE ALSO: f_johnson_z2y, randn

% -----Author:-----
% by David L. Jones, Mar-2014
%
% This file is part of the 'JOHNSON CURVE TOOLBOX FOR MATLAB'
% and is released under the BSD 2-clause license.

% -----Set defaults & check input:-----
if (nargin < 4), m = 1; end % default return 1 colum
% -------------------------------------

% Generate random Johnson variates:
R = f_johnson_z2y(randn(n,m),coef,type);
