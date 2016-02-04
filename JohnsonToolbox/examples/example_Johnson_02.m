% Examples of using the Johnson System of Distributions (Johnson Curves)
% 
% by David L. Jones, Apr-2014
%
% This file is part of the 'JOHNSON CURVE TOOLBOX FOR MATLAB'
% and is released under the BSD 2-clause license.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Weekly S&P Index return rates:                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example considers weekly opening and closing prices of the S&P stock
% market index for the period Apr-1995 to Mar-2005 and follows that presented in
% http://www.financialwisdomforum.org/gummy-stuff/johnson.htm.

% Clear the workspace:
clc; clear all; close all;

% Load the data (weekly opening and weekly closing price are listed in columns
% 1 & 2, respectively):
X = load('SPindex.dat');

% Plot a times series of the closing price:
figure; set(gcf,'color','w');
plot(X(:,2),'b-')
title({'S&P Index'});
xlabel('Week')
ylabel('Closing Price');

% Calculate percent weekly returns:
Y = (X(:,2)./X(:,1)-1)*100;

% Visualize weekly return data:
figure; set(gcf,'color','w');
hist(Y,100);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[1 1 1]*0.65,'EdgeColor','none')
title('S&P Index');
xlabel('% return')
ylabel('Frequency (weeks)');

% Fit a Johnson Curve:
jsn_m = f_johnson_fit(Y,'M',2); % fit using moments
jsn_q = f_johnson_fit(Y,'Q',0); % fit using quantiles
[jsn_m.coef' jsn_q.coef']       % compare estimates
% 
% ans =
%        0.2733      -1.7913
%        1.8624       4.2128
%       0.83092      -3.6022
%        3.7515       8.4897

% What the probability your return will be less than 2%?
prb = f_johnson_cdf(1.5,jsn_m.coef,jsn_m.type);
h   = text(1.5,prb,{sprintf('%2.0f%%',prb*100)});
set(h,'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center',...
   'VerticalAlignment','middle', 'Color',[0 0 0],'BackgroundColor',...
   [1 1 1],'Margin',0.5);

% What is the probability your return will be 4%?
pval = 1 - f_johnson_cdf(4,jsn_m.coef,jsn_m.type)
% pval = 0.04436

% Compare the fit of the Johnson 'SU' distribution with that of, say, a two
% parameter normal distribution defined by mu = 1038.2 and sigma = 257.68:
% 
% Show measures of fit for Normal distribution (calculated elswehere):
% normal.nLL  = 7262.9
% normal.AIC  = 14530
% normal.AICc = 14530
% normal.BIC  = 14540

% Calculate measures of fit for Johnson 'SU' distribution:
nLL            = f_johnson_lik(Y,jsn_m.coef,jsn_m.type);
[AIC,AICc,BIC] = f_johnson_aic(Y,jsn_m.coef,jsn_m.type);
[nLL AIC AICc BIC]'
% ans =  1175.2
%        2358.5
%        2358.5
%        2375.5

% -> Smaller is better so all of these goodness-of-fit measures indicate the
%    Johnson 'SU' distribution provides a better fit to the data than the
%    Normal distribution.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Example presented in Hill et al. (1976):                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This example demonstrates the accuracy of using Hill et al.'s (1976) method
% based on moments to (1) find the parameters of a Johnson curve; (2) transform
% Johnson variates to normal variates; and (3) use the transformed variates to
% assess statistical significance.
% 
% Hill, I. D., R. Hill, and R. L. Holder, 1976. Algorithm AS 99: Fitting Johnson
%  curves by moments. Journal of the Royal Statistical Society. Series C
%  (Applied Statistics) 25(2): 180-189.

% Clear the workspace:
clc; clear all; close all;

% Specify range of degrees-of-freedom:
df = (1:4)';

% Find parametes for Johnson curves using moments:
for i=1:numel(df);
   jsn{i} = f_johnson_M(df(i),sqrt(df(i)+df(i)),sqrt(8/df(i)),12/df(i)+3);
end

% Define target percentage points
pp = [0.50 0.10 0.01];

% Get Chi-square values for df's at target percentage points:
C = nan(size(df,1),size(pp,2));
for i=1:numel(df);
   C(i,:) = chi2inv(1-pp,df(i));
end
C
% C =
%       0.45494       2.7055       6.6349
%        1.3863       4.6052       9.2103
%         2.366       6.2514       11.345
%        3.3567       7.7794       13.277

% Confirm values match target percentage points for a Chi-square distribution:
C_pp = nan(size(C));
for i=1:size(df,1);
   C_pp(i,:) = 1-chi2cdf(C(i,:),df(i));
end
C_pp
% C_pp =
%           0.5          0.1         0.01
%           0.5          0.1         0.01
%           0.5          0.1         0.01
%           0.5          0.1         0.01

% Transform Chi-square variates to normal variates using Johnson transform:
Z = nan(size(C));
for i=1:size(C,1);
   Z(i,:) = f_johnson_y2z(C(i,:)',jsn{i}.coef,jsn{i}.type)';
end
Z
% Z =
%      -0.09689       1.3092       2.3094
%     -0.030035       1.2975       2.3069
%     -0.011945        1.291       2.3098
%    -0.0052176       1.2874       2.3127

% Show percentage points for a normal distribution match the target values:
1-normcdf(Z)
% ans =
%       0.53859     0.095235     0.010461
%       0.51198     0.097229     0.010529
%       0.50477      0.09836     0.010449
%       0.50208     0.098979     0.010371

% Alternative method:
P = nan(size(C));
for i=1:size(P,1);
   P(i,:) = 1- f_johnson_cdf(C(i,:)',jsn{i}.coef,jsn{i}.type)';
end
P
% P =
%       0.53859     0.095235     0.010461
%       0.51198     0.097229     0.010529
%       0.50477      0.09836     0.010449
%       0.50208     0.098979     0.010371
