function f = LL_DistFit(bounds, dist, B)

% save CDF_WTP_tmp
% return

if dist==0
    f = log(normcdf(bounds(:,2),B(1),B(2)) - normcdf(bounds(:,1),B(1),B(2)));
elseif dist==1
    f = log(logncdf(bounds(:,2),B(1),B(2)) - logncdf(bounds(:,1),B(1),B(2)));
elseif dist==2
    f = log(cdf('Exponential',bounds(:,2),B) - cdf('Exponential',bounds(:,1),B));
elseif dist==3
    f = log(cdf('Logistic',bounds(:,2),B(1),B(2)) - cdf('Logistic',bounds(:,1),B(1),B(2)));
elseif dist==4
    f = log(cdf('LogLogistic',bounds(:,2),B(1),B(2)) - cdf('LogLogistic',bounds(:,1),B(1),B(2)));
elseif dist==5 
%     f = log(cdf('Poisson',bounds(:,2),B) - cdf('Poisson',bounds(:,1),B));
    f = log(max(cdf('Poisson',bounds(:,2),B) - cdf('Poisson',bounds(:,1),B),realmin));
elseif dist==6
    f = log(cdf('Weibull',bounds(:,2),B(1),B(2)) - cdf('Weibull',bounds(:,1),B(1),B(2)));
elseif dist==7
    f = log(cdf('Gamma',bounds(:,2),B(1),B(2)) - cdf('Gamma',bounds(:,1),B(1),B(2)));
elseif dist==8
%     f = log(cdf('Rayleigh',bounds(:,2),B) - cdf('Rayleigh',bounds(:,1),B));
    f = log(max(cdf('Rayleigh',bounds(:,2),B) - cdf('Rayleigh',bounds(:,1),B),realmin));
end


