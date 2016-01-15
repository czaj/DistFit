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
elseif dist==9
    f = log(cdf('BirnbaumSaunders',bounds(:,2),B(1),B(2)) - cdf('BirnbaumSaunders',bounds(:,1),B(1),B(2)));
elseif dist==10
    f = log(cdf('ExtremeValue',bounds(:,2),B(1),B(2)) - cdf('ExtremeValue',bounds(:,1),B(1),B(2)));
elseif dist==11
    f = log(cdf('GeneralizedExtremeValue',bounds(:,2),B(1),B(2),B(3)) - cdf('GeneralizedExtremeValue',bounds(:,1),B(1),B(2),B(3)));
elseif dist==12
    f = log(cdf('GeneralizedPareto',bounds(:,2),B(1),B(2)) - cdf('GeneralizedPareto',bounds(:,1),B(1),B(2)));
elseif dist==13
    f = log(cdf('InverseGaussian',bounds(:,2),B(1),B(2)) - cdf('InverseGaussian',bounds(:,1),B(1),B(2)));
elseif dist==14
    f = log(cdf('Nakagami',bounds(:,2),B(1),B(2)) - cdf('Nakagami',bounds(:,1),B(1),B(2)));
elseif dist==15
    f = log(cdf('Rician',bounds(:,2),B(1),B(2)) - cdf('Rician',bounds(:,1),B(1),B(2)));
elseif dist==16
    f = log(cdf('tLocationScale',bounds(:,2),B(1),B(2),B(3)) - cdf('tLocationScale',bounds(:,1),B(1),B(2),B(3)));
elseif dist==17 % Johnson SU
    f = log(f_johnson_cdf(bounds(:,2),B, 'SU') - f_johnson_cdf(bounds(:,1),B, 'SU'));
elseif dist==18 % Johnson SB - nie dzia³a, gdy any(u <= 0) || any(u>=1.0)
    f = log(f_johnson_cdf(bounds(:,2),B, 'SB') - f_johnson_cdf(bounds(:,1),B, 'SB'));
elseif dist==19 % Burr
    f = log(cdf('Burr',bounds(:,2),B(1),B(2),B(3)) - cdf('Burr',bounds(:,1),B(1),B(2),B(3)));
elseif dist==20 % F
    f = log(cdf('F',bounds(:,2),B(1),B(2)) - cdf('F',bounds(:,1),B(1),B(2)));
elseif dist==21 % negative binomial
    f = log(cdf('Negative Binomial',bounds(:,2),B(1),B(2)) - cdf('Negative Binomial',bounds(:,1),B(1),B(2)));
elseif dist==22 % uniform
    f = log(cdf('Uniform',bounds(:,2),B(1),B(2)) - cdf('Uniform',bounds(:,1),B(1),B(2)));
elseif dist==23 % generalized inverse Gaussian
    f = 
elseif dist==24 % sinh-arcsinh
    f = 
end

