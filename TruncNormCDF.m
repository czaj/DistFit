function p = TruncNormCDF(x, mu, sigma, a, b)

% save TruncNormCDF_tmp
% return

%         mu - mean
%         sigma - st.dev. 
%         a - lower bound
%         b - upper bound


p = (cdf('Normal',(x-mu)./sigma,mu,sigma) - cdf('Normal',(a-mu)./sigma,mu,sigma)) ./ ...
    (cdf('Normal',(b-mu)./sigma,mu,sigma) - cdf('Normal',(a-mu)./sigma,mu,sigma));

