function [LL]=LogLikBayesNash(par0,log_Cr,Nr)

N = max(Nr);
ll = 0;
for i = 1:length(log_Cr)
    Pr = bayesNash(log_Cr(i),par0(1),par0(2),par0(3),N);
    lr = factorial(N)/(factorial(Nr(i))*factorial(N-Nr(i)))*Pr^Nr(i)*(1-Pr)^(N-Nr(i));       % binomial distr
    ll = ll+log(lr);
end

LL = -ll;

end