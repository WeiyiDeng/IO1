function [Pj]=entryProb(log_Cr,Pr,pi,kappa,tau,N)

Pj = normcdf((log_Cr+log(pi)-log((1+(N-1)*Pr)*kappa))/tau);

end

