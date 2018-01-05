function [p]=bayesNash(log_Cr,delta,kappa,tau,N)

pi = @(n) delta^(n-1);

entryProbFun = @(Pr) normcdf((log_Cr+log(pi(1+(N-1)*Pr))-log((1+(N-1)*Pr)*kappa))/tau)-Pr;

% % bisection method
% p_L = 0;
% p_U = 1;
% p_mid = 0.5;
% % while abs(entryProbFun(p_mid)) > eps
% % while abs(entryProbFun(p_mid)) > 1.0000e-10
% while abs(entryProbFun(p_mid)) > 1.0000e-5
%     if entryProbFun(p_L)*entryProbFun(p_mid)<0
%         p_U = p_mid;
%     elseif entryProbFun(p_U)*entryProbFun(p_mid)<0
%         p_L = p_mid;
%     end
%     p_mid = (p_L+p_U)/2;
% end
% p = p_mid;

% use zeroin algorithm by matlab (much faster !!)
p = fzero(entryProbFun,[0 1]);

end