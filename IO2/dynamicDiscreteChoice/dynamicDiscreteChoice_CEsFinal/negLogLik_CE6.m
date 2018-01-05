%% negLogLik.m
% The function computes minus the log partial likelihood and, optionally, the corresponding minus score and information matrix estimate for conditional choice data on the basic firm entry and exit model used as an example in CentER's Empirical Industrial Organization II.

%{
The function |negLogLik| computes minus the log partial likelihood for the conditional choice part of the data. Optionally, it also returns minus the corresponding score vector and an estimate of the information matrix for the parameter (sub)vector $\theta\equiv(\beta_0,\beta_1,\delta_1)'$ (the scores are specific to the estimation example in Section \ref{script}'s script and should be adapted for inference on other parameters).
%}
function nll = ...
         negLogLik_CE6(choices,iX,supportX,capPi,beta,delta,theta1,rho,...
         flowpayoffs,bellman,fixedPoint,tolFixedPoint)

nSuppX = size(supportX,1);
%{
Next, it computes the flow payoffs $u_0$ (|u0|) and $u_1$ (|u1|), the choice-specific net expected discounted values $U_0$ (|capU0|) and $U_1$ (|capU1|), their contrast $\Delta U$ (|deltaU|), and the implied probabilities $1/\left[1+\exp(\Delta U)\right]$ of not serving the market (|pExit|) for the inputted parameter values. Note that this implements the NFXP procedure's inner loop.
%}
delta_type1 = delta([1 2]);
delta_type2 = delta([1 3]);

[u0,u1] = flowpayoffs(supportX,beta,delta_type1);
[capU0,capU1] = fixedPoint(u0,u1,capPi,rho,tolFixedPoint,bellman,[],[]);
deltaU_type1 = capU1-capU0;
pExit_type1 = 1./(1+exp(deltaU_type1));

% clearvars u0 u1 capU0 capU1 
[u0,u1] = flowpayoffs(supportX,beta,delta_type2);
[capU0,capU1] = fixedPoint(u0,u1,capPi,rho,tolFixedPoint,bellman,[],[]);
deltaU_type2 = capU1-capU0;
pExit_type2 = 1./(1+exp(deltaU_type2));

%{
\paragraph{Log Partial Likelihood}
%}
laggedChoices = [zeros(1,size(choices,2));choices(1:end-1,:)];
p_type1 = choices + (1-2*choices).*pExit_type1(iX+nSuppX*laggedChoices);         % pExit(iX+nSuppX*laggedChoices) is T*N matrix   100*1000
p_type2 = choices + (1-2*choices).*pExit_type2(iX+nSuppX*laggedChoices);
% w: same as following:
% p = choices.*(1-pExit(iX+nSuppX*laggedChoices)) + (1-choices).*pExit(iX+nSuppX*laggedChoices);
% nll = -sum(sum(log(p)));
nll = sum(log(theta1.*exp(sum(log(p_type1),1))+(1-theta1).*exp(sum(log(p_type2),1))));
% nll = sum(theta1.*sum(log(p_type1),1)+(1-theta1).*sum(log(p_type2)));          % not the correct one

end
%{
When evaluated at an estimate of $\theta$, this provides an estimate of the expected partial likelihood information matrix for $\theta$. This estimate can in turn be used to estimate the variance-covariance matrix, and thus the standard errors, of the maximum partial likelihood estimator $\hat\theta$ of $\theta$. 
%}