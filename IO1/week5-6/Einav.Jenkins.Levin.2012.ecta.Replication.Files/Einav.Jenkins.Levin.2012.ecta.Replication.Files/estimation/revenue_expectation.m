%--------------------------------------------------------------------------
% revenue_expectation.m
%--------------------------------------------------------------------------

function [saleind downtemp mindind fractemp deftind pvpmts pvrecovery conrev concost conprof uncrev uncprof pr_less pr_sale pr_nosl pr_deft fracsum] = revenue_expectation(pKtemp,dKtemp,par,shadowcost)

global G X2 X3 X4 xi ep et;
global SLTXRT docfeeK slsfeeK tlsfeeK cK;
global S APRUSE TERM X5 X6;
global mxi vxi mep vep met vet;

%--------------------------------------------------------------------------
% define paramters and variables for probability calculations
%--------------------------------------------------------------------------

% define parameters
[Lambda Delta Beta Gamma] = define_parameters(par);

% generate regression variables that depend on pK and dK
dKtemp2 = dKtemp.*dKtemp;
ptilda = pKtemp.*(1+SLTXRT) + docfeeK + slsfeeK + tlsfeeK;

%calculate down payment as a function of minimum down
XBeta = [X2 pKtemp]*Beta;
downtemp = max(XBeta + ep, dKtemp);

%calculate minimum down payment indicator
mindind = downtemp <= dKtemp;

%calculate sale indicator
XDelta = [X4 pKtemp dKtemp dKtemp2]*Delta;
saleind = (XDelta + xi > 0);

%--------------------------------------------------------------------------
% calculate simulated payment and profit outcomes
%--------------------------------------------------------------------------
%calculate fraction of payments made
XGamma = [X3 ptilda-downtemp]*Gamma;
fractemp = exp(XGamma + et);
fractemp(fractemp<1/(G*4)) = 0;
fractemp(fractemp>1) = 1;

%calculate default indicator
deftind = fractemp<1;

%calculate revenues conditional on fracpaid
[pvlpratio pvrecovery] = revenue_coefficients(fractemp,S,APRUSE,TERM,X5,X6);
pvpmts = pvlpratio.*(ptilda-downtemp);
conrev = downtemp + pvpmts + pvrecovery;

%calculate profits conditional on sale
concost = cK + pKtemp.*SLTXRT + tlsfeeK;
conprof = conrev - concost - shadowcost;

%calculate unconditional profits
uncrev  = saleind.*conrev;
uncprof = saleind.*conprof;

%calculate probability of sale without simulation (not currently used)
pr_less = 1e-99 + normcdf(dKtemp-XBeta-mep,0,vep^0.5);
pr_sale = 1e-99 + normcdf(XDelta+mxi,0,vxi^0.5);
pr_nosl = 1e-99 + (1-normcdf([X4 pKtemp dKtemp dKtemp2]*Delta+mxi,0,vxi^0.5));
pr_deft = 1e-99 + normcdf(log(1)-XGamma-met,0,vet^0.5);

%calculate expected profit from default
fracsum = 0;
for g=1:G-1,
    frachi = (g+1)/G;
    fraclo = g/G;
    prfrac = normcdf(log(frachi)-XGamma-met,0,vet^0.5)-normcdf(log(fraclo)-XGamma-met,0,vet^0.5);
    [pvlpratio pvrecovery] = revenue_coefficients((frachi+fraclo)/2,S,APRUSE,TERM,X5,X6);
    pvpmts = pvlpratio.*(ptilda-downtemp);
    conrev = downtemp + pvpmts + pvrecovery;
    conprof = conrev - concost - shadowcost;
    fracsum = fracsum + prfrac.*conprof;
end

%calculate expected profit from full payment
prpaid = 1-normcdf(log(1)-XGamma-met,0,vet^0.5);
[pvlpratio pvrecovery] = revenue_coefficients(1,S,APRUSE,TERM,X5,X6);
pvpmts = pvlpratio.*(ptilda-downtemp);
conrev = downtemp + pvpmts + pvrecovery;
conprof = conrev - concost - shadowcost;
fracsum = fracsum + prpaid.*conprof;

%calculate expected profit from no payments
prnone = normcdf(log(1/G)-XGamma-met,0,vet^0.5);
[pvlpratio pvrecovery] = revenue_coefficients(0,S,APRUSE,TERM,X5,X6);
pvpmts = pvlpratio.*(ptilda-downtemp);
conrev = downtemp + pvpmts + pvrecovery;
conprof = conrev - concost - shadowcost;
fracsum = fracsum + prnone.*conprof;

%--------------------------------------------------------------------------
% end of function
%--------------------------------------------------------------------------
