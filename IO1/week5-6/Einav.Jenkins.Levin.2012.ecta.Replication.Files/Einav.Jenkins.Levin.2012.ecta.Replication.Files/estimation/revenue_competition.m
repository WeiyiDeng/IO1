%--------------------------------------------------------------------------
% revenue_competition.m
%--------------------------------------------------------------------------

function [conrev conprof urevinc uprfinc urevent uprfent nonsale incsale entsale divsale] = revenue_competition(pKtemp,dKinc,dKent,par,shadowcost)

global N G lK X1 X2 X3 X4 nu xi ep et;
global SLTXRT docfeeK slsfeeK tlsfeeK cK;
global S APRUSE TERM X5 X6;
global mxi vxi mep vep met vet;

%--------------------------------------------------------------------------
% define paramters and variables for probability calculations
%--------------------------------------------------------------------------

% define parameters
[Lambda Delta Beta Gamma] = define_parameters(par);

% generate regression variables that depend on pK and dK
dKinc2 = dKinc.*dKinc;
dKent2 = dKent.*dKent;
ptilda = pKtemp.*(1+SLTXRT) + docfeeK + slsfeeK + tlsfeeK;

%calculate down payment as a function of minimum down
XBeta = [X2 pKtemp]*Beta;
downinc = max(XBeta + ep, dKinc);
downent = max(XBeta + ep, dKent);

%calculate minimum down payment indicator
mindinc = downinc <= dKinc;
mindent = downent <= dKent;

%calculate sale indicator
XDinc = [X4 pKtemp mindinc.*dKinc mindinc.*dKinc2]*Delta;
XDent = [X4 pKtemp mindent.*dKent mindent.*dKent2]*Delta;
saleinc = (XDinc + xi > 0) & (XDinc >= XDent);
saleent = (XDent + xi > 0) & (XDinc <= XDent);

%determine sale outcomes
divsale = saleinc & saleent;
nonsale = ~saleinc & ~saleent;
incsale = saleinc & ~saleent;
entsale = ~saleinc & saleent;

%determine down payments
downtemp=zeros(N,1);
downtemp(incsale) = downinc(incsale);
downtemp(entsale) = downent(entsale);
downtemp(divsale) = downinc(divsale);
downtemp(nonsale) = -max(-downinc(nonsale),-downent(nonsale));

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
urevinc = incsale.*conrev  + 0.5*divsale.*conrev;
uprfinc = incsale.*conprof + 0.5*divsale.*conprof;
urevent = entsale.*conrev  + 0.5*divsale.*conrev;
uprfent = entsale.*conprof + 0.5*divsale.*conprof;

%--------------------------------------------------------------------------
% end of function
%--------------------------------------------------------------------------
