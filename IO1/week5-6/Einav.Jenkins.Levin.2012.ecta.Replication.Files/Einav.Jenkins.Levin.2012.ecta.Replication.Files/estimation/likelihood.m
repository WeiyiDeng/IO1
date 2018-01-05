%--------------------------------------------------------------------------
% likelihood.m
%--------------------------------------------------------------------------

function [loglike] = likelihood(par)

global N K1 K2 K3 K4 X1 X2 X3 X4 lK dK dK2 afK pK tK;
global I less more none cenind obsdef fracpaid cenpoint G;

%--------------------------------------------------------------------------
% define parameters and fitted values
%--------------------------------------------------------------------------
%define parameters
Lambda   = par(1:K1+1);
Delta    = par(K1+2:K1+K4+4);
Beta     = par(K1+K4+5:K1+K4+K2+5);
Gamma    = par(K1+K4+K2+6:K1+K4+K2+K3+6);
Signu2   = par(K1+K2+K3+K4+7);
Sigxi2   = 1.0000;
Sigep2   = par(K1+K2+K3+K4+8);
Siget2   = par(K1+K2+K3+K4+9);
Rho_nuxi = par(K1+K2+K3+K4+10);
Rho_nuep = par(K1+K2+K3+K4+11);
Rho_nuet = par(K1+K2+K3+K4+12);
Rho_xiep = par(K1+K2+K3+K4+13);
Rho_xiet = par(K1+K2+K3+K4+14);
Rho_epet = par(K1+K2+K3+K4+15);

%ensure variances greater than 0
Signu2 = max(1e-9,Signu2);
Sigxi2 = max(1e-9,Sigxi2);
Sigep2 = max(1e-9,Sigep2);
Siget2 = max(1e-9,Siget2);

%ensure correlations between 0 and 1
Rho_nuxi = min(0.999,max(Rho_nuxi,-0.999));
Rho_nuep = min(0.999,max(Rho_nuep,-0.999));
Rho_nuet = min(0.999,max(Rho_nuet,-0.999));
Rho_xiep = min(0.999,max(Rho_xiep,-0.999));
Rho_xiet = min(0.999,max(Rho_xiet,-0.999));
Rho_epet = min(0.999,max(Rho_epet,-0.999));

%--------------------------------------------------------------------------
% simulate down unobservable (ep), conditional on price unobservable (nu)
%--------------------------------------------------------------------------
XLambda = [X1 lK]*Lambda; nu = pK - XLambda;

%compute conditional moments of epsilon distribution
A = Rho_nuep*(Signu2*Sigep2)^0.5;
B = Signu2;
C = nu';
mep = (A'*(B\C))';
vep = max(1e-9,Sigep2 - A'*(B\A));

%no assumption about non-sale down payments
XBeta  = [X2 pK]*Beta;
trunc1 = -999;
trunc2 = less.*dK-XBeta + ~less.*999;

%draw unobservables
stream = RandStream('mt19937ar','Seed',2);
RandStream.setDefaultStream(stream);
a = normcdf(trunc1,mep,vep^0.5);
b = normcdf(trunc2,mep,vep^0.5);
temp = a + (b-a).*rand(N,1);
ep = more.*(tK-XBeta) + ~more.*norminv(temp,mep,vep^0.5);
simless = XBeta + ep <= dK;

%--------------------------------------------------------------------------
% simulate sale unobservables (xi), conditional on nu and ep
%--------------------------------------------------------------------------
%compute conditional moments of xi distribution
A = [Rho_nuxi*(Signu2*Sigxi2)^0.5; Rho_xiep*(Sigxi2*Sigep2)^0.5];
B = [Signu2 Rho_nuep*(Signu2*Sigep2)^0.5; Rho_nuep*(Signu2*Sigep2)^0.5 Sigep2];
C = [nu'; ep'];
mxi = (A'*(B\C))';
vxi = max(1e-9,Sigxi2 - A'*(B\A));

%determine truncation points for draws (only used if rho_xiet ~= 0)
XDelta = [X4 pK dK dK2]*Delta;
trunc1 = I.*-XDelta + (1-I).*-999;
trunc2 = I.*+999 + (1-I).*-XDelta;

%draw unobservables
stream = RandStream('mt19937ar','Seed',1);
RandStream.setDefaultStream(stream);
a = normcdf(trunc1,mxi,vxi^0.5);
b = normcdf(trunc2,mxi,vxi^0.5);
temp = a + (b-a).*rand(N,1);
xi = norminv(temp,mxi,vxi^0.5);

%--------------------------------------------------------------------------
% compute conditional moments of eta distribution
%--------------------------------------------------------------------------
%compute conditional moments of epsilon distribution
A = [Rho_nuet*(Signu2*Siget2)^0.5;...
     Rho_xiet*(Sigxi2*Siget2)^0.5;...
     Rho_epet*(Sigep2*Siget2)^0.5];
B = [Signu2 Rho_nuxi*(Signu2*Sigxi2)^0.5 Rho_nuep*(Signu2*Sigep2)^0.5;...
     Rho_nuxi*(Signu2*Sigxi2)^0.5 Sigxi2 Rho_xiep*(Sigxi2*Sigep2)^0.5;...
     Rho_nuep*(Signu2*Sigep2)^0.5 Rho_xiep*(Sigxi2*Sigep2)^0.5 Sigep2];
C = [nu'; xi'; ep'];
met = (A'*(B\C))';
vet = max(1e-9,Siget2 - A'*(B\A));
XGamma  = [X3 afK]*Gamma;

%--------------------------------------------------------------------------
% compute event probabilities and likelihood function
%--------------------------------------------------------------------------
%compute purchase probabilities 
pr_pric = 1e-99 + normpdf(pK-XLambda,0,Signu2^0.5);
pr_less = 1e-99 + normcdf(dK-XBeta-mep,0,vep^0.5);
pr_more = 1e-99 + normpdf(tK-XBeta-mep,0,vep^0.5);
pr_sale = 1e-99 + normcdf(XDelta+mxi,0,vxi^0.5);
pr_nosl = 1e-99 + (1-normcdf([X4 pK dK dK2]*Delta+mxi,0,vxi^0.5));

%compute fracpaid probabilities
fraclo  = floor(fracpaid*G)/G;
frachi  = fraclo + 1/G;
pr_frac = 1e-99 + normcdf(log(frachi)-XGamma-met,0,vet^0.5)-normcdf(log(fraclo)-XGamma-met,0,vet^0.5); %observed fracpaid
pr_none = 1e-99 + normcdf(log(1/(G*4))-XGamma-met,0,vet^0.5); %no payments
pr_paid = 1e-99 + normcdf(-log(cenpoint)+XGamma+met,0,vet^0.5); %censored

%compute likelihood function
logL = log(pr_pric) + I.*log(pr_sale) + (1-I).*log(pr_nosl) + less.*log(pr_less) + more.*log(pr_more) + obsdef.*log(pr_frac) + none.*log(pr_none) + cenind.*log(pr_paid);
logL(logL==-Inf) = -999;
logL(isnan(logL))= -999;
loglike = - sum(logL);

%--------------------------------------------------------------------------
% output parameter and likelihood function values
%--------------------------------------------------------------------------
fid = fopen('likefull.txt','wt');
fprintf(fid,'L = %6.3f\n',loglike);
fprintf(fid,'pr(sale) = %6.3f\n',mean(pr_sale));
fprintf(fid,'%6.3f\n',par);
fclose(fid);

%--------------------------------------------------------------------------
% end of function
%--------------------------------------------------------------------------
