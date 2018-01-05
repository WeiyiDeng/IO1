%--------------------------------------------------------------------------
% simulate_unobservables.m
%--------------------------------------------------------------------------

function [nu xi ep et] = simulate_unobservables(par)

global N X1 X2 X3 X4 lK pK dK dK2 tK afK fracpaid;
global I less more none cenind obsdef cenpoint G;
global mxi vxi mep vep met vet;

%define parameters
[Lambda Delta Beta Gamma Signu2 Sigxi2 Sigep2 Siget2 Rho_nuxi Rho_nuep Rho_nuet Rho_xiep Rho_xiet Rho_epet] = define_parameters(par);
XLambda = [X1 lK]*Lambda; nu = pK - XLambda;

%--------------------------------------------------------------------------
% simulate down unobservable (ep), conditional on nu
%--------------------------------------------------------------------------
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
% simulate repayment unobservables (et), conditional on nu, xi, and ep
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

%determine truncation points for less than minimum draws (1=low,2=high)
XGamma  = [X3 afK]*Gamma;
trunc1 = cenind.*(log(cenpoint)-XGamma); trunc1(~cenind) = -999;
trunc2 = none.*(log(1/(G*4))-XGamma) + ~none.*999;

%draw unobservables
stream = RandStream('mt19937ar','Seed',3);
RandStream.setDefaultStream(stream);
a = normcdf(trunc1,met,vet^0.5);
b = normcdf(trunc2,met,vet^0.5);
temp = a + (b-a).*rand(N,1);
draw = norminv(temp,met,vet^0.5); draw(draw==Inf) = 1-1e-9;
logf = log(fracpaid); logf(isnan(logf))=-998; logf(logf==-Inf)=-999;
et = obsdef.*(logf-XGamma) + ~obsdef.*draw;

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------
