%--------------------------------------------------------------------------
% simulate_for_model_fit.m
%--------------------------------------------------------------------------

function [nu xi ep et] = simulate_for_model_fit(par)

global N X1 lK pK;
global mxi vxi mep vep met vet;

%define parameters
[Lambda Delta Beta Gamma Signu2 Sigxi2 Sigeps2 Sigeta2 Rho_nuxi Rho_nueps Rho_nueta Rho_xieps Rho_xieta Rho_epseta] = define_parameters(par);
XLambda = [X1 lK]*Lambda; nu = pK - XLambda;

%--------------------------------------------------------------------------
% simulate sale unobservables (xi), conditional on nu
%--------------------------------------------------------------------------
%compute conditional moments of xi distribution
A = Rho_nuxi*(Signu2*Sigxi2)^0.5;
B = Signu2;
C = nu';
mxi = (A'*inv(B)*C)';
vxi = max(1e-9,Sigxi2 - A'*inv(B)*A);

%draw unobservables
stream = RandStream('mt19937ar','Seed',1);
RandStream.setDefaultStream(stream);
xi = mxi + (vxi^0.5).*randn(N,1);
%rand('state',1);
%xi = random('normal',mxi,vxi^0.5);

%--------------------------------------------------------------------------
% simulate down unobservables (ep), conditional on nu and xi
%--------------------------------------------------------------------------
%compute conditional moments of epsilon distribution
A = [Rho_nueps*(Signu2*Sigeps2)^0.5; Rho_xieps*(Sigxi2*Sigeps2)^0.5];
B = [Signu2 Rho_nuxi*(Signu2*Sigxi2)^0.5; Rho_nuxi*(Signu2*Sigxi2)^0.5 Sigxi2];
C = [nu'; xi'];
mep = (A'*inv(B)*C)';
vep = max(1e-9,Sigeps2 - A'*inv(B)*A);

%draw unobservables
stream = RandStream('mt19937ar','Seed',2);
RandStream.setDefaultStream(stream);
ep = mep + (vep^0.5).*randn(N,1);

%--------------------------------------------------------------------------
% simulate repayment unobservables (et), conditional on nu, xi, and ep
%--------------------------------------------------------------------------
%compute conditional moments of epsilon distribution
A = [Rho_nueta *(Signu2 *Sigeta2)^0.5;...
     Rho_xieta *(Sigxi2 *Sigeta2)^0.5;...
     Rho_epseta*(Sigeps2*Sigeta2)^0.5];
B = [Signu2 Rho_nuxi*(Signu2*Sigxi2)^0.5 Rho_nueps*(Signu2*Sigeps2)^0.5;...
     Rho_nuxi*(Signu2*Sigxi2)^0.5 Sigxi2 Rho_xieps*(Sigxi2*Sigeps2)^0.5;...
     Rho_nueps*(Signu2*Sigeps2)^0.5 Rho_xieps*(Sigxi2*Sigeps2)^0.5 Sigeps2];
C = [nu'; xi'; ep'];
met = (A'*inv(B)*C)';
vet = max(1e-9,Sigeta2 - A'*inv(B)*A);

%draw unobservables
stream = RandStream('mt19937ar','Seed',3);
RandStream.setDefaultStream(stream);
et = met + (vet^0.5).*randn(N,1);

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------
