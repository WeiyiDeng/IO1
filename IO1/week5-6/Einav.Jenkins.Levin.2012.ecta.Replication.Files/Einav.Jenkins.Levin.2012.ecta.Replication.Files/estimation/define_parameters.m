%--------------------------------------------------------------------------
% define_parameters.m
%--------------------------------------------------------------------------

function [Lambda Delta Beta Gamma Signu2 Sigxi2 Sigep2 Siget2 Rho_nuxi Rho_nuep Rho_nuet Rho_xiep Rho_xiet Rho_epet] = define_parameters(par)

global K1 K2 K3 K4;

% define parameters (copied from SimulatedLikelihood.m)
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

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------
