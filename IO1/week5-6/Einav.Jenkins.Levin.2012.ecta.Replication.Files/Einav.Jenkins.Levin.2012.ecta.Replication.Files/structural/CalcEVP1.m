%--------------------------------------------------------------------------
% Inputs vector of liquidity and car value indices, outputs future values
%--------------------------------------------------------------------------

function [output] = CalcEVP1(vi,yi,evp1)

global Nsim lnum;

vii = vi*ones(1,lnum);                                  %NxL
yii = yi*ones(1,lnum);                                  %NxL
lii = ones(Nsim,1)*(1:lnum);                            %NxL
index = sub2ind(size(evp1),yii,lii,vii);
output = evp1(index);                                   %NxL

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------
