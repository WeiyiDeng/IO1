%--------------------------------------------------------------------------
% Inputs vector of liquidity and car value indices, outputs future values
%--------------------------------------------------------------------------

function [output] = CalcEVPt(li,vi,yit,X) %X = EVP, EVD, or D

global T;

lii = li*ones(1,T);
vii = vi*ones(1,T);
yii = yit;
tii = ones(size(yit,1),1)*(1:T);
index = sub2ind(size(X),yii,lii,vii,tii);
output = X(index);

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------
