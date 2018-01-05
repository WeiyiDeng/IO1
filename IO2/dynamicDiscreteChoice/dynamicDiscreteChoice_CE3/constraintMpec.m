function [diff, outeq] = constraintMpec(parameter,delta, rho, capPi,supportX)

beta = parameter(1:2);
delta = [delta(1);parameter(3)];

nSuppX = length(supportX);
capU0 = reshape(parameter(4:13),nSuppX,2);
capU1 = reshape(parameter(14:23),nSuppX,2);
inU0 = capU0;
inU1 = capU1;
[u0,u1] = flowpayoffs(supportX,beta,delta);
[capU0,capU1] = bellman(inU0,inU1,u0,u1,capPi,rho);

diff = max(max(abs([inU0;inU1]-[capU0;capU1])));
outeq = [];
% diff = [inU0;inU1]-[capU0;capU1];
