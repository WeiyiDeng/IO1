%% fixedPoint.m

function [capU0,capU1] = fixedPoint_CE2(u0,u1,capPi,rho,tolFixedPoint,bellman,capU0,capU1)

rand_ind = rand(1,1);
rng(rand_ind)

nSuppX = size(capPi,1);

if isempty(capU0) 
% 	capU0 = zeros(nSuppX,2);
%     rng(rand(1,1))
    capU0 = rand(nSuppX,2);                   % w: changed here, still makes no difference?
end
if isempty(capU1)
% 	capU1 = zeros(nSuppX,2);
    capU1 = rand(nSuppX,2);                   % w: changed here, still makes no difference?
end

inU0 = capU0+2*tolFixedPoint;
inU1 = capU1+2*tolFixedPoint;
while (max(max(abs(inU0-capU0)))>tolFixedPoint) || (max(max(abs(inU1-capU1)))>tolFixedPoint);
	inU0 = capU0;
	inU1 = capU1;
	[capU0,capU1] = bellman(inU0,inU1,u0,u1,capPi,rho);
end
