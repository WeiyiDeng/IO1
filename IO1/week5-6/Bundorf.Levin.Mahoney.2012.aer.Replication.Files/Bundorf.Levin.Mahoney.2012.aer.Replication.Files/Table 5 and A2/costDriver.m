clear all; clc;

global W M bidData costData...
     zBid zCost eBid eCost dropped; %#ok<NUSED>

load structureBidDataDrop;
load structureCostDataDrop;
load costOut

theta0 = thetaHat;

W = inv(ST);

%set options for objective function 
options = optimset('GradObj','on','DerivativeCheck', 'off', 'TolX', 10^(-6),'TolFun',10^(-6),...
    'MaxFunEvals', 25, 'Display','iter');

[thetaHat t exitflag] = fminunc(@costMoment,theta0, options);

tBid = length(eBid);
sBid = 0;

for i=1:tBid
    ZZ = zBid(i,:)'*zBid(i,:);
    ee = eBid(i)^2;
    sBid = sBid + ee*ZZ;
end

SBidT = sBid/tBid;
SBidTT = SBidT/tBid;

tCost = length(eCost);
sCost = 0;

for i=1:tCost
    ZZ = zCost(i,:)'*zCost(i,:);
    ee = eCost(i)^2;
    sCost = sCost + ee*ZZ;
end

SCostT = sCost/tCost;
SCostTT = SCostT/tCost;

V1 = inv(M'*W*M);

ST = zeros(length(thetaHat));
STT = zeros(length(thetaHat));

ST(1:12,1:12) = SBidT;
ST(13:24,13:24) = SCostT;

STT(1:12,1:12) = SBidTT;
STT(13:24,13:24) = SCostTT;

V2 = M'*W*STT*W*M;

V = V1*V2*V1;

se = sqrt(diag(V));

save costOut thetaHat W M  eBid zBid eCost zCost ...
     V se ST;


