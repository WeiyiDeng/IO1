clear all; clc;

global NS W M choiceData bidData costData...
    zChoice zBid zCost eChoice eBid eCost; %#ok<NUSED>

load structureChoiceData;
load structureBidData;
load structureCostData;
load jointOutMu20;

NS = 20;

W = inv(ST);

theta0 = thetaHat; %#ok<USENS>

%set options for function 
options = optimset('GradObj','on','DerivativeCheck', 'off', 'TolX', 10^(-4),'TolFun',10^(-4),...
    'MaxFunEvals', 5, 'Display','iter');

[thetaHat t exitflag] = fminunc(@jointMomentMuHighRisk,theta0, options);

% coefficients on plan 3 are not identified
% subtract coefficients on plan 3 from other plans

tH1 = [thetaHat(1:2) thetaHat(4)] - thetaHat(3);
tH2 = [thetaHat(5:6) thetaHat(8)] - thetaHat(7);
tH3 = [thetaHat(9:10) thetaHat(12)] - thetaHat(11);
tH4 = [thetaHat(13:14) thetaHat(16)] - thetaHat(15);
tH5 = [thetaHat(17:18) thetaHat(20)] - thetaHat(19);
tH6 = [thetaHat(21:22) thetaHat(24)] - thetaHat(23);

beta = [tH1 tH2 tH3 tH4 tH5 tH6 thetaHat(25:39)];

% average across draws

M = mean(M,3);
eChoice = mean(eChoice,2);
eBid = mean(eBid,2);
eCost = mean(eCost,2);
zChoice = mean(zChoice,3);
zCost = mean(zCost,3);
zBid = mean(zBid,3);

for i=1:44
    mO1(i,:) = [M(i,1:2) M(i,4)] - M(i,3); %#ok<AGROW>
    mO2(i,:) = [M(i,5:6) M(i,8)] - M(i,7); %#ok<AGROW>
    mO3(i,:) = [M(i,9:10) M(i,12)] - M(i,11); %#ok<AGROW>
    mO4(i,:) = [M(i,13:14) M(i,16)] - M(i,15); %#ok<AGROW>
    mO5(i,:) = [M(i,17:18) M(i,20)] - M(i,19); %#ok<AGROW>
    mO6(i,:) = [M(i,21:22) M(i,24)] - M(i,23); %#ok<AGROW>
end

Momit = [mO1 mO2 mO3 mO4 mO5 mO6 M(:,25:39)];

t = length(eChoice);
s = 0;

for i=1:t
    ZZ = zChoice(i,:)'*zChoice(i,:);
    esquared = (1)*eChoice(i)^2;
    s = s + esquared*ZZ;
end

SChoiceT = s/t;
SChoiceTT = SChoiceT/t;

A = Momit(1:28,1:22);
Around = round(A*10000);
WA = W(1:28,1:28);
VChoice = inv(Around'*WA*Around)*100000000;

tBid = length(eBid);
sBid = 0;

for i=1:tBid
    ZZ = zBid(i,:)'*zBid(i,:);
    ee = (1)*eBid(i)^2;
    sBid = sBid + ee*ZZ;
end

SBidT = sBid/tBid;
SBidTT = SBidT/tBid;

tCost = length(eCost);
sCost = 0;

for i=1:tCost
    ZZ = zCost(i,:)'*zCost(i,:);
    ee = (1)*eCost(i)^2;
    sCost = sCost + ee*ZZ;
end

SCostT = sCost/tCost;
SCostTT = SCostT/tCost;

MBC = M(29:44,29:39);

WBC = W(29:44,29:44);

VBC = inv(MBC'*WBC*MBC);

V1 = zeros(33);

V1(1:22,1:22) = VChoice;
V1(23:33,23:33) = VBC;

ST = zeros(length(thetaHat));
STT = zeros(length(thetaHat));

ST(1:28,1:28) = SChoiceT;
ST(29:36,29:36) = SBidT;
ST(37:44,37:44) = SCostT;

STT(1:28,1:28) = SChoiceTT;
STT(29:36,29:36) = SBidTT;
STT(37:44,37:44) = SCostTT;

V2 = Momit'*W*STT*W*Momit;

V = V1*V2*V1;

se = sqrt(diag(V));

save jointOutMuNoRisk thetaHat W M eChoice zChoice eBid zBid eCost zCost ...
    beta V se ST;
