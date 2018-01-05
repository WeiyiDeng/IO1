function [t,grad] = jointMomentMuHighRisk(theta)

global W M NS;

% clear all; clc;
% 
% global bidData costData choiceData NS %#ok<NUSED>
% 
% load structureChoiceData;
% load structureBidData;
% load structureCostData;
% 
% load jointOutChoice;
% load jointOutBC;
% 
% W = eye(44);
% 
% NS = 1;
% 
% theta = [thetaChoice thetaBC 1]; %#ok<USENS>

[mChoice MChoice] = jointChoiceMuHighRisk(theta);

[mBid MBid] = jointBidMuHighRisk(theta);

[mCost MCost] = jointCostMuHighRisk(theta);

tVector = zeros(1,NS);
gradVector = zeros(39,NS);

M = zeros(44,39,NS);

for ns = 1:NS
    
    m = [mChoice(:,ns); mBid; mCost(:,ns)];

    M(1:28,1:28,ns) = MChoice(1:28,1:28,ns);
    M(1:28,39,ns) = MChoice(1:28,29,ns);
    M(29:36,29:38,ns) = MBid;
    M(37:44,29:39,ns) = MCost(:,:,ns);

    tVector(ns) = m'*W*m;

    gradVector(:,ns) = 2*M(:,:,ns)'*W*m;

end

t = mean(tVector);

grad = mean(gradVector,2);

