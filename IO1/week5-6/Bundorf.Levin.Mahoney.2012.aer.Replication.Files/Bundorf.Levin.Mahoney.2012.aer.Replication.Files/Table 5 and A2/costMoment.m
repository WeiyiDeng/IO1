function [t,grad] = costMoment(theta)

global W M;

% clear all; clc;
% 
% global bidData costData choiceData %#ok<NUSED>
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
% theta = [thetaChoice thetaBC]; %#ok<USENS>

[mBid MBid] = costBid(theta);

[mCost MCost] = costCost(theta);

M = zeros(24:14);

m = [mBid; mCost];

M(1:12,1:14) = MBid;
M(13:24,1:14) = MCost;

t = m'*W*m;

grad = 2*M'*W*m;



