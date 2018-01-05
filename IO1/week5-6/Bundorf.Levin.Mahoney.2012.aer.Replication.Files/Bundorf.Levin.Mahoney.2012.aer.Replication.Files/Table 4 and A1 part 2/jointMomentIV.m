function [t,grad] = jointMomentIV(theta)

global W M;

% clear all; clc;
% 
% global choiceData %#ok<NUSED>
% 
% load structureChoiceData;
% 
% W = eye(28);
% 
% theta = [thetaChoice thetaBC]; %#ok<USENS>

[mChoice MChoice] = jointChoiceIV_10(theta);

M = zeros(28,28);

m = mChoice;

M(1:28,1:28) = MChoice;

t = m'*W*m;

grad = 2*M'*W*m;



