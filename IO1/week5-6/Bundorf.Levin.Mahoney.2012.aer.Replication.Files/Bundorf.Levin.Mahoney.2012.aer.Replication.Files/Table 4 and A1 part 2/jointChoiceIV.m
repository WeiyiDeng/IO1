function [mChoice, MChoice] = jointChoiceIV(theta)

global choiceData eChoice zChoice dropped;

% clear all; clc;
% load structureChoiceData;
% %load jointOut;
% 
% theta = .001*ones(1,28);
% W = eye(28);

betaPhi = theta(25:28);
 
betaXi = theta(1:24);

%the code loops over each individual in the data set.

%for each individaul the code loops over the choice set to determine the
%denominator of the logit probability: denominator = sum of exp(Xbeta)

%the code loops over the choices again to determine the choice probability
%pr = exp(Xbeta) / denominator

%the moment condidition is g = (dChoice-pr)'*z, where z are the exogenous
%variables

%the function returns t = g*I*g'

N = length(choiceData); %#ok<USENS>
K = 28;
e = zeros(N,1);
z = zeros(N,K);
de = zeros(N,K);
k = 1;
dropped = 0;

for n = 1:N
    
    if max(choiceData(n).bidHatDemAT)==99.99
        
        dropped = dropped + 1;
        
        continue
        
    end

    I = length(choiceData(n).group);  %#ok<NASGU>
   
    %plan characteristics

    phi = [choiceData(n).premiumAT' choiceData(n).co' choiceData(n).ded' choiceData(n).ph'];

    %instruments
    
    phiIV = [choiceData(n).bidHatDemATAdjust' phi(:,2:4)];
    
    %plan dummies and interactions with pro, pro90, age, male and inc

    xi = [choiceData(n).planX choiceData(n).proX choiceData(n).pro90X...
        choiceData(n).ageX choiceData(n).maleX choiceData(n).incX];

    %xBeta

    xBeta = phi*betaPhi' + xi*betaXi';

    numerator = exp(xBeta);

    denominator = sum(numerator);

    %calcualte the logit choice probability
    
    pr = numerator/denominator;

    e(k:(k+I-1),1) = pr - choiceData(n).dchoice';
    
    %construct instrument vector z

    z(k:(k+I-1),:) = [xi phiIV];

    %construct gradient input matrix dMDBeta
    
    for i = 1:I

        A = denominator*numerator(i)*[xi(i,:) phi(i,:) ];
        B = numerator(i)*numerator'*[xi phi];
        C = denominator*denominator;
        
        de(k+i-1,:) = (A-B)/C;
        
    end    

    k = k+I;
        
end

J = length(e);

mChoice = z'*e/J;

MChoice = z'*de/J;

zChoice = z;

eChoice = e;







