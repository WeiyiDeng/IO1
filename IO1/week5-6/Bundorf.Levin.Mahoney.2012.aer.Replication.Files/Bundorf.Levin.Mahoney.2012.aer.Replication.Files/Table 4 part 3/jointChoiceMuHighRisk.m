function [mChoice, MChoice] = jointChoiceMuHighRisk(theta)

global NS choiceData eChoice zChoice;

% clear all; clc;
% load structureChoiceData;
% 
% NS = 1;
% 
% theta = .001*ones(1,39);

betaPhi = theta(25:28);

betaXi = theta(1:24);

rho = theta(39);

%the code loops over each individual in the data set.

%for each individaul the code loops over the choice set to determine the
%denominator of the logit probability: denominator = sum of exp(Xbeta)

%the code loops over the choices again to determine the choice probability
%pr = exp(Xbeta) / denominator

%the moment condidition is g = (dChoice-pr)'*z, where z are the exogenous
%variables

N = length(choiceData); %20; %%%#ok<USENS>
t = 13281; %77; 
K = 28;
mChoice = zeros(K,NS);
MChoice = zeros(K,K+1,NS);
eChoice = zeros(t,NS);
zChoice = zeros(t,K,NS);

for ns = 1:NS
    
    e = zeros(N,1);
    z = zeros(N,K);
    de = zeros(N,K+1);
    k = 1;

    for n = 1:N

        I = length(choiceData(n).group);  %#ok<NASGU>

        %plan characteristics

        phi = [choiceData(n).premium' choiceData(n).co' choiceData(n).ded' choiceData(n).ph'];

        %plan dummies and interactions with pro, pro90, age, male and inc

        proMu = choiceData(n).proX + rho*choiceData(n).mu(:,:,ns);
        
        maxPro = choiceData(n).maxPro(1) + rho*choiceData(n).maxMu(1,ns);
        
        if maxPro>2.25
            pro90 = 1;
        else
            pro90 = 0;
        end
        
        pro90Mu = choiceData(n).planX*pro90;

        xi = [choiceData(n).planX proMu pro90Mu ...
            choiceData(n).ageX choiceData(n).maleX choiceData(n).incX ];

        %xBeta

        xBeta = phi*betaPhi' + xi*betaXi';

        numerator = exp(xBeta);

        denominator = sum(numerator);

        %calcualte the logit choice probability

        pr = numerator/denominator;

        e(k:(k+I-1),1) = pr - choiceData(n).dchoice';

        %construct instrument vector z

        z(k:(k+I-1),:) = [xi phi];

        %construct gradient input matrix dMDBeta
        
       dXBetadrho = xi(:,1:4)*betaXi(5:8)'.*sum(choiceData(n).mu(:,:,ns),2);

        for i = 1:I

            A = denominator*numerator(i)*[xi(i,:) phi(i,:) dXBetadrho(i)];
            B = numerator(i)*numerator'*[xi phi dXBetadrho];
            C = denominator*denominator;

            de(k+i-1,:) = (A-B)/C;

        end

        k = k+I;

    end

    J = length(e);

    mChoice(:,ns) = z'*e/J;

    MChoice(:,:,ns) = z'*de/J;

    zChoice(:,:,ns) = z;

    eChoice(:,ns) = e;

end







