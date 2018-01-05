clear all; clc;

global W M choiceData ...
    zChoice eChoice dropped; %#ok<NUSED>

load structureChoiceDataIV;

load jointOut

theta0 = thetaHat; %#ok<NASGU>

W = inv(ST);

%set options for objective function 
options = optimset('GradObj','on','DerivativeCheck', 'on', 'TolX', 10^(-6),'TolFun',10^(-6),...
    'MaxFunEvals', 1, 'Display','iter');

[thetaHat t exitflag] = fminunc(@jointMomentIV_10,theta0, options);

save jointOut_10 thetaHat W M eChoice zChoice;

tH1 = [thetaHat(1:2) thetaHat(4)] - thetaHat(3);
tH2 = [thetaHat(5:6) thetaHat(8)] - thetaHat(7);
tH3 = [thetaHat(9:10) thetaHat(12)] - thetaHat(11);
tH4 = [thetaHat(13:14) thetaHat(16)] - thetaHat(15);
tH5 = [thetaHat(17:18) thetaHat(20)] - thetaHat(19);
tH6 = [thetaHat(21:22) thetaHat(24)] - thetaHat(23);

beta = [tH1 tH2 tH3 tH4 tH5 tH6 thetaHat(25:28)];

for i=1:28
    mO1(i,:) = [M(i,1:2) M(i,4)] - M(i,3); %#ok<AGROW>
    mO2(i,:) = [M(i,5:6) M(i,8)] - M(i,7); %#ok<AGROW>
    mO3(i,:) = [M(i,9:10) M(i,12)] - M(i,11); %#ok<AGROW>
    mO4(i,:) = [M(i,13:14) M(i,16)] - M(i,15); %#ok<AGROW>
    mO5(i,:) = [M(i,17:18) M(i,20)] - M(i,19); %#ok<AGROW>
    mO6(i,:) = [M(i,21:22) M(i,24)] - M(i,23); %#ok<AGROW>
end

Momit = [mO1 mO2 mO3 mO4 mO5 mO6 M(:,25:28)];

t = length(eChoice);
s = 0;

for i=1:t
    ZZ = zChoice(i,:)'*zChoice(i,:);
    esquared = eChoice(i)^2;
    s = s + esquared*ZZ;
end

SChoiceT = s/(t-dropped);
SChoiceTT = SChoiceT/(t-dropped);

A = Momit(1:28,1:22);
Around = round(A*10000);
WA = W(1:28,1:28);
VChoice = inv(Around'*WA*Around)*100000000;

V1 = VChoice;

ST = SChoiceT;
STT = SChoiceTT;

V2 = Momit'*W*STT*W*Momit;

V = V1*V2*V1;

se = sqrt(diag(V));

save jointOut thetaHat W M eChoice zChoice ...
    beta V se ST;


