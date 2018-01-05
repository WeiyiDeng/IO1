clear all; clc;

load structureChoiceData;
load jointOut;

N = length(choiceData); %#ok<USENS>

betaPhi = thetaHat(25:28);
betaXi = thetaHat(1:24);
alpha = thetaHat(25);

count = zeros(4,1);
phi = zeros(4,4);
xi = zeros(4,6);
d = zeros(4,1);

phiBar = zeros(4,4);
xiBar = zeros(4,24);
mfx = zeros(4,1);
mfxSE = zeros(4,1);

for n = 1:N

    if choiceData(n).fam(1)>1
        continue
    end
    
    I = length(choiceData(n).group);  %#ok<NASGU>

    for i = 1:I

        %plan characteristics

        phi(choiceData(n).plan(i),:) =  phi(choiceData(n).plan(i),:) +...
            [choiceData(n).premium(i) choiceData(n).co(i) choiceData(n).ded(i) choiceData(n).ph(i)];

        %plan dummies and interactions with pro, pro90, age, male and inc

        xi(choiceData(n).plan(i),:) = xi(choiceData(n).plan(i),:) + ...
            [choiceData(n).planX(i,choiceData(n).plan(i)) choiceData(n).proX(i,choiceData(n).plan(i))...
            choiceData(n).pro90X(i,choiceData(n).plan(i)) choiceData(n).ageX(i,choiceData(n).plan(i)) ...
            choiceData(n).maleX(i,choiceData(n).plan(i)) choiceData(n).incX(i,choiceData(n).plan(i))];

        d(choiceData(n).plan(i),1) = d(choiceData(n).plan(i),1) + ...
            choiceData(n).dchoice(i);
        
        count(choiceData(n).plan(i),1) = count(choiceData(n).plan(i),1) + 1;

    end

end


for i = 1:4

    phiBar(i,:) = phi(i,:)/count(i);
    xi(i,:) = xi(i,:)/count(i); %#ok<AGROW>

end

for i=1:4
    xiBar(i,:) = [xi(1,:) xi(2,:) xi(3,:) xi(4,:)];
end

xBeta = phiBar*betaPhi' + xiBar*betaXi';

numerator = exp(xBeta);

denominator = sum(numerator);

for i = 1:4

    A = denominator*numerator(i)*alpha;
    B = numerator(i)*numerator(i)*alpha;
    C = denominator*denominator;

    mfx(i) = (A-B)/C;
    mfxSE(i) = (A-B)/(alpha*C)*se(19);

end

%mfx as is is percentage point change in marketshare for a $100 increase in
%the monthly price

share = mean(d./count);

mfxAdjust = mfx/(12*share);

mfxAdjustSE = se/(12*share);
