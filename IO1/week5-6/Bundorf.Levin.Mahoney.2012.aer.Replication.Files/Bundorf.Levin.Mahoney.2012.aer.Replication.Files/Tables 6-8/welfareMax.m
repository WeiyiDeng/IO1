% this code determines allocations and the prices that decentralize these
% allocations under the following nondiscrimination constraints
% 1 observed
% 2 efficient (price - cost);
% 3 uninform
% 4 uniform by firm
% 5 uniform by tier
% 6 uniform by tier by firm

% the data are saved in the vectors
% choiceData.choiceSim
% choiceData.price

clear all; clc;

load structureWelfareData;
load welfareOpt;

% define variables

N = length(choiceData);

uniform = zeros(gridsize,gridsize,gridsize);

for g = 1:14
    group(g).mat = zeros(gridsize,gridsize,gridsize); %#ok<AGROW>
end

for t = 1:6
    tier(t).mat = zeros(gridsize,gridsize,gridsize); %#ok<AGROW>
end

for g = 1:14
    for t = 1:6
        grouptier(g,t).mat = zeros(gridsize,gridsize,gridsize); %#ok<AGROW>
    end
end


for n = 1:N

    %*************************************************************************
    %generate gridsize matrixes of social surplus

    uniform = uniform + opt(n).ss;

    group(choiceData(n).group(1)).mat = group(choiceData(n).group(1)).mat + opt(n).ss; %#ok<AGROW>

    tier(choiceData(n).tier(1)).mat = tier(choiceData(n).tier(1)).mat + opt(n).ss; %#ok<AGROW>

    grouptier(choiceData(n).group(1),choiceData(n).tier(1)).mat = ...
        grouptier(choiceData(n).group(1), choiceData(n).tier(1)).mat + opt(n).ss; %#ok<AGROW>

end

% find social surplus maximizing element of each gridsize by gridsize by
% gridsize matrix

price.uniform = zeros(1,4);
price.group = zeros(14,4);
price.tier = zeros(6,4);
price.grouptier =zeros(14,6,4);

maximum.uniform = 0;
maximum.group = zeros(14,1);
maximum.tier = zeros(6,1);
maximum.grouptier = zeros(14,6);

for i = 1:gridsize
    for j = 1:gridsize
        for k = 1:gridsize

            if uniform(i,j,k) > maximum.uniform
                maximum.uniform = uniform(i,j,k);
                price.uniform = [priceVector(i) priceVector(j) 0 priceVector(k)];
            end

            for g = 1:14
                if group(g).mat(i,j,k) > maximum.group(g)
                    maximum.group(g) = group(g).mat(i,j,k);
                    price.group(g,:) = [priceVector(i) priceVector(j) 0 priceVector(k)];
                end
            end

            for t = 1:6
                if tier(t).mat(i,j,k) > maximum.tier(t)
                    maximum.tier(t) = tier(t).mat(i,j,k);
                    price.tier(t,:) = [priceVector(i) priceVector(j) 0 priceVector(k)];
                end
            end

            for g = 1:14
                for t = 1:6
                    if grouptier(g,t).mat(i,j,k) > maximum.grouptier(g,t)
                        maximum.grouptier(g,t) = grouptier(g,t).mat(i,j,k);
                        price.grouptier(g,t,:) = [priceVector(i) priceVector(j) 0 priceVector(k)];
                    end
                end
            end


        end
    end
end

% assign prices and allocations to choiceData.price and
% choiceData.choiceSim respectively

for n = 1:N

    I = length(choiceData(n).group);

    planIndex = zeros(I,1);

    for i = 1:I

        planIndex(i) = find(choiceData(n).planX(i,:)==1);

    end

    %*****************************************************************
    %efficient prices and allocations

    choiceData(n).price(:,2) = choiceData(n).cost'; %#ok<AGROW>

    %optimal prices and allocations

    for i = 1:I

        choiceData(n).price(i,3) = price.uniform(planIndex(i)); %#ok<AGROW>

        choiceData(n).price(i,4) = price.group(choiceData(n).group(1),planIndex(i)); %#ok<AGROW>

        choiceData(n).price(i,5) = price.tier(choiceData(n).tier(1),planIndex(i)); %#ok<AGROW>

        choiceData(n).price(i,6) = price.grouptier(choiceData(n).group(1),choiceData(n).tier(1),planIndex(i)); %#ok<AGROW>

    end

    %*****************************************************************
    %observed prices and allocations

    choiceData(n).price(:,1) = choiceData(n).premium; %#ok<AGROW>

    [valuel choiceSim] = max(choiceData(n).dchoice);

    choiceData(n).choiceSim(1,1) = choiceSim; %#ok<AGROW>

    %****************************************************************
    %optimal prices and allocations

    for sim = 2:6

        %plan characteristics

        phi = [choiceData(n).price(:,sim) choiceData(n).co' choiceData(n).ded' choiceData(n).ph'];

        % plan dummies and interactions with pro, pro90, age, male and inc

        proMu = choiceData(n).proX; %+ rho*choiceData(n).mu(:,:,ns);

        xi = [choiceData(n).planX proMu choiceData(n).pro90X...
            choiceData(n).ageX choiceData(n).maleX choiceData(n).incX ];

        %indirect utility

        v = phi*betaPhi' + xi*betaXi' + choiceData(n).epsilon; %#ok<USENS>

        [value choiceSim] = max(v);

        choiceData(n).choiceSim(1,sim) = choiceSim; %#ok<AGROW>

    end

end

save structureWelfareData choiceData betaPhi betaXi;


