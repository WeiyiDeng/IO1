clear all; clc;

load structureWelfareData;

gridmax = 2;
step = 1;

priceVector = -gridmax:step:gridmax;
gridsize = length(priceVector);
N = length(choiceData);

for n = 1:N
    
    if mod(n,100) == 0
        disp(n);
    end
    
    fam = choiceData(n).fam(1);
    
    %generate a mapping from rows to plans

    L = length(choiceData(n).group);  

    planIndex = zeros(L,1);

    for l = 1:L

        planIndex(l) = find(choiceData(n).planX(l,:)==1);

    end
    
    opt(n).year = choiceData(n).year(1); %#ok<AGROW>
    
    opt(n).group = choiceData(n).group(1); %#ok<AGROW>
    
    opt(n).tier = choiceData(n).tier(1); %#ok<AGROW>

    for i = 1:gridsize
        for j = 1:gridsize
            for k = 1:gridsize

                price = zeros(L,1);

                for l = 1:L

                    if planIndex(l) == 1
                        price(l) = priceVector(i);
                    elseif planIndex(l) == 2
                        price(l) = priceVector(j);
                    elseif planIndex(l) == 4
                        price(l) = priceVector(k);
                    end

                end

                %plan characteristics

                phi = [price choiceData(n).co' choiceData(n).ded' choiceData(n).ph'];

                % plan dummies and interactions with pro, pro90, age, male and inc

                proMu = choiceData(n).proX; %+ rho*choiceData(n).mu(:,:,ns);

                xi = [choiceData(n).planX proMu choiceData(n).pro90X...
                    choiceData(n).ageX choiceData(n).maleX choiceData(n).incX ];

                %indirect utility

                v = phi*betaPhi' + xi*betaXi' + choiceData(n).epsilon;

                [value choice] = max(v); 

                ss = v(choice) + price(choice) - choiceData(n).cost(choice);

                opt(n).ss(i,j,k) = ss*fam; %#ok<AGROW>

            end
        end
    end
end

save welfareOpt priceVector gridsize opt;



