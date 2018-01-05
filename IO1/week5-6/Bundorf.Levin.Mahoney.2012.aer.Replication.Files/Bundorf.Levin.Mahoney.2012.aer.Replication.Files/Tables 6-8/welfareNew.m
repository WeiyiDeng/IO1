% this script calculates welfare from single plan offerings
% 7 is AHMO
% 8 is APPO
% 9 is BHMO
% 10 is BPOS
% 11 is optimal single plan by firm

clear all; clc;

load structureWelfareData;

%zero vectors for each firm i and plan j

for i = 1:12
    for j = 1:2
        new(i,j).ss = zeros(4,1);%#ok<AGROW>
    end
end

% for each indvidual

N = length(choiceData);

index = 0;

for n = 1:N

    I = length(choiceData(n).group);  

    %planIndex maps from the row in the data to the plan's index

    planIndex = zeros(I,1);

    for i = 1:I

        planIndex(i) = find(choiceData(n).planX(i,:)==1);

    end

    % define group and year

    group = choiceData(n).group(1);

    year = choiceData(n).year(1)-2003;

    fam = choiceData(n).fam(1);

    %plan characteristics

    phi = [choiceData(n).cost' choiceData(n).co' choiceData(n).ded' choiceData(n).ph'];

    %plan dummies and interactions with pro, pro90, age, male and inc

    xi = [choiceData(n).planX choiceData(n).proX choiceData(n).pro90X...
        choiceData(n).ageX choiceData(n).maleX choiceData(n).incX ];

    %calculate social surplus

    ss = phi*betaPhi' + xi*betaXi' + choiceData(n).epsilon; %#ok<USENS>

    new(group,year).offered = planIndex; %#ok<AGROW>

    for i = 1:I

        new(group,year).ss(planIndex(i)) = new(group,year).ss(planIndex(i)) + ss(i)*fam; %#ok<AGROW>

    end

end

%find optimal plan offering

for i=1:12
    for j=1:2

        ss = - 9999;

        for k = 1:length(new(i,j).offered)

            if new(group,year).ss(new(i,j).offered(k)) > ss

                ss = new(group,year).ss(new(i,j).offered(k));

                new(i,j).plan = new(i,j).offered(k); %#ok<AGROW>

            end
        end
    end
end

for n=1:N

    I = length(choiceData(n).group);  

    planIndex = zeros(I,1);

    for i = 1:I

        planIndex(i) = find(choiceData(n).planX(i,:)==1);

    end

    group = choiceData(n).group(1);

    year = choiceData(n).year(1)-2003;
    
    %*************************************************************************
    %single plan
    
    for  j = 1:4
        
        indexPlan = find(planIndex == j);
        
        if isempty(indexPlan) == 1,
            
            choiceData(n).choiceSim(1,6+j) = 99; %#ok<AGROW>
            
        else 
            
            choiceData(n).choiceSim(1,6+j) = indexPlan; %#ok<AGROW>  
            
        end
        
                choiceData(n).price(i,6+j) = 0; %#ok<AGROW>      
    end
            

    %*************************************************************************
    %optimal single plan by firm
    
    indexOne = find(planIndex == new(group,year).plan);

    choiceData(n).choiceSim(1,11) = indexOne; %#ok<AGROW>

    choiceData(n).price(i,11) = 0; %#ok<AGROW>

end

save structureWelfareData choiceData betaPhi betaXi;






