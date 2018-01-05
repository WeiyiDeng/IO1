% this script calculates allocations and welfare under observed alternative
% pricing scenarios

% efficient pricing is implemented using prices equal to estimated costs

% the code adds the following variables to the choiceData stucture
% welfare.gs(1,1:2) is v_ij - price_ij
% welfare.ss(1,1:2) is  v_ij - beta*price_ij + beta*cost_ij
% welfare.cost(1,1:2) is cost_ij
% welfare.price(1,1:2)
% welfare.dsim(1:,1:2)

clear all; clc;

load structureWelfareData;

N = length(choiceData);
SIM = 11;

for n = 1:N;

    I = length(choiceData(n).group);  

    planIndex = zeros(I,1);

    for i = 1:I

        planIndex(i) = find(choiceData(n).planX(i,:)==1);

    end

    offered = zeros(4,1);

    for i = 1:I
        offered(planIndex(i)) = 1;
    end
    
    %plans offered

    welfare(n).offered = offered; %#ok<AGROW>

    welfare(n).fam = choiceData(n).fam(1); %#ok<AGROW>

    welfare(n).pro = max(max(choiceData(n).proX)); %#ok<AGROW>

    %plan characteristics 

    phi = [choiceData(n).cost' choiceData(n).co' choiceData(n).ded' choiceData(n).ph'];

    %plan dummies and interactions with pro, pro90, age, male and inc

    xi = [choiceData(n).planX choiceData(n).proX choiceData(n).pro90X...
        choiceData(n).ageX choiceData(n).maleX choiceData(n).incX];

    %gs sets contributions equal to zero

    gs = phi(:,2:4)*betaPhi(2:4)' + xi*betaXi' + choiceData(n).epsilon;

    %social surplus sets contributions equal to cost

    ss = phi*betaPhi' + xi*betaXi' + choiceData(n).epsilon;
    
    for sim = 1:SIM;

        choice = choiceData(n).choiceSim(1,sim);
        
        if choice == 99
            
            welfare(n).plan(sim) = 99; %#ok<AGROW>
            
        else

            welfare(n).plan(sim) = planIndex(choice);  %#ok<AGROW>

            welfare(n).gs(sim) = 100*gs(choice); %#ok<AGROW>

            welfare(n).ss(sim) = 100*ss(choice); %#ok<AGROW>

            welfare(n).cost(sim) = 100*choiceData(n).cost(choice); %#ok<AGROW>

            % assign price by planIndex

            welfare(n).price(:,sim) = zeros(4,1); %#ok<AGROW>

            for i=1:I

                welfare(n).price(planIndex(i),sim) = 100*choiceData(n).price(i,sim); %#ok<AGROW>

            end
        
        end

    end

end

save welfareStats welfare;



