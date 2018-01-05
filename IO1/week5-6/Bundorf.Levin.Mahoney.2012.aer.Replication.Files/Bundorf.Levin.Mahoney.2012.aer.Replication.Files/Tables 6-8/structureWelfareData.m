clear all; clc;

load structureChoiceData;

afterTaxPPP = -.012894;

alpha_premium = 100*afterTaxPPP; 

premium = -1;

coinsurance = 2.6;

deductible = 0.02; 

phnonstandard = -21.52;

betaPhi = [-1 coinsurance deductible*10 phnonstandard/100];

constant = [-111.43 -99.79 0 -86.32];

risk = [-5.59 -0.34 0 -4.25];

risk90 = zeros(1,4);

age = [1.69 0.16 0 1.39];

male = [14.67 -41.71 0 -34.28];

income = zeros(1,4);

betaXi = [constant/100 risk/100 risk90/100 age male/100 income];

alpha = [219.9085 246.2935 236.5411 247.9076];
 
beta = [171.6085 176.4882 73.8935 127.6227];

N = length(choiceData);

%seed random number generator
rand('twister', sum(100*clock))

for n = 1:N
    
    %I is size of choice set

    I = length(choiceData(n).group); 

    %generate type I extreme value random variables

    mu = rand(I,1);

    epsilon = -log(-log(mu));

    choiceData(n).epsilon = epsilon/alpha_premium; %#ok<AGROW>

    % cost_ij = alpha_j*choiceData.planX + beta_j*(choiceData.proX-1)

    cost = choiceData(n).planX*alpha'+(choiceData(n).proX-choiceData(n).planX)*beta';

    choiceData(n).cost = cost'/100; %#ok<AGROW>

end

save structureWelfareData choiceData betaPhi betaXi;
