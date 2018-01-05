function [mCost MCost] = costCost(theta)

global costData eCost zCost; %#ok<NUSED>

%debug
% clear all; clc;
% load structureCostData;
% load costOut_2;
% theta = thetaHat;

%define parameters
alpha = theta(3:6);
beta = theta(7:10);
gamma = theta(11:14);

%mean coinsurance
meanCo = [85.53 83.96 95.14 77.7];

%define length of vectors
J = length(costData);
L = 12;
K = 14;

%zero vectors
e = zeros(J,1);
z = zeros(J,L);
de = zeros(J,K);

%loop over costIndex
for j=1:J

    I = length(costData(j).group);

    %zero sum variables
    costSum = 0;
    pmpmSum = 0;
    z14Sum = zeros(1,4);
    z58Sum = zeros(1,4);
    z912Sum = zeros(1,4);

    %loop over individuals
    for i=1:I

        %sum individual costs
        costSum = costSum...
            +alpha(costData(j).plan(i))...
            +beta(costData(j).plan(i))*(costData(j).pro(i)-1)...
            +gamma(costData(j).plan(i))*(costData(j).co(i)-meanCo(costData(j).plan(i)));

        %instruments
        z14Sum(1,costData(j).plan(i)) = z14Sum(1,costData(j).plan(i))+1;
        z58Sum(1,costData(j).plan(i)) = z58Sum(1,costData(j).plan(i))+costData(j).pro(i)-1;
        z912Sum(1,costData(j).plan(i)) = z912Sum(1,costData(j).plan(i))+(costData(j).co(i)-meanCo(costData(j).plan(i)));

        %sum pmpm
        pmpmSum = pmpmSum + costData(j).cost(i);

    end
    
    %moment
    e(j,1) = pmpmSum/I-costSum/I;

    %instruments
    z14 = z14Sum(1,:)/I;
    z58 = z58Sum(1,:)/I;
    z912 = z912Sum(1,:)/I;

    z(j,:) = [z14 z58 z912];

    %derivative
    d36 = -z14Sum/I;
    d710 = -z58Sum/I;
    d1114 = -z912Sum/I;

    de(j,:) = [0 0 d36 d710 d1114];

end

mCost = z'*e/J;

MCost = z'*de/J;

zCost = z;

eCost = e;









