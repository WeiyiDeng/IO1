function [mCost MCost] = jointCostMuHighRisk(theta)

global NS costData eCost zCost; %#ok<NUSED>

% %debug
% clear all; clc;
% load structureCostData;
% 
% theta = ones(1,39);

%define parameters
alpha = theta(31:34);
beta = theta(35:38);
rho = theta(39);

%define length of vectors
J = length(costData);
L = 8;
K = 11;

mCost = zeros(L,NS);
MCost = zeros(L,K,NS);
eCost = zeros(J,NS);
zCost = zeros(J,L,NS);

for ns = 1:NS

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
        d9Sum = 0;
        
        %loop over individuals
        for i=1:I

            proMu = costData(j).pro(i)+rho*costData(j).mu(i,ns)-1;

            %sum individual costs
            costSum = costSum...
                +alpha(costData(j).plan(i))...
                +beta(costData(j).plan(i))*proMu;

            %instruments
            z14Sum(1,costData(j).plan(i)) = z14Sum(1,costData(j).plan(i))+1;
            z58Sum(1,costData(j).plan(i)) = z58Sum(1,costData(j).plan(i))+proMu;
            
            %derivate
            d9Sum = d9Sum + beta(costData(j).plan(i))*costData(j).mu(i,ns);

            %sum pmpm
            pmpmSum = pmpmSum + costData(j).cost(i);

        end

        %moment
        e(j,1) = pmpmSum/I-costSum/I;

        z14 = z14Sum(1,:)/I;
        z58 = z58Sum(1,:)/I;

        z(j,:) = [z14 z58];

        %derivative
        d36 = -z14Sum/I;
        d710 = -z58Sum/I;
        d11 = -d9Sum/I;

        de(j,:) = [0 0 d36 d710 d11];

    end

    mCost(:,ns) = z'*e/J;

    MCost(:,:,ns) = z'*de/J;

    zCost(:,:,ns) = z;

    eCost(:,ns) = e;

end









