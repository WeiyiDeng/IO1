function [mBid MBid] = jointBidMuHighRisk(theta)

global bidData eBid zBid; %#ok<NUSED>

% %debug
% clear all; clc;
% load structureBidData;
% load jointOut;

%define parameters
delta = theta(29:30);
alpha = theta(31:34);
beta = theta(35:38);

%define length of vectors
J = length(bidData);
L = 8;
K = 10;

%zero vectors
e = zeros(J,1);
z = zeros(J,L);
de = zeros(J,K);

%loop over bidIndex
for j=1:J

    I = length(bidData(j).group);

    %zero sum variables
    costSum = 0;
    z14Sum = zeros(1,4);
    z58Sum = zeros(1,4);

    %loop over individuals
    for i=1:I

        %sum individual costs
        costSum = costSum...
            +alpha(bidData(j).plan(i))...
            +beta(bidData(j).plan(i))*(bidData(j).proPredict(i)-1);

        %instruments
        z14Sum(1,bidData(j).plan(i)) = z14Sum(1,bidData(j).plan(i))+1;
        z58Sum(1,bidData(j).plan(i)) = z58Sum(1,bidData(j).plan(i))+(bidData(j).proPredict(i)-1);

    end

    %moment
    e(j,1) = bidData(j).bid(i)-delta(1,bidData(j).carrier(i))*costSum/I;

    %instruments    
    z14 = z14Sum(1,:)/I;
    z58 = z58Sum(1,:)/I;
  
    z(j,:) = [z14 z58];

    %derivative
    d12 = zeros(1,2);
    d12(1,bidData(j).carrier(i)) = -costSum/I;
    d36 = -delta(1,bidData(j).carrier(i))*z14Sum/I;
    d710 = -delta(1,bidData(j).carrier(i))*z58Sum/I;
    
    de(j,:) = [d12 d36 d710];

end

mBid = z'*e/J;

MBid = z'*de/J;

zBid = z;

eBid = e;






