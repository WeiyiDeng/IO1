function [mBid MBid] = costBid(theta)

global bidData eBid zBid; %#ok<NUSED>

% %debug
% clear all; clc;
% load structureBidData;
% load jointOut;

%define parameters
delta = theta(1:2);
alpha = theta(3:6);
beta = theta(7:10);
gamma = theta(11:14);

%mean coinsurance
meanCo = [85.53 83.96 95.14 77.7];

%define length of vectors
J = length(bidData);
L = 12;
K = 14;

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
    z912Sum = zeros(1,4);

    %loop over individuals 
    for i=1:I
        
        %sum individual costs
        costSum = costSum...
            +alpha(bidData(j).plan(i))...
            +beta(bidData(j).plan(i))*(bidData(j).proPredict(i)-1)...
            +gamma(bidData(j).plan(i))*(bidData(j).co(i)-meanCo(bidData(j).plan(i)));

        %instruments
        z14Sum(1,bidData(j).plan(i)) = z14Sum(1,bidData(j).plan(i))+1;
        z58Sum(1,bidData(j).plan(i)) = z58Sum(1,bidData(j).plan(i))+(bidData(j).proPredict(i)-1);
        z912Sum(1,bidData(j).plan(i)) = z912Sum(1,bidData(j).plan(i))+(bidData(j).co(i)-meanCo(bidData(j).plan(i)));
        
    end

    %moment
    e(j,1) = bidData(j).bid(i)-delta(1,bidData(j).carrier(i))*costSum/I;
    
    %instruments    
    z14 = z14Sum(1,:)/I;
    z58 = z58Sum(1,:)/I;
    z912 = z912Sum(1,:)/I;
  
    z(j,:) = [z14 z58 z912];

    %derivative
    d12 = zeros(1,2);
    d12(1,bidData(j).carrier(i)) = -costSum/I;
    d36 = -delta(1,bidData(j).carrier(i))*z14Sum/I;
    d710 = -delta(1,bidData(j).carrier(i))*z58Sum/I;
    d1114 = -delta(1,bidData(j).carrier(i))*z912Sum/I;
    
    de(j,:) = [d12 d36 d710 d1114];
      
end

mBid = z'*e/J;

MBid = z'*de/J;

zBid = z;

eBid = e;






