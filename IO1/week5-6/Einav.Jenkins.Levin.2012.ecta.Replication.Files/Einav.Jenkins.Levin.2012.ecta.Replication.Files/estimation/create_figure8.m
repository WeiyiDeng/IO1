%--------------------------------------------------------------------------
% create_figure8.m
%--------------------------------------------------------------------------

function [] = create_figure8(pKtemp,dKtemp,par,shadowcost)

global N X1 X2 X3 X4 gnum;

%--------------------------------------------------------------------------
% Panels (a) and (b): No down payment
%--------------------------------------------------------------------------

pKmean = pKtemp;
x = 0*ones(N,1); %set lower limit for minimum down payment
[s downtemp mindind fractemp deftind] = revenue_expectation(pKmean,x,par,shadowcost);

t=1;
SERIES = 3;
FigureData1 = zeros(9,SERIES);
FigureData1(1,1:SERIES) = [mean(s) mean(downtemp(s)) mean(deftind(s))];
for g=1:8,
        x1 = s(gnum==g & t);
        %x2 = mindind(gnum==g & t);
        x3 = deftind(gnum==g & t);
        x4 = downtemp(gnum==g & t);
        %x5 = fractemp(gnum==g & t); 
        %x6 = pvpmts(gnum==g & t);
        %x7 = crev(gnum==g & t);
        %x8 = cprof(gnum==g & t);
        %x9 = urev(gnum==g & t);
        %x10 = uprof(gnum==g & t);
        %x11 = placeholder(gnum==g & t);
        FigureData1(g+1,1:SERIES) = [mean(x1) mean(x4(x1)) mean(x3(x1))];
end;
FigureData1

%--------------------------------------------------------------------------
% Panels (b) and (c): Average down payments
%--------------------------------------------------------------------------

%roundedmin = [0.0; 0.4; 0.5; 0.6; 0.8; 1.0; 1.2; 1.4];
gpct = zeros(8,1);
x = zeros(N,1);
for i=1:8,
    gpct(i) = mean(gnum==i);
    x(gnum==i) = mean(dKtemp(gnum==i));
    %x(gnum==i) = roundedmin(i);
end;

[s downtemp mindind fractemp deftind] = revenue_expectation(pKmean,x,par,shadowcost);

SERIES = 3;
FigureData2 = zeros(9,SERIES);
FigureData2(1,1:SERIES) = [mean(s) mean(downtemp(s)) mean(deftind(s))];
for g=1:8,
        x1 = s(gnum==g & t);
        %x2 = mindind(gnum==g & t);
        x3 = deftind(gnum==g & t);
        x4 = downtemp(gnum==g & t);
        %x5 = fractemp(gnum==g & t); 
        %x6 = pvpmts(gnum==g & t);
        %x7 = crev(gnum==g & t);
        %x8 = cprof(gnum==g & t);
        %x9 = urev(gnum==g & t);
        %x10 = uprof(gnum==g & t);
        %x11 = placeholder(gnum==g & t);
        FigureData2(g+1,1:SERIES) = [mean(x1) mean(x4(x1)) mean(x3(x1))];
end;
FigureData2

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------