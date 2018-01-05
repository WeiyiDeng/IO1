%--------------------------------------------------------------------------
% create_figure7.m
%--------------------------------------------------------------------------

function [] = create_figure7(pK,dK,par,shadowcost)

global N gnum appmonth;

%--------------------------------------------------------------------------
% calculate expected revenues at alternative prices
%--------------------------------------------------------------------------

GRAPHPOINTS = 121;
STEPSIZE = 100;
SERIES = 11;
TOTALSERIES = 9*SERIES; %9 = 8 grades + all grades

t=(appmonth<14);
MarginChart = zeros(GRAPHPOINTS,TOTALSERIES);
for i=1:GRAPHPOINTS,
    tic
    x = pK + (STEPSIZE/1000)*(i-(GRAPHPOINTS+1)/2)+1.300;
    [s downtemp mindind fractemp deftind pvpmts pvrec crev ccost cprof urev uprof pr_less pr_sale pr_nosl pr_deft] = revenue_expectation(x,dK,par,shadowcost);
    %placeholder = zeros(N,1);
    MarginChart(i,1:SERIES) = [mean(s(t)) mean(mindind(s & t)) mean(deftind(s & t)) mean(downtemp(s & t)) mean(fractemp(s & t)) mean(pvpmts(s & t)) mean(crev(s & t)) mean(cprof(s & t)) mean(urev(t)) mean(uprof(t)) mean(pr_nosl(t))];
    for g=1:8,
        start = 1 + SERIES*g;
        finish = start + (SERIES-1);
        x1 = s(gnum==g & t);
        x2 = mindind(gnum==g & t);
        x3 = deftind(gnum==g & t);
        x4 = downtemp(gnum==g & t);
        x5 = fractemp(gnum==g & t); 
        x6 = pvpmts(gnum==g & t);
        x7 = crev(gnum==g & t);
        x8 = cprof(gnum==g & t);
        x9 = urev(gnum==g & t);
        x10 = uprof(gnum==g & t);
        x11 = pr_nosl(gnum==g & t);
        MarginChart(i,start:finish) = [mean(x1) mean(x2(x1)) mean(x3(x1)) mean(x4(x1)) mean(x5(x1)) mean(x6(x1)) mean(x7(x1)) mean(x8(x1)) mean(x9) mean(x10) mean(x11)];
    end;
    toc
end;

%--------------------------------------------------------------------------
% print output
%--------------------------------------------------------------------------

fid = fopen('figure6.txt','wt');
for i=1:GRAPHPOINTS,
    fprintf(fid,'%6.4f  ',MarginChart(i,:));
    fprintf(fid,'\n');
end
fclose(fid);

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------
