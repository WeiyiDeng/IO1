%--------------------------------------------------------------------------
% create_table6.m
%--------------------------------------------------------------------------

function [] = create_table6(pKtemp,dKtemp,par,shadowcost)

global N gnum appmonth;

P = 2; %revised for dpperiod pricing
CloseRates = zeros(P,9);
ExpProfits = zeros(P,9);
PreProfits = zeros(P,9);

%--------------------------------------------------------------------------
% Compute optimal downs
%--------------------------------------------------------------------------
t = appmonth==2;
OptimalDowns = optimal_mindown(pKtemp,dKtemp,par,shadowcost,t,1);

%--------------------------------------------------------------------------
% Pricing under observed prices
%--------------------------------------------------------------------------

[s downtemp mindind fractemp deftind pvpmts pvrecovery crev ccost cprof urev uprof] = revenue_expectation(pKtemp,dKtemp,par,shadowcost);
%onmargin = s & abovemar<=0.1;

CloseRates(10,1) = mean(s);
ExpProfits(10,1) = mean(cprof(s));
PreProfits(10,1) = mean(uprof);
for g=1:8,
    CloseRates(10,g+1) = mean(s(gnum==g));
    ExpProfits(10,g+1) = mean(cprof(s & gnum==g));
    PreProfits(10,g+1) = mean(uprof(gnum==g));
end;

%--------------------------------------------------------------------------
% Pricing under model optimal prices by credit grade
%--------------------------------------------------------------------------
for i=1:8,
    dKtemp(gnum==i) = OptimalDowns(i+1,1); %original
    %dKtemp(gnum==i) = OptimalDowns(i+1,3); %by dpperiod
end;

[s downtemp mindind fractemp deftind pvpmts pvrecovery crev ccost cprof urev uprof] = revenue_expectation(pKtemp,dKtemp,par,shadowcost);
%onmargin = s & abovemar<=0.1;

CloseRates(1,1) = mean(s);
ExpProfits(1,1) = mean(cprof(s));
PreProfits(1,1) = mean(uprof);
for g=1:8,
    CloseRates(1,g+1) = mean(s(gnum==g));
    ExpProfits(1,g+1) = mean(cprof(s & gnum==g));
    PreProfits(1,g+1) = mean(uprof(gnum==g));
end;
%{
% tax season only - optimal pricing by grade
CloseRates(6,1) = mean(s(t));
ExpProfits(6,1) = mean(cprof(s & t));
PreProfits(6,1) = mean(uprof(t));
for g=1:8,
    CloseRates(6,g+1) = mean(s(gnum==g & t));
    ExpProfits(6,g+1) = mean(cprof(s & gnum==g & t));
    PreProfits(6,g+1) = mean(uprof(gnum==g & t));
end;
%}
%--------------------------------------------------------------------------
% Pricing under model optimal price for all applicants
%--------------------------------------------------------------------------

dKtemp = OptimalDowns(1,1)*ones(N,1); %original
%dKtemp = OptimalDowns(1,3)*ones(N,1); %by dpperiod

[s downtemp mindind fractemp deftind pvpmts pvrecovery crev ccost cprof urev uprof] = revenue_expectation(pKtemp,dKtemp,par,shadowcost);
%onmargin = s & abovemar<=0.1;

CloseRates(2,1) = mean(s);
ExpProfits(2,1) = mean(cprof(s));
PreProfits(2,1) = mean(uprof);
for g=1:8,
    CloseRates(2,g+1) = mean(s(gnum==g));
    ExpProfits(2,g+1) = mean(cprof(s & gnum==g));
    PreProfits(2,g+1) = mean(uprof(gnum==g));
end;
%{
% tax season only - optimal uniform pricing
CloseRates(7,1) = mean(s(t));
ExpProfits(7,1) = mean(cprof(s & t));
PreProfits(7,1) = mean(uprof(t));
for g=1:8,
    CloseRates(7,g+1) = mean(s(gnum==g & t));
    ExpProfits(7,g+1) = mean(cprof(s & gnum==g & t));
    PreProfits(7,g+1) = mean(uprof(gnum==g & t));
end;
%}
%--------------------------------------------------------------------------
% Pricing with perfect price discrimination
%--------------------------------------------------------------------------
%{
% define parameters
[Lambda Delta Beta Gamma] = define_parameters(par);
dcoef = -Delta(K4+2);

%calculate threshold down payment
dKstar = max([X2 pKtemp]*Beta + ep + 0.00001,0);
XDtemp = [X4 pKtemp zeros(N,1) zeros(N,1)]*Delta;
dKsale = max((XDtemp+xi)/dcoef - 0.0001,0);

%calculate expected profit
[s1 down1 mind1 frac1 deft1 pvpmts1 pvrec1 crev1 ccost1 cprof1 urev1 uprof1 pr_less pr_sale pr_nosl pr_deft fracsum1] = revenue_expectation(pKtemp,dKstar,par,shadowcost);
[s2 down2 mind2 frac2 deft2 pvpmts2 pvrec2 crev2 ccost2 cprof2 urev2 uprof2 pr_less pr_sale pr_nosl pr_deft fracsum2] = revenue_expectation(pKtemp,dKsale,par,shadowcost);
cprofmax = max(s1.*fracsum1,s2.*fracsum2);
s=cprofmax>0;

CloseRates(3,1) = mean(s);
ExpProfits(3,1) = mean(cprofmax(s));
PreProfits(3,1) = CloseRates(3,1)*ExpProfits(3,1);
for g=1:8,
    CloseRates(3,g+1) = mean(s(gnum==g));
    ExpProfits(3,g+1) = mean(cprofmax(s & gnum==g));
    PreProfits(3,g+1) = CloseRates(3,g+1)*ExpProfits(3,g+1);
end;
%}
%--------------------------------------------------------------------------
% Optimal tax season pricing, by grade (not used)
%--------------------------------------------------------------------------
%{
for i=1:8,
    dKtemp(gnum==i & ~t) = OptimalDowns(i+1,2);
    dKtemp(gnum==i &  t) = OptimalDowns(i+1,3);
end;

[s downtemp mindind fractemp deftind pvpmts pvrecovery crev ccost cprof urev uprof] = revenue_expectation(pKtemp,dKtemp,par,shadowcost);
%onmargin = s & abovemar<=0.1;

CloseRates(4,1) = mean(s);
ExpProfits(4,1) = mean(cprof(s));
PreProfits(4,1) = mean(uprof);
for g=1:8,
    CloseRates(4,g+1) = mean(s(gnum==g));
    ExpProfits(4,g+1) = mean(cprof(s & gnum==g));
    PreProfits(4,g+1) = mean(uprof(gnum==g));
end;

% tax season only - optimal tax sesons pricing by grade
CloseRates(8,1) = mean(s(t));
ExpProfits(8,1) = mean(cprof(s & t));
PreProfits(8,1) = mean(uprof(t));
for g=1:8,
    CloseRates(8,g+1) = mean(s(gnum==g & t));
    ExpProfits(8,g+1) = mean(cprof(s & gnum==g & t));
    PreProfits(8,g+1) = mean(uprof(gnum==g & t));
end;

%--------------------------------------------------------------------------
% Optimal tax season pricing, uniform (not used)
%--------------------------------------------------------------------------

dKtemp(~t) = OptimalDowns(1,2);
dKtemp( t) = OptimalDowns(1,3);

[s downtemp mindind fractemp deftind pvpmts pvrecovery crev ccost cprof urev uprof] = revenue_expectation(pKtemp,dKtemp,par,shadowcost);
%onmargin = s & abovemar<=0.1;

CloseRates(5,1) = mean(s);
ExpProfits(5,1) = mean(cprof(s));
PreProfits(5,1) = mean(uprof);
for g=1:8,
    CloseRates(5,g+1) = mean(s(gnum==g));
    ExpProfits(5,g+1) = mean(cprof(s & gnum==g));
    PreProfits(5,g+1) = mean(uprof(gnum==g));
end;

% non tax season only - optimal uniform pricing
CloseRates(9,1) = mean(s(t));
ExpProfits(9,1) = mean(cprof(s & t));
PreProfits(9,1) = mean(uprof(t));
for g=1:8,
    CloseRates(9,g+1) = mean(s(gnum==g & t));
    ExpProfits(9,g+1) = mean(cprof(s & gnum==g & t));
    PreProfits(9,g+1) = mean(uprof(gnum==g & t));
end;
%}
%--------------------------------------------------------------------------
% output data
%--------------------------------------------------------------------------
panel1 = [CloseRates(1,:); ExpProfits(1,:)];
panel2 = [CloseRates(2,:); ExpProfits(2,:)];
Table5Data = [panel1; panel2];

% print to file
fid = fopen('table5.txt','a');
for i=1:2*P,
    fprintf(fid,'%6.4f  ',Table5Data(i,:));
    fprintf(fid,'\n');
end;
fclose(fid);
end;

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------
