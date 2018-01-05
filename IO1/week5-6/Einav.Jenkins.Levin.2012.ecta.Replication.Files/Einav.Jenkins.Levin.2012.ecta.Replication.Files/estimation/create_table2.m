%--------------------------------------------------------------------------
% create_table2.m
%--------------------------------------------------------------------------

function [] = create_table2(pKtemp,dKtemp,par)

global N SLTXRT docfeeK slsfeeK tlsfeeK gnum I less tK fracpaid censored default;

%--------------------------------------------------------------------------
% Table 2: Model Fit - Actual Data
%--------------------------------------------------------------------------

ActData=zeros(6,9); ModData=zeros(6,9);

%create comparable measure of loan size
ptilda = pKtemp.*(1+SLTXRT) + docfeeK + slsfeeK + tlsfeeK;
actloan = ptilda - tK;
modloan = ptilda - down;

ActData(1,1) = mean(I);
ActData(2,1) = mean(less(I));
ActData(3,1) = mean(actloan(I));
ActData(4,1) = 1 - mean(actfrac<1);
ActData(5,1) = mean(fracpaid(I & ~censored));
ActData(6,1) = 0;
for i=1:8,
    ActData(1,i+1) = mean(I(gnum==i));
    ActData(2,i+1) = mean(less(I & gnum==i));
    ActData(3,i+1) = mean(actloan(I & gnum==i));
    ActData(4,i+1) = 1 - mean(actdeft(I & ~censored & gnum==i));
    ActData(5,i+1) = mean(fracpaid(I & ~censored & gnum==i));
    ActData(6,i+1) = 0;
end;

%--------------------------------------------------------------------------
% Table 2: Model Fit - Model Data
%--------------------------------------------------------------------------
[s down min fractemp deftemp pvpmt pvrec conrev concost conprof uncrev uncprof pr_less pr_sale pr_nosl pr_deft] = revenue_expectation(pKtemp,dKtemp,par,0);

ModData(1,1) = mean(s);
ModData(2,1) = mean(min(s));
ModData(3,1) = mean(modloan(s));
ModData(4,1) = 1 - mean(deftemp(s & ~censored));
ModData(5,1) = mean(fractemp(s & ~censored));
ModData(6,1) = 0;
for i=1:8,
    ModData(1,i+1) = mean(s(gnum==i));
    ModData(2,i+1) = mean(min(s & gnum==i));
    ModData(3,i+1) = mean(modloan(s & gnum==i));
    ModData(4,i+1) = 1 - mean(deftemp(s & ~censored & gnum==i));
    ModData(5,i+1) = mean(fractemp(s & ~censored & gnum==i));
    ModData(6,i+1) = 0;
end;

%--------------------------------------------------------------------------
% output data
%--------------------------------------------------------------------------

% print to file
fid = fopen('table2.txt','wt');
for i=1:6,
    fprintf(fid,'%6.4f  ',ActData(i,:));
    fprintf(fid,'\n');
end
for i=1:6,
    fprintf(fid,'%6.4f  ',ModData(i,:));
    fprintf(fid,'\n');
end
fclose(fid);

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------