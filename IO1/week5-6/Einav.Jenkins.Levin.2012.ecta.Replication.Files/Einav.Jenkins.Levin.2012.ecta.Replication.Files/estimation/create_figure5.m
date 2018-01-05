%--------------------------------------------------------------------------
% create_figure5.m
%--------------------------------------------------------------------------

function [] = create_figure5(pKtemp,dKtemp,par)

global N SLTXRT docfeeK slsfeeK tlsfeeK gnum I less tK fracpaid censored default;

%--------------------------------------------------------------------------
%compute data for model fit assessment 
%--------------------------------------------------------------------------
[s down min fractemp deftemp pvpmt pvrec conrev concost conprof uncrev uncprof pr_less pr_sale pr_nosl pr_deft] = revenue_expectation(pKtemp,dKtemp,par,0);

%--------------------------------------------------------------------------
% Figure 5(a): Extra Down Payment Histogram
%--------------------------------------------------------------------------
actextra = tK - dKtemp;
modextra = down - dKtemp;
EDGES = [-Inf 0.1:0.1:2.0 Inf];
acthist = histc(actextra(I),EDGES)/sum(I);
modhist = histc(modextra(s),EDGES)/sum(s);
downout = [acthist(1:end-1); 1-mean(I); modhist(1:end-1); 1-mean(s)];
fid = fopen('figure4a.txt','wt');
for i=1:size(downout,1),
    fprintf(fid,'%6.4f  ',downout(i,:));
    fprintf(fid,'\n');
end
fclose(fid);

%--------------------------------------------------------------------------
% Figure 5(b): Fraction of Payments Histogram
%--------------------------------------------------------------------------
actfrac = fracpaid(~censored & I);
actdeft = fracpaid<1;
modfrac = fractemp(~censored & s);
EDGES = 0.00:0.05:1.0;
acthist = histc(actfrac,EDGES)/sum(~censored & I);
modhist = histc(modfrac,EDGES)/sum(~censored & s);
fracout = [acthist; modhist];
fid = fopen('figure4b.txt','wt');
for i=1:size(fracout,1),
    fprintf(fid,'%6.4f  ',fracout(i,:));
    fprintf(fid,'\n');
end
fclose(fid);

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------
