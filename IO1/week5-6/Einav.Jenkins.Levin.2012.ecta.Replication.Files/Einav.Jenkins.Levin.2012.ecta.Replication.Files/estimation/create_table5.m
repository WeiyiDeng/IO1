%--------------------------------------------------------------------------
% create_table5.m
%--------------------------------------------------------------------------
% 
% This program calculates a supply-side objective function that is
% maximized over the firm's shadow cost of selling the marginal vehicle.
%
%--------------------------------------------------------------------------
% variable declarations
%--------------------------------------------------------------------------

function [optimal_cost optimal_rate] = create_table5(dpar)

global N APRUSE TERM X5 X6 pK dK gnum appmonth;
global mdK dpperiod lpK mperiod lK;
global S; %should be removed when S incorporated into revenue_expectation

%--------------------------------------------------------------------------
% define variables
%--------------------------------------------------------------------------
% grade classes
gclass = ones(N,1);
gclass(gnum>3) = 2;
gclass(gnum>5) = 3;

% tax season dummy
tax = zeros(N,1);
tax(appmonth==2) = 1;

T = size(mdK,2);
conprof = zeros(N,T);
saleind = zeros(N,T);

%--------------------------------------------------------------------------
% define deviations for objective function
%--------------------------------------------------------------------------

%PLEV = 0.100; %in $1000s
DLEV = 0.100; %in $1000s
%PPCT = 0.010;
%DPCT = 0.100;

%use uniform dollar level changes for estimation
%phigh  = pK + PLEV*ones(N,1);
%plow   = pK - PLEV*ones(N,1);
dhigh  = dK + DLEV*ones(N,1);
dlow   = dK - DLEV*ones(N,1);
        
%use uniform percent changes for estimation
%phigh  = pK*(1 + PPCT);
%plow   = pK*(1 - PPCT);
%dhigh  = dK*(1 + DPCT);
%dlow   = dK*(1 - DPCT);
            
%--------------------------------------------------------------------------
% initialize objective function matrices
%--------------------------------------------------------------------------

MAXS = 1;
MAXC = 100;

objfn_mindp = zeros(MAXS,MAXC);
objfn_class = zeros(MAXS,MAXC,3);
objfn_grade = zeros(MAXS,MAXC,8);
objfn_tax   = zeros(MAXS,MAXC,2);
%objfn_price = zeros(MAXS,MAXC);
%objfn_total = zeros(MAXS,MAXC);
objfn_learn = zeros(MAXS,MAXC,4);
%objfn_marg  = zeros(MAXS,MAXC,4);

minobj = 1e99 * ones(18,1);
mini = ones(18,1);
minj = ones(18,1);

%--------------------------------------------------------------------------
% search for minimum objective function over grid of S and shadowcost
%--------------------------------------------------------------------------
SHADOWSTEP = 50;
for i=1:MAXS,
    S = 10*i; shadowcost = 0.0;
    [s1 down min frac def pvpmt pvrec crev1 ccost cprof1 urev uprof] = revenue_expectation(pK,dK,dpar,shadowcost);
    [s2 down min frac def pvpmt pvrec crev2 ccost cprof2 urev uprof] = revenue_expectation(pK,dhigh,dpar,shadowcost);
    [s3 down min frac def pvpmt pvrec crev3 ccost cprof3 urev uprof] = revenue_expectation(pK,dlow,dpar,shadowcost);
    
    
    for t=1:T, % loop over DP periods
        dKtemp = mdK(:,t);
        [s down min frac def pvpmt pvrec crev ccost cprof urev uprof] = revenue_expectation(pK,dKtemp,dpar,shadowcost);
        conprof(:,t) = cprof;
        saleind(:,t) = s;
    end;
    %{
    for t=1:M, % loop over margin periods
        pKtemp = pK + (lpK(:,t) - lK);
        [s down min frac def pvpmt pvrec crev ccost conprof urev uprof] = revenue_expectation(pKtemp,dK,dpar,shadowcost);
        conprofM(:,t) = cprof;
        saleindM(:,t) = s;
    end;
    %}
    for j=1:MAXC,
        shadowcost = (SHADOWSTEP/1000)*j;
        
        %calculate ex ante expected profits at each pricing deviation
        act = s1.*(cprof1 - shadowcost);
        dhi = s2.*(cprof2 - shadowcost);
        dlo = s3.*(cprof3 - shadowcost);
        %phi = s4.*(cprof4 - shadowcost);
        %plo = s5.*(cprof5 - shadowcost);
        
        %------------------------------------------------------------------
        % calculate value of objective function using various moments
        %------------------------------------------------------------------
        
        % all grades
        mact = mean(act)*1000;
        mdhi = mean(dhi)*1000;
        mdlo = mean(dlo)*1000;
        objfn_mindp(i,j) = (max(0,mdhi-mact))^2+(max(0,mdlo-mact))^2;

        %mphi = mean(phi)*1000;
        %mplo = mean(plo)*1000;
        %objfn_price(i,j) = (min(0,mact-mphi))^2+(min(0,mact-mplo))^2;
        %objfn_total(i,j) = (min(0,mact-mphi))^2+(min(0,mact-mplo))^2+(min(0,mact-mdhi))^2+(min(0,mact-mdlo))^2;

        % grade classes (A/B, C, and D)
        for g=1:3,
            mact = mean(act(gclass==g))*1000;
            mdhi = mean(dhi(gclass==g))*1000;
            mdlo = mean(dlo(gclass==g))*1000;
            objfn_class(i,j,g) = (max(0,mdhi-mact))^2+(max(0,mdlo-mact))^2;
        end;

        % eight grades
        for g=1:8,
            mact = mean(act(gnum==g))*1000;
            mdhi = mean(dhi(gnum==g))*1000;
            mdlo = mean(dlo(gnum==g))*1000;
            objfn_grade(i,j,g) = (max(0,mdhi-mact))^2+(max(0,mdlo-mact))^2;
        end;

        % tax season
        for t=1:2,
            mact = mean(act(tax==t-1))*1000;
            mdhi = mean(dhi(tax==t-1))*1000;
            mdlo = mean(dlo(tax==t-1))*1000;
            objfn_tax(i,j,t) = (max(0,mdhi-mact))^2+(max(0,mdlo-mact))^2;
        end;

        %------------------------------------------------------------------
        % calculate objective based on minimum down improvement inequalities
        %------------------------------------------------------------------
        
        prev = saleind(:,1).*(conprof(:,1) - shadowcost);
        for t=2:T,
            appst = sum(dpperiod==t)/N;
            curr = saleind(:,t).*(conprof(:,t) - shadowcost);
            mcurr = mean(curr(dpperiod==t))*1000;
            mprev = mean(prev(dpperiod==t))*1000;
            mcur2 = mean(curr(dpperiod==t-1))*1000;
            mpre2 = mean(prev(dpperiod==t-1))*1000;
            objfn_learn(i,j,1) = objfn_learn(i,j,1) + (max(0,mprev-mcurr))^2;
            objfn_learn(i,j,2) = objfn_learn(i,j,2) + ((max(0,mprev-mcurr))^2)*appst;
            objfn_learn(i,j,3) = objfn_learn(i,j,3) + (max(0,mpre2-mcur2))^2;
            objfn_learn(i,j,4) = objfn_learn(i,j,4) + ((max(0,mpre2-mcur2))^2)*appst;
            prev = curr;
        end;
        %}
        %------------------------------------------------------------------
        % calculate objective based on margin improvement inequalities
        %------------------------------------------------------------------
        %{
        prev = saleindM(:,1).*(conprofM(:,1) - shadowcost);
        for t=2:M,
            appst = sum(mperiod==t)/N;
            curr = saleindM(:,t).*(conprofM(:,t) - shadowcost);
            mcurr = mean(curr(mperiod==t))*1000;
            mprev = mean(prev(mperiod==t))*1000;
            mcur2 = mean(curr(mperiod==t-1))*1000;
            mpre2 = mean(prev(mperiod==t-1))*1000;
            objfn_marg(i,j,1) = objfn_marg(i,j,1) + (min(0,mcurr-mprev))^2+(min(0,mcurr-mprev))^2;
            objfn_marg(i,j,2) = objfn_marg(i,j,2) + ((min(0,mcurr-mprev))^2+(min(0,mcurr-mprev))^2)*appst;
            objfn_marg(i,j,3) = objfn_marg(i,j,3) + (min(0,mcur2-mpre2))^2+(min(0,mcur2-mpre2))^2;
            objfn_marg(i,j,4) = objfn_marg(i,j,4) + ((min(0,mcur2-mpre2))^2+(min(0,mcur2-mpre2))^2)*appst;
            prev = curr;
       end;
       %}
        %------------------------------------------------------------------
        % store shadow costs that minimize respective objective functions
        %------------------------------------------------------------------
        
        curobj(1)  = objfn_mindp(i,j);
        curobj(2)  = objfn_class(i,j,1);
        curobj(3)  = objfn_class(i,j,2);
        curobj(4)  = objfn_class(i,j,3);
        curobj(5)  = objfn_grade(i,j,1);
        curobj(6)  = objfn_grade(i,j,2);
        curobj(7)  = objfn_grade(i,j,3);
        curobj(8)  = objfn_grade(i,j,4);
        curobj(9)  = objfn_grade(i,j,5);
        curobj(10) = objfn_grade(i,j,6);
        curobj(11) = objfn_grade(i,j,7);
        curobj(12) = objfn_grade(i,j,8);
        curobj(13) = objfn_tax(i,j,1);
        curobj(14) = objfn_tax(i,j,2);
        curobj(15) = objfn_learn(i,j,1);
        curobj(16) = objfn_learn(i,j,2);
        curobj(17) = objfn_learn(i,j,3);
        curobj(18) = objfn_learn(i,j,4);
        
        for k=1:18,
            if (curobj(k)<minobj(k))
                minobj(k)=curobj(k);
                mini(k)=i;
                minj(k)=j;
            end;
        end;

    end; %shadow cost loop
end; %discount rate loop

shadow_vector = SHADOWSTEP*minj;
fid = fopen('table4.txt','wt');
for i=1:size(shadow_vector,1),
    fprintf(fid,'%6.4f  ',shadow_vector(i,:));
    fprintf(fid,'\n');
end;
fclose(fid);

optimal_rate = 10.0;
optimal_cost = shadow_vector(1)/1000;

%--------------------------------------------------------------------------
% Output objective function values over grid
%--------------------------------------------------------------------------
%{
fid = fopen('supply_grid.txt','wt');
for i=1:MAXC,
    %fprintf(fid,'%6.4f  ',objfn_price(i,:));
    %fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_mindp(:,i));
    fprintf(fid,'\n');
    %{
    fprintf(fid,'%6.4f  ',objfn_total(i,:));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_class(i,:,1));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_class(i,:,2));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_class(i,:,3));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_grade(i,:,1));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_grade(i,:,2));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_grade(i,:,3));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_grade(i,:,4));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_grade(i,:,5));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_grade(i,:,6));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_grade(i,:,7));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_grade(i,:,8));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_tax(i,:,1));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_tax(i,:,2));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_learn(i,:,1));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_learn(i,:,2));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_learn(i,:,3));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_learn(i,:,4));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_marg(i,:,1));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_marg(i,:,2));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_marg(i,:,3));
    fprintf(fid,'\n');
    fprintf(fid,'%6.4f  ',objfn_marg(i,:,4));
    fprintf(fid,'\n');
    %}
end;
fclose(fid);
%}
%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------
