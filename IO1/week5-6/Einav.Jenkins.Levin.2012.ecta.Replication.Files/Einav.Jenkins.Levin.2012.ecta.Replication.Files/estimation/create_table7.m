%--------------------------------------------------------------------------
% create_table7.m
%--------------------------------------------------------------------------

function [] = create_table7(pKtemp,dKtemp,par,shadowcost)

global N gnum;

%--------------------------------------------------------------------------
% Set parameters
%--------------------------------------------------------------------------

GRAPHPOINTS = 41;
STEPSIZE = 50;
IncProfUniform = zeros(GRAPHPOINTS,GRAPHPOINTS,9);
EntProfUniform = zeros(GRAPHPOINTS,GRAPHPOINTS,9);

%(optional) allow entrant to undercut incumbent by epsilon
EPSILON = 0.000;

%--------------------------------------------------------------------------
% Compute incumbent prices under model optimal price for all applicants
%--------------------------------------------------------------------------

% calculate profits at grid of pricing strategies
for i=1:GRAPHPOINTS,
    dKinc = ((STEPSIZE/1000)*(i-1))*ones(N,1);
    for j=1:GRAPHPOINTS,
        dKent = ((STEPSIZE/1000)*(j-1))*ones(N,1);
        [cr cp uri upi ure upe non inc ent div] = revenue_competition(pKtemp,dKinc,dKent,par,shadowcost);
        IncProfUniform(i,j,1) = mean(upi);
        EntProfUniform(i,j,1) = mean(upe);
        for g=1:8,
            IncProfUniform(i,j,g+1) = mean(upi(gnum==g));
            EntProfUniform(i,j,g+1) = mean(upe(gnum==g));
        end;    
    end;
end;

%--------------------------------------------------------------------------
% find equilibrium
%--------------------------------------------------------------------------
apps          = zeros(8,1);
IncBestMindwn = zeros(GRAPHPOINTS,9);
IncBestProfit = zeros(GRAPHPOINTS,9);
EntProfits    = zeros(GRAPHPOINTS,10);
EntBestProfit = zeros(9,1);
EntBestMindwn = zeros(9,1);
BestTotal     = 0;

for j=1:GRAPHPOINTS,
    
    % for each price of the entrant find the optimal response of incumbent
    IncBestProfit(j,1) = max(IncProfUniform(:,j,1));
    for i=1:GRAPHPOINTS,
       if (IncBestProfit(j,1)==IncProfUniform(i,j,1)),
           IncBestMindwn(j,1) = i;
        end;
    end;

    % for this price, the entrant gets profits based on inc optimal choice
    EntProfits(j,1) = EntProfUniform(IncBestMindwn(j,1),j,1);
    if (EntProfits(j,1)  > EntBestProfit(1)),
        EntBestProfit(1) = EntProfits(j,1);
        EntBestMindwn(1)  = j;
    end;

    % repeat for each grade
    for g=1:8,
        IncBestProfit(j,g+1) = max(IncProfUniform(:,j,g+1));
        for i=1:GRAPHPOINTS,
           if (IncBestProfit(j,g+1)==IncProfUniform(i,j,g+1)),
               IncBestMindwn(j,g+1) = i;
            end;
        end;
        EntProfits(j,g+1) = EntProfUniform(IncBestMindwn(j,g+1),j,g+1);
        if (EntProfits(j,g+1) > EntBestProfit(g+1)),
            EntBestProfit(g+1) = EntProfits(j,g+1);
            EntBestMindwn(g+1) = j;
        end;
        apps(g) = mean(gnum==g);
    end;

    % calculate weighted (by apps) total profit for this entrant price
    EntProfits(j,10) = apps'*EntProfits(j,2:9)';
    if (EntProfits(j,10) > BestTotal),
        BestTotal = EntProfits(j,10);
        BestVsOpt = j;
    end;
end;
%{
IncBestMindwn
IncBestProfit
EntProfits
EntBestMindwn
EntBestProfit
BestVsOpt
%}
%--------------------------------------------------------------------------
% Output results
%--------------------------------------------------------------------------
% print to file
fid = fopen('table6I.txt','wt');
for i=1:GRAPHPOINTS,
    fprintf(fid,'%6.4f  ',IncProfUniform(i,:,1));
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('table6E.txt','wt');
for i=1:GRAPHPOINTS,
    fprintf(fid,'%6.4f  ',EntProfUniform(i,:,1));
    fprintf(fid,'\n');
end
fclose(fid);

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------
