%--------------------------------------------------------------------------
% optimal_mindown.m
%--------------------------------------------------------------------------
%
% Compute optimal minimum down payment requirements given estimated
% demand parameters
%
%--------------------------------------------------------------------------
function [optimal_down_matrix] = optimal_mindown(pK,dK,par,shadowcost,tax,type)

global N gnum;

%--------------------------------------------------------------------------
% calculate expected revenues at alternative prices
%--------------------------------------------------------------------------

optimal_prof_abs = -99*ones(9,3);
optimal_prof_rel = -99*ones(9,3);
optimal_down_abs = zeros(9,3);
optimal_down_rel = zeros(9,3);

GRAPHPOINTS = 61;
STEPSIZE = 50;

for i=1:GRAPHPOINTS,
tic
    % minimum down payments are set uniformly
    if (type==1),
    x = (STEPSIZE/1000)*(i-1);
    [s downtemp mindind fractemp deftind pvpmts pvrecovery crev ccost cprof urev uprof] = revenue_expectation(pK,x*ones(N,1),par,shadowcost);
    if (mean(uprof)>optimal_prof_abs(1,1)),
        optimal_prof_abs(1,1) = mean(uprof);
        optimal_down_abs(1,1) = x;
    end;
    for g=1:8,
        if (mean(uprof(gnum==g))>optimal_prof_abs(g+1,1)),
            optimal_prof_abs(g+1,1) = mean(uprof(gnum==g));
            optimal_down_abs(g+1,1) = x;
        end;   
    end;
    for t=0:1,
        if (mean(uprof(tax==t))>optimal_prof_abs(1,t+2)),
            optimal_prof_abs(1,t+2) = mean(uprof(tax==t));
            optimal_down_abs(1,t+2) = x;
        end;
        for g=1:8,
            if (mean(uprof(gnum==g & tax==t))>optimal_prof_abs(g+1,t+2)),
                optimal_prof_abs(g+1,t+2) = mean(uprof(gnum==g & tax==t));
                optimal_down_abs(g+1,t+2) = x;
            end;   
        end;
    end;
    end; % end if statement
    
    % minimum down payments are set in relation to observed minimum downs
    if (type==2),
    x = (STEPSIZE/1000)*(i-(GRAPHPOINTS+1)/2);
    [s downtemp mindind fractemp deftind pvpmts pvrecovery crev ccost cprof urev uprof] = revenue_expectation(pK,dK+x,par,shadowcost);
    if (mean(uprof)>optimal_prof_rel(1,1)),
        optimal_prof_rel(1,1) = mean(uprof);
        optimal_down_rel(1,1) = x;
    end;
    for g=1:8,
        if (mean(uprof(gnum==g))>optimal_prof_rel(g+1,1)),
            optimal_prof_rel(g+1,1) = mean(uprof(gnum==g));
            optimal_down_rel(g+1,1) = x;
        end;   
    end;
    for t=0:1,
        if (mean(uprof(tax==t))>optimal_prof_rel(1,t+2)),
            optimal_prof_rel(1,t+2) = mean(uprof(tax==t));
            optimal_down_rel(1,t+2) = x;
        end;
        for g=1:8,
            if (mean(uprof(gnum==g & tax==t))>optimal_prof_rel(g+1,t+2)),
                optimal_prof_rel(g+1,t+2) = mean(uprof(gnum==g & tax==t));
                optimal_down_rel(g+1,t+2) = x;
            end;   
        end;
    end;
    end; % end if statement
toc
end; % gridpoints loop

if (type==1),
    optimal_down_matrix = optimal_down_abs;
end;
if (type==2),
    optimal_down_matrix = optimal_down_rel;
end;

%--------------------------------------------------------------------------
% print output
%--------------------------------------------------------------------------
 
fid = fopen('optimal_down.txt','a');
for i=1:9,
    fprintf(fid,'%6.4f  ',optimal_down_abs(i,:));
    fprintf(fid,'\n');
end
fclose(fid);
%{
fid = fopen('Output/output_opt_down_rel.txt','wt');
for i=1:9,
    fprintf(fid,'%6.4f  ',optimal_down_rel(i,:));
    fprintf(fid,'\n');
end
fclose(fid);
%}
%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------
