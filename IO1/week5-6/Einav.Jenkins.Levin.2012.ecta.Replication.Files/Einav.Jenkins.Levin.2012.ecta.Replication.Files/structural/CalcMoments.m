%--------------------------------------------------------------------------
% CalcMoments.m
%--------------------------------------------------------------------------

function [MOMENTS SIMMOMENTS OBJFN] = CalcMoments(s,simless,simxtra,simdef,simfrac,s2,s3)

global dscale pscale;

MOMENTS = [0.6550; 0.4273; 0.5771; 0.6082; 0.3435; 0.02244; 0.00152];

m1=mean(~s);                        %probability of non purchase
m2=mean(simless(s));                %probability of minimum down
m3=mean(simxtra(s & ~simless));     %mean of extra down if extra
m4=mean(simdef);                    %probability of default if sale
m5=mean(simfrac(simdef));           %mean of payments if default
m6=(mean(s)-mean(s2))/dscale;       %change in close w.r.t. min down
m7=(mean(s)-mean(s3))/pscale;       %change in close w.r.t. price

%compute objective function
SIMMOMENTS = [m1; m2; m3; m4; m5; m6; m7];
W = 10^5*eye(size(SIMMOMENTS,1)); 
W(6,6)=W(6,6)*100; 
W(7,7)=W(7,7)*10000;
OBJFN = (MOMENTS-SIMMOMENTS)'*W*(MOMENTS-SIMMOMENTS);

%--------------------------------------------------------------------------
% end of function
%--------------------------------------------------------------------------
