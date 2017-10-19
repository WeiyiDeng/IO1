function [cbapmmom, cbapmgrad] = cbapm(theta, datas, instrum);

% [cbapmmom, cbapmgrad] = cbapm(theta, datas, instrum);
% This function calculates:
% The sample moments (cbapmmom) and their gradient (cbapmgrad) in the 
% case of a CBAPM model with a CRRA utility function; see Hansen and Singleton (1982)
%
% INPUT
% theta: The parametrs of interest ion the CBAPM (dimension 2x1)
% datas: The dataset used for estimation (dimension Txp)
% instrum: The instruments used in estimation (dimension Txq)
%
% OUTPOUT
% cbapmmom: A Txq vector of the sample moments
% cbapmgrad: A qxp matrix, calculating the gradient of the moment conditions 



theta1 = theta(1,1);
theta2 = theta(2,1);
x1 = datas(:,1);
x2 = datas(:,2);

[T, numinstr] = size(instrum);

mom = theta2*x2.*(x1.^(theta1 - 1)) - 1; 
for i=1:numinstr
    cbapmmom(:,i) = mom.*instrum(:,i);
end

if nargout>1
    logx1 = log(x1);
    dmom1 = theta2*x2.*logx1.*(x1.^(theta1 - 1));
    dmom2 = x2.*x1.^(theta1 - 1);
    
    for j=1:numinstr
        cbapmgrad(j,1) = sum(dmom1.*instrum(:,j))/T;
        cbapmgrad(j,2) = sum(dmom2.*instrum(:,j))/T;
    end
end