%--------------------------------------------------------------------------
% computes value and policy functions for all states (i.e., all y and l)
%--------------------------------------------------------------------------

function [P EVP EVD D] = optimal_repayment(sig2,delta,beta,psi,phi,theta)

global T ynum lnum vnum ystep ystates vstates up ud v;

% define state transition matrix
P = zeros(ynum,ynum);
points = [-Inf ystates(1:end-1)+ystep/2 Inf];
for s=1:ynum,
    temp = normcdf(points,ystates(s),sig2^0.5);
    P(s,:) = diff(temp);
end

% generate terminal values for each state (all YxLxV md arrays)
VTY = repmat((beta/(1-beta))*ystates',[1,lnum,vnum]);
VTV = ones(1,1,vnum); VTV(1,1,:) = vstates;
VTV = repmat(1/(1-beta*delta)*VTV,[ynum,lnum,1]);
VTP = VTY + VTV;
VTD = VTY + psi*VTV;

% generate matrices for choice and value function storage
D   = zeros(ynum,lnum,vnum,T);
EVP = zeros(ynum,lnum,vnum,T+1); EVP(:,:,:,T+1) = VTP;
EVD = zeros(ynum,lnum,vnum,T+1); EVD(:,:,:,T+1) = VTD;

% solve dynamic programming problem for all states and times
for t=T:-1:1,

    vp = theta*up +     delta^t*v + beta*EVP(:,:,:,t+1);
    vd = theta*ud + psi*delta^t*v + beta*EVD(:,:,:,t+1) - phi;
    D(:,:,:,t) = vp > vd;

    %calculate expected value from payment
    temp = D(:,:,:,t).*vp + (1-D(:,:,:,t)).*vd;
    temp = reshape(temp,ynum,lnum*vnum); temp = P*temp;
    EVP(:,:,:,t) = reshape(temp,ynum,lnum,vnum);

    %calculate expected value from default
    temp = vd;
    temp = reshape(temp,ynum,lnum*vnum); temp = P*temp;
    EVD(:,:,:,t) = reshape(temp,ynum,lnum,vnum);
    clear vp vd temp;
end;

%--------------------------------------------------------------------------
% end of function
%--------------------------------------------------------------------------
