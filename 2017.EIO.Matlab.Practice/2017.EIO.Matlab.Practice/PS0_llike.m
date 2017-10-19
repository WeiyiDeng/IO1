function [llike,grad,hess] = PS0_llike(beta)

global X Y

param = size(X,2);
XB    = X*beta;
Pr    = normcdf(XB);
pdf   = normpdf(XB);
llike = -sum(Y.*log(Pr)+(1-Y).*log(1-Pr));
A     = (Y./Pr).*pdf.*(-pdf./Pr-XB) + ((1-Y)./(1-Pr)).*(-pdf).*(pdf./(1-Pr)-XB);
B     = (Y./Pr).*pdf + ((1-Y)./(1-Pr)).*(-pdf);

% w: provide analytical hessian and gradients computed by hand
% Gradient:                                    % w: derivative check?
grad = zeros(param,1);
for k=1:param
    grad(k,1) = -sum(B.*X(:,k));               % w: how to compute gradient and hessian?
end


% Hessian:
hess  = zeros(param);
for j=1:param
    for k=1:param
        hess(j,k)  = -sum(A.*X(:,j).*X(:,k));
    end
end


%{
% For Logit model:
Pr    = exp(X*beta)./(1+exp(X*beta));
llike = -sum(Y.*log(Pr)+(1-Y).*log(1-Pr));
grad  = -sum(repmat(Y-Pr,1,size(X,2)).*X);
hess  = zeros(size(X,2));
for i=1:size(X,2)
    for j=1:size(X,2)
        hess(i,j) = sum(X(:,i).*X(:,j).*Pr.*(1-Pr));
    end
end
%}

end