function f = GMMObjFun_Start(x0)

global XX y IV W tol_inner Kbeta Ktheta ns v share TM prods Total T g

theta1 = x0(1:Kbeta,1);                       % Mean of tastes
theta2 = x0(Kbeta+1:Kbeta+Ktheta,1);          % Std deviation of tastes

theta2_mat = diag(theta2);                    % Ktheta*Ktheta           4*4

sigma_jm = exp(y);                            % initial guess of sigma using DV of logit model

ii          = 0;
norm_max    = 1;
while norm_max > tol_inner && ii < 10000      % Loop until convergence
    s_jm_prep_mat = exp(XX*theta2_mat*v+repmat(sigma_jm,1,ns));         % obtain the numerator for computing the mkt share integral sjm    % Total*ns  970*1000      
    sum_s_jm_expand = zeros(Total,ns);                                  % obtain the denominator for computing the mkt share integral sjm
    for j=1:TM
        sum_s_jm = sum(s_jm_prep_mat(T(j,1):T(j,2),:),1)+1;
        temp = repmat(sum_s_jm,prods(j),1);
        sum_s_jm_expand(T(j,1):T(j,2),:) = temp;
    end
    sjm = mean(s_jm_prep_mat./sum_s_jm_expand, 2);          % compute mkt share integral sjm         % Total*1     970*1
    
    sigma_jm_previous = sigma_jm;
    sigma_jm = sigma_jm_previous+log(share)-log(sjm);       % contraction mapping
    norm_max = max(abs(sigma_jm-sigma_jm_previous));

    ii           = ii + 1;
end

xi_jm = sigma_jm - [ones(Total,1) XX]*theta1;               % Total*1     970*1

g = IV'*xi_jm;                                              % nIV*1       13*1

f       = g'*W*g;                                           % Objective function

end