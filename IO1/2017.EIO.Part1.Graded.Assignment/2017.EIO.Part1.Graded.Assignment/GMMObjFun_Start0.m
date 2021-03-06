function f = GMMObjFun_Start(x0)

global XX y IV W tol_inner Kbeta Ktheta ns v share TM prods Total T g exp_delta_store prep_shareSum

theta1 = x0(1:Kbeta,1);                       % Mean of tastes
theta2 = x0(Kbeta+1:Kbeta+Ktheta,1);          % Std deviation of tastes

theta2_mat = diag(theta2);                    % Ktheta*Ktheta           4*4

% delta_jm = exp(y);                          % initial guess of sigma using DV of logit model
% delta_jm = y;                               % initial guess of sigma using DV of logit model
exp_delta_jm = exp_delta_store;                       % store and use new value of delta_jm computed from last loop

ii          = 0;
norm_max    = 1;
exp_mu_ijm = exp(XX*theta2_mat*v);
while norm_max > tol_inner && ii < 10000      % Loop until convergence
    s_jm_prep_mat = exp_mu_ijm.*(repmat(exp_delta_jm,1,ns));         % obtain the numerator for computing the mkt share integral sjm    % Total*ns  970*1000      
%     s_jm_prep_mat = exp(XX*theta2_mat*v+repmat(delta_jm,1,ns));         % obtain the numerator for computing the mkt share integral sjm    % Total*ns  970*1000      
    % NP: The first part of your RHS expression can be calculated outside
    % (before) the loop, and then you can just have exp(mu)*exp(delta)
    
    shareSum = prep_shareSum*s_jm_prep_mat;                           % TM*ns    50*1000
    shareSum_expand = prep_shareSum'*shareSum+1;                      % obtain the denominator for computing the mkt share integral sjm        % Total*ns  970*1000 
        
%     for j=1:TM
%         sum_s_jm = sum(s_jm_prep_mat(T(j,1):T(j,2),:),1)+1;
%         temp = repmat(sum_s_jm,prods(j),1);
%         sum_s_jm_expand(T(j,1):T(j,2),:) = temp;
%     end
    % NP: This loop inside a loop will slow you down, you should vectorize
    % this calculation by pre-generating a (50*970) matrix with appropriate
    % ones to make the summation for the denominator faster (just 1 line of matrix multiplication)

    
%     sjm = mean(s_jm_prep_mat./sum_s_jm_expand, 2);          % compute mkt share integral sjm         % Total*1     970*1
    sjm = mean(s_jm_prep_mat./shareSum_expand, 2);          % compute mkt share integral sjm         % Total*1     970*1
    
%     delta_jm_previous = delta_jm;
%     delta_jm = delta_jm_previous+log(share)-log(sjm);       % contraction mapping
    exp_delta_jm_previous = exp_delta_jm;
    exp_delta_jm = exp_delta_jm_previous.*share./sjm;       % contraction mapping
    % NP: It should be faster to use exponentials instead of logs in this
    % formula: exp(delta_jm') = exp(delta_jm).*share./sjm
    
%     norm_max = max(abs(delta_jm-delta_jm_previous));
    norm_max = max(abs(log(exp_delta_jm./exp_delta_jm_previous)));

    ii           = ii + 1;
end

xi_jm = log(exp_delta_jm) - [ones(Total,1) XX]*theta1;               % Total*1     970*1

g = IV'*xi_jm;                                              % nIV*1       13*1

f       = g'*W*g;                                           % Objective function

exp_delta_store = exp_delta_jm;

end