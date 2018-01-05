function f = MergerObjFun(p_star)

global XX y IV W tol_inner Kbeta Ktheta ns v share TM prods Total T g exp_delta_store prep_shareSum...
    Kgamma price z xi_jm_store xi_jm_star A_star omega_matrices mc_m m theta2 theta1

XX_star = [A_star p_star];
exp_mu_ijm_star = exp(XX_star*diag(theta2)*v);     % prods*1000
exp_delta_jm_star = exp([ones(prods(m),1) XX_star]*theta1 + xi_jm_star);    % prods*1
s_jm_prep_mat_star = exp_mu_ijm_star.*(repmat(exp_delta_jm_star,1,ns));     % prods*1000
sum_s_jm_star = sum(s_jm_prep_mat_star,1)+1;
sum_s_jm_expand_star = repmat(sum_s_jm_star,prods(m),1);
sjm_star = mean(s_jm_prep_mat_star./sum_s_jm_expand_star, 2);    % prods*1

f = p_star-mc_m-inv(omega_matrices{m})*sjm_star;

end
