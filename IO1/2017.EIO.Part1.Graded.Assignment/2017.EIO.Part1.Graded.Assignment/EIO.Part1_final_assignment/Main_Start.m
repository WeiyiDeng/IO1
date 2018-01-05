%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMPIRICAL INDUSTRIAL ORGANIZATION (230319) - PART 1 - FINAL ASSIGNMENT %
% by Weiyi Deng (revised version, much faster!)         Nov 20 2017
% w: var-cov problem not solved yet 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
global XX y IV W tol_inner Kbeta Ktheta ns v share TM prods Total T g exp_delta_store prep_shareSum

% This gets you started for the following models
% (1) Logit with and without IV
% (2) Nested Logit
% (3) Random Coefficients Logit


%%%%%%%%%%%%%%%%
%%% SETTINGS %%%
%%%%%%%%%%%%%%%%

DATA       = csvread('Data.csv');
IDmkt      = DATA(:,1);                 % Market identifier
IDprod     = DATA(:,2);                 % Product identifier
share      = DATA(:,3);                 % Market share
A          = DATA(:,4:6);               % Product characteristics
price      = DATA(:,7);                 % Price
z          = DATA(:,8:10);              % Instruments
% group      = DATA(:,11);                % Group identifier for Nested Logit
TM         = max(IDmkt);                % Number of markets
prods      = zeros(TM,1);               % # of products in each market
for m=1:TM
    prods(m,1) = max(IDprod(IDmkt==m,1));
end
T          = zeros(TM,2);               % w: T is the end row of products in each market
T(1,1)     = 1;
T(1,2)     = prods(1,1); 
for i=2:TM
    T(i,1) = T(i-1,2)+1;                % 1st Column market starting point
    T(i,2) = T(i,1)+prods(i,1)-1;       % 2nd Column market ending point
end
Total      = T(TM,2);                   % # of obsevations
TotalProd  = max(prods);                % Max # of products in a given market
% Ngroups    = max(group);                % # of groups

%% Q1 BLP ESTIMATION OF RANDOM COEFFICIENTS LOGIT MODEL
% w: use y (not exp(y)) (y is DV of logit) as starting guess of delta
% dim 970*1000

tol_inner = 1.e-14;                         % Tolerance for inner loop (NFXP)
tol_outer = 1.e-6;                          % Tolerance for outer loop (Estimation)
randn('seed',1);                           % Reset normal random number generator
rand('seed',1);                            % Reset uniform random number generator
Kbeta   = 2+size(A,2);                      % # of parameters in mean utility
Ktheta  = 1+size(A,2);                      % # of parameters with random coefficient (w: # of SE to be estimated)
% ns      = 1000;                             % # draws to simulate shares
ns      = 200;                             % # draws to simulate shares
v       = randn(Ktheta,ns);                 % Draws for share integrals during estimation
IV      = [ones(Total,1) A z A.^2 z.^2];    % Instruments
nIV     = size(IV,2);                       % # instrumental variables
W       = (IV'*IV)\eye(nIV);                % Starting GMM weighting matrix
x0      = rand(Kbeta+Ktheta,1);             % Starting values for all parameters
opts    = optimset('Display','iter','TolCon',1E-6,'TolFun',1E-10,'TolX',1E-10);
x_L     = [-Inf*ones(Kbeta,1);zeros(Ktheta,1)];     % Lower bounds is zero for standard deviations of random coefficients
x_U     = Inf.*ones(Kbeta+Ktheta,1);                % Upper bounds for standard deviations of random coefficients

XX = [A price];
exp_delta_store = exp(y);
prep_shareSum = zeros(TM,Total);
for r = 1:TM
    prep_shareSum(r,T(r,1):T(r,2))=1;
end
    
tic
[X,fval_rep,exitflag,output,lambda,grad,hessian] = fmincon(@GMMObjFun_Start0,x0,[],[],[],[],x_L,x_U,[],opts);
toc
theta1 = X(1:Kbeta,1);
theta2 = X(Kbeta+1:Kbeta+Ktheta,1);

standard_errors = sqrt(diag(inv(hessian)));         % standard errors incorrect? how to compute?
% ref https://www.mathworks.com/matlabcentral/answers/153414-estimator-standard-errors-using-fmincon-portfolio-optimization-context

disp('**********************************');
disp('   GMM estimates of RC Logit:');
disp('**********************************');
disp(['      Coeff','   ','Std Err']);
disp([theta1 standard_errors(1:Kbeta,1)]);
disp(['   Deviations','  ','Std Err']);
disp([theta2 standard_errors(Kbeta+1:Kbeta+Ktheta,1)]);

%% clear Q1
clear;
global XX y IV W tol_inner Kbeta Ktheta ns v share TM prods Total T g exp_delta_store prep_shareSum...
    Kgamma price z xi_jm_store xi_jm_star A_star omega_matrices mc_m m theta2 theta1

% This gets you started for the following models
% (1) Logit with and without IV
% (2) Nested Logit
% (3) Random Coefficients Logit


%%%%%%%%%%%%%%%%
%%% SETTINGS %%%
%%%%%%%%%%%%%%%%

DATA       = csvread('Data.csv');
IDmkt      = DATA(:,1);                 % Market identifier
IDprod     = DATA(:,2);                 % Product identifier
share      = DATA(:,3);                 % Market share
A          = DATA(:,4:6);               % Product characteristics
price      = DATA(:,7);                 % Price
z          = DATA(:,8:10);              % Instruments
% group      = DATA(:,11);                % Group identifier for Nested Logit
TM         = max(IDmkt);                % Number of markets
prods      = zeros(TM,1);               % # of products in each market
for m=1:TM
    prods(m,1) = max(IDprod(IDmkt==m,1));
end
T          = zeros(TM,2);               % w: T is the end row of products in each market
T(1,1)     = 1;
T(1,2)     = prods(1,1); 
for i=2:TM
    T(i,1) = T(i-1,2)+1;                % 1st Column market starting point
    T(i,2) = T(i,1)+prods(i,1)-1;       % 2nd Column market ending point
end
Total      = T(TM,2);                   % # of obsevations
TotalProd  = max(prods);                % Max # of products in a given market
% Ngroups    = max(group);                % # of groups

%% Q2 Jointly estimate demand and supply sides using RC logit 
% w: use y (not exp(y)) (y is DV of logit) as starting guess of delta
% dim 970*1000

tol_inner = 1.e-14;                         % Tolerance for inner loop (NFXP)
tol_outer = 1.e-6;                          % Tolerance for outer loop (Estimation)
randn('seed',1);                           % Reset normal random number generator
rand('seed',1);                            % Reset uniform random number generator
Kbeta   = 2+size(A,2);                      % # of parameters in mean utility
Ktheta  = 1+size(A,2);                      % # of parameters with random coefficient (w: # of SE to be estimated)

Kgamma = 1+size(z,2);                       % # of parameters for observed costs

% ns      = 1000;                             % # draws to simulate shares
ns      = 200;                             % # draws to simulate shares
v       = randn(Ktheta,ns);                 % Draws for share integrals during estimation
IV      = [ones(Total,1) A z A.^2 z.^2];    % Instruments
nIV     = size(IV,2);                       % # instrumental variables
% W       = (IV'*IV)\eye(nIV);                % Starting GMM weighting matrix
% x0      = rand(Kbeta+Ktheta,1);             % Starting values for all parameters

% prepare instruments for demand and supply side moments
IV(size(IV,1)+1:2*size(IV,1),size(IV,2)+1:size(IV,2)+1+size(z,2)*2)=[ones(Total,1) z z.^2];         % 2*Totalx(nIV+Kgamma)      1940*(13+7)
W       = (IV'*IV)\eye(size(IV,2));                % Starting GMM weighting matrix            % 2*nIVx(nIV+Kgamma)      (13+7)*(13+7)
x0      = rand(Kbeta+Ktheta+Kgamma,1);             % Starting values for all parameters

opts    = optimset('Display','iter','TolCon',1E-10,'TolFun',1E-10,'TolX',1E-10);
x_L     = [-Inf*ones(Kbeta,1);zeros(Ktheta,1);-Inf*ones(Kgamma,1)];     % Lower bounds is zero for standard deviations of random coefficients
x_U     = Inf.*ones(Kbeta+Ktheta+Kgamma,1);                % Upper bounds for standard deviations of random coefficients

%
sum_mktShares = zeros(TM,1);
for j=1:TM
    sum_mktShares(j) = sum(share(T(j,1):T(j,2)));
end
outsideProduct = 1- sum_mktShares;

outsideProduct_expand = [];
for j=1:TM
    temp = repmat(outsideProduct(j),prods(j),1);
    outsideProduct_expand = [outsideProduct_expand; temp];
end

n = size(DATA,1);
k = 5;

y = log(share./outsideProduct_expand);

XX = [A price];
exp_delta_store = exp(y);         % set initial values for exp_delta_store
prep_shareSum = zeros(TM,Total);
for r = 1:TM
    prep_shareSum(r,T(r,1):T(r,2))=1;
end
    
tic
[X,fval_rep,exitflag,output,lambda,grad,hessian] = fmincon(@GMMObjFun_Start,x0,[],[],[],[],x_L,x_U,[],opts);
toc
theta1 = X(1:Kbeta,1);
theta2 = X(Kbeta+1:Kbeta+Ktheta,1);
gamma = X(Kbeta+Ktheta+1:end,1);

standard_errors = sqrt(diag(inv(hessian)));         % standard errors incorrect? how to compute?
% ref https://www.mathworks.com/matlabcentral/answers/153414-estimator-standard-errors-using-fmincon-portfolio-optimization-context

disp('****************************************');
disp('RC Logit demand & supply joint estimates:');
disp('****************************************');
disp([' Demand coeff','   ','Std Err']);
disp([theta1 standard_errors(1:Kbeta,1)]);
disp(['Demand deviations','  ','Std Err']);
disp([theta2 standard_errors(Kbeta+1:Kbeta+Ktheta,1)]);
disp([' Supply coeff','  ','Std Err']);
disp([gamma standard_errors(Kbeta+Ktheta+1:end,1)]);

%% Q3 merger counterfactual analysis
% back out marginal costs based on demand and supply joint estimation results
alpha_hat = -2;
mc = price+1./(alpha_hat.*(1-share));                               % Total*1        970*1

exp_mu_ijm_hat = exp(XX*diag(theta2)*v);
s_jm_prep_mat_hat = exp_mu_ijm_hat.*(repmat(exp_delta_store,1,ns));
shareSum_hat = prep_shareSum*s_jm_prep_mat_hat;
shareSum_expand_hat = prep_shareSum'*shareSum_hat+1;
sijm = s_jm_prep_mat_hat./shareSum_expand_hat;                      % Total*ns       970*1000

alpha_i = alpha_hat + theta2(end).*v(end,:);                        % 1*ns           1*1000

% compute price elasticties matrices based on demand and supply joint estimates
eta_matrices = cell(TM,1);
for m = 1:TM
    total_products_market_m = prods(m);
    elasticities_mat = zeros(total_products_market_m,total_products_market_m);
    for j = 1:total_products_market_m
        for l = 1:total_products_market_m
            j_ID = find(IDmkt==m & IDprod==j);
            l_ID = find(IDmkt==m & IDprod==l);
            if j~=l
                elasticities_mat(j,l) = -price(l_ID)/share(j_ID)...
                    *mean(alpha_i.*sijm(j_ID,:).*sijm(l_ID,:),2);
            else
                elasticities_mat(j,l) = price(j_ID)/share(j_ID)...
                    *mean(alpha_i.*sijm(j_ID,:).*(1-sijm(j_ID,:)),2);
            end
        end
    end
    eta_matrices{m} = elasticities_mat;
end

% % test
% sum_own_price_elasticities = [];
% sum_cross_price_elasticities = [];
% for m = 1:TM
%     own_price_elasticities = diag(eta_matrices{m});
%     cross_price_elasticities = eta_matrices{m}-eye(prods(m)).*eta_matrices{m};
%     cross_price_elasticities = reshape(cross_price_elasticities, numel(cross_price_elasticities),1);
%     sum_own_price_elasticities = [sum_own_price_elasticities; own_price_elasticities];
%     sum_cross_price_elasticities = [sum_cross_price_elasticities; cross_price_elasticities];
% end
% elast1 = mean(sum_own_price_elasticities)      % mean own price elasticities across markets
% elast2 = mean(sum_cross_price_elasticities)    % mean cross price elasticities across markets


% compute new ownership structure
omega_star_matrices = cell(TM,1);
for m = 1:TM
    share_m = share(IDmkt==m);
    max_product1_ind = find(share_m==max(share_m));
    share_m(max_product1_ind)=0;
    max_product2_ind = find(share_m==max(share_m));
    omega_star_mat = diag(ones(size(share_m)));
    omega_star_mat(max_product1_ind,max_product2_ind)=1;
    omega_star_mat(max_product2_ind,max_product1_ind)=1;
    omega_star_matrices{m} = omega_star_mat;                % prods*prods
end

% compute capital omega matrices for finding counterfactual prices and shares
omega_matrices = cell(TM,1);
for m = 1:TM
    omega_star_mat_m = omega_star_matrices{m};
    eta_mat_m = eta_matrices{m};
    omega_matrices{m} = omega_star_mat_m.*(-eta_mat_m);     % prods*prods
end    

xi_jm_hat = xi_jm_store;                                    % Total*1        970*1          % get the estimated xi_jm from demand and supply estimation            

% compute new prices and mkt shares under FOC that satisfies profit maximization
p_star_store = cell(TM,1);
s_star_store = cell(TM,1);
profit_star_store = cell(TM,1);
profit_original_store = cell(TM,1);
sum_profit_star = zeros(TM,1);
sum_profit_original = zeros(TM,1);
change_consumer_surplus = zeros(TM,1);
change_welfare = zeros(TM,1);
for m = 1:TM
    % solve for new p that satisfies FOC for profit maximization 
    mc_m = mc(IDmkt==m);
    xi_jm_star = xi_jm_hat(IDmkt==m);
    A_star = A(IDmkt==m,:);
    p0=1.*ones(prods(m),1);
    opts = optimset('Display','off');
    p = fsolve(@MergerObjFun,p0,opts);                            % prods*1        % solve for new price under FOC
    p_star_store{m} = p;
    
    % compute new mkt shares under new price
    XX_star = [A_star p];
    exp_mu_ijm_star = exp(XX_star*diag(theta2)*v);           % prods*1000
    exp_delta_jm_star = exp([ones(prods(m),1) XX_star]*theta1 + xi_jm_star);    % prods*1
    s_jm_prep_mat_star = exp_mu_ijm_star.*(repmat(exp_delta_jm_star,1,ns));     % prods*1000
    sum_s_jm_star = sum(s_jm_prep_mat_star,1)+1;
    sum_s_jm_expand_star = repmat(sum_s_jm_star,prods(m),1);
    sjm_star = mean(s_jm_prep_mat_star./sum_s_jm_expand_star, 2);               % prods*1
    s_star_store{m} = sjm_star;
    
    % compute profit in counterfactual scenario
    profit_star_store{m} = (p-mc_m).*sjm_star;
    sum_profit_star(m) = sum(profit_star_store{m});
    
    % compute profit in baseline scenario
    profit_original_store{m} = (price(IDmkt==m)-mc_m).*share(IDmkt==m);
    sum_profit_original(m) = sum(profit_original_store{m});
    
    % compute consumer surplus in counterfactual scenario
    expect_CSi = 1/alpha_hat.*log(sum(s_jm_prep_mat_star,1));        % 1*1000            % plus an unobserved constant
    
    % compute consumer surplus in baseline scenario
    XX_original = [A_star price(IDmkt==m)];
    exp_mu_ijm_original = exp(XX_original*diag(theta2)*v);           % prods*1000
    exp_delta_jm_original = exp([ones(prods(m),1) XX_original]*theta1 + xi_jm_star);    % prods*1
    s_jm_prep_mat_original = exp_mu_ijm_original.*(repmat(exp_delta_jm_original,1,ns));     % prods*1000
    expect_CSi_original = 1/alpha_hat.*log(sum(s_jm_prep_mat_original,1));        % 1*1000            % plus an unobserved constant
    
    % compute change in consumer surplus
    change_consumer_surplus(m) = sum(1/alpha_hat.*(log(sum(s_jm_prep_mat_star,1))-log(sum(s_jm_prep_mat_original,1))));

    % compute change in total welfare
    change_welfare(m) = sum_profit_star(m)-sum_profit_original(m)+change_consumer_surplus(m);
    
end

% compute changes in different scenarios for merged and unmerged products
change_price_merged_products = zeros(TM,1);
change_share_merged_products = zeros(TM,1);
change_profit_merged_products = zeros(TM,1);
change_price_un_merged_products = zeros(TM,1);
change_share_un_merged_products = zeros(TM,1);
change_profit_un_merged_products = zeros(TM,1);
for m = 1:TM
    p_star = p_star_store{m};
    p_original = price(IDmkt==m);
    change_price = p_star-p_original;
    s_star = s_star_store{m};
    s_original = share(IDmkt==m);
    change_share = s_star-s_original;
    profit_star = profit_star_store{m};
    profit_original = profit_original_store{m};
    change_profit = profit_star-profit_original;
    merged_products = sum(omega_star_matrices{m},1)-1;
    un_merged_products = 1-merged_products;
    change_price_merged_products(m) = sum(change_price.*merged_products');
    change_share_merged_products(m) = sum(change_share.*merged_products');
    change_profit_merged_products(m) = sum(change_profit.*merged_products');
    change_price_un_merged_products(m) = sum(change_price.*un_merged_products');
    change_share_un_merged_products(m) = sum(change_share.*un_merged_products');
    change_profit_un_merged_products(m) = sum(change_profit.*un_merged_products');
end

% mean(change_price_merged_products)
% mean(change_share_merged_products)
% mean(change_profit_merged_products)
% mean(change_price_un_merged_products)
% mean(change_share_un_merged_products)
% mean(change_profit_un_merged_products)
% 
% mean(change_consumer_surplus)
% mean(change_welfare)

disp(['mean changes in prices for merged products across markets: ' num2str(mean(change_price_merged_products))])
disp(['mean changes in shares for merged products across markets: ' num2str(mean(change_share_merged_products))])
disp(['mean changes in profits for merged products across markets: ' num2str(mean(change_profit_merged_products))])
disp(['mean changes in prices for un-merged products across markets: ' num2str(mean(change_price_un_merged_products))])
disp(['mean changes in shares for un-merged products across markets: ' num2str(mean(change_share_un_merged_products))])
disp(['mean changes in profits for un-merged products across markets: ' num2str(mean(change_profit_un_merged_products))])

disp(['mean changes in consumer surplus across markets: ' num2str(mean(change_consumer_surplus))])
disp(['mean changes in total welfare across markets: ' num2str(mean(mean(change_welfare)))])

disp('Prices and profits for both merged and un-merged products drop, while market shares for both increases.')
disp('Both consumer surplus and total welfare drop significantly after merger.')
disp('w: The numbers seems unreliable given such heavy drops. There may be bugs somewhere in the computation which lead to incorrect results.')

