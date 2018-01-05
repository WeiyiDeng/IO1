%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMPIRICAL INDUSTRIAL ORGANIZATION (230319) - PART 1 - DEMAND ESTIMATION %
% by Weiyi Deng (revised version, much faster!)         Oct 12 2017
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
group      = DATA(:,11);                % Group identifier for Nested Logit
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
Ngroups    = max(group);                % # of groups


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

X = [ones(n,1) A price];

b_ols = inv(X'*X)*X'*y;

% sigma_squared = 1/(n-k)*sum((y-X*b_ols).^2);           % ???
% Var_b_ols = sigma_squared*inv(X'*X)
% std_b_ols = sqrt(diag(Var_b_ols))

std_b_ols = sqrt(diag(mean((y-X*b_ols).^2)*((X'*X)\eye(size(X,2)))));

disp('*************************');
disp('OLS estimates:');
disp('*************************');
disp(['    Coeff','     ','Std Err']);
disp([b_ols std_b_ols]);

RSS    = sum((y-X*b_ols).^2,1);
TSS    = sum((y-mean(y,1)).^2,1);
Rsq    = 1 - RSS/TSS
AdjRsq = 1-(1-Rsq)*(Total-1)/(Total-size(b_ols,1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OLS ESTIMATION OF FIRST STAGE IV %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1st stage IV:
X_stage1 = [ones(n,1) A z];
b_stage1 = inv(X_stage1'*X_stage1)*X_stage1'*price;
price_hat = X_stage1*b_stage1;

std_b_stage1 = sqrt(diag(mean((price-X_stage1*b_stage1).^2)*((X_stage1'*X_stage1)\eye(size(X_stage1,2)))));

disp('*************************');
disp('IV 1st stage estimates:');
disp('*************************');
disp(['    Coeff','     ','Std Err']);
disp([b_stage1 std_b_stage1]);

RSS    = sum((price-X_stage1*b_stage1).^2,1);
TSS    = sum((price-mean(price,1)).^2,1);
Rsq    = 1 - RSS/TSS
AdjRsq = 1-(1-Rsq)*(Total-1)/(Total-size(b_stage1,1))

% F-statistic for instruments:
X_exog = [ones(n,1) A];
RIVfs  = (X_exog'*X_exog)\(X_exog'*price);
RRSS   = sum((price-X_exog*RIVfs).^2,1);
F      = ((RRSS-RSS)/(size(b_stage1,1)-size(RIVfs,1)))/(RSS/(Total-size(b_stage1,1)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2SLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2nd stage IV:
X_stage2 = [ones(n,1) A price_hat];
b_stage2 = (X_stage2'*X_stage2)\(X_stage2'*y);

% w: or use 2sls estimator directly (see Ectrics 2 Lecture5 slide 21 & 25) 
XX = [ones(n,1) A price];
IV = [ones(n,1) A z];
PZ       = IV*inv(IV'*IV)*IV';
Xhat     = PZ*XX;
beta2sls = (Xhat'*Xhat)\(Xhat'*y);
se2sls   = sqrt(diag(mean((y-XX*beta2sls).^2)*((XX'*PZ*XX)\eye(size(XX,2)))));   % same as inv here

disp('*************************');
disp('IV 2SLS estimates:');
disp('*************************');
disp(['    Coeff','     ','Std Err']);
disp([beta2sls se2sls]);

% absolute value of 2SLS estimated price coefficient increased than OLS, IV
% method helps solve the endogeneity issue, but SE values also increased
% as 2SLS gives less precise estimates than OLS ?

% compute price elasticties matrices
price_coef = beta2sls(5);
eta_matrices = cell(TM,1);
for m = 1:TM
    total_products_market_m = prods(m);
    elasticities_mat = zeros(total_products_market_m,total_products_market_m);
    for j = 1:total_products_market_m
        for l = 1:total_products_market_m
            j_ID = find(IDmkt==m & IDprod==j);
            l_ID = find(IDmkt==m & IDprod==l);
            if j~=l
                elasticities_mat(j,l) = -price_coef*price(l_ID)*share(l_ID);
            else
                elasticities_mat(j,l) = price_coef*price(j_ID)*(1-share(j_ID));
            end
        end
    end
    eta_matrices{m} = elasticities_mat;
end

sum_own_price_elasticities = [];
sum_cross_price_elasticities = [];
for m = 1:TM
    own_price_elasticities = diag(eta_matrices{m});
    cross_price_elasticities = eta_matrices{m}-eye(prods(m)).*eta_matrices{m};
    cross_price_elasticities = reshape(cross_price_elasticities, numel(cross_price_elasticities),1);
    sum_own_price_elasticities = [sum_own_price_elasticities; own_price_elasticities];
    sum_cross_price_elasticities = [sum_cross_price_elasticities; cross_price_elasticities];
end
elast1 = mean(sum_own_price_elasticities);      % mean own price elasticities across markets
elast2 = mean(sum_cross_price_elasticities);    % mean cross price elasticities across markets

disp(['mean own price elasticities across markets: ' num2str(elast1)])
disp(['mean cross price elasticities across markets: ' num2str(elast2)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2SLS ESTIMATION OF NESTED LOGIT MODEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % test groups
% store_num_groups = [];
% for m = 1:TM
%     groups_mkt_m = group(find(IDmkt==m));
%     num_groups_m = length(unique(groups_mkt_m));
%     store_num_groups = [store_num_groups num_groups_m];
% end
% mean(store_num_groups)==Ngroups
% % each market has exactly 3 groups

group_shares = zeros(TM,Ngroups);
for m = 1:TM
    m_share_rows = share(T(m,1):T(m,2));
    m_group_rows = group(T(m,1):T(m,2));
    for gp = 1:Ngroups
        g_group_share_m = sum(m_share_rows(m_group_rows==gp));
        group_shares(m,gp) = g_group_share_m;
    end
end

group_shares_vec_temp = zeros(Total,1);
for i = 1:Total
    group_shares_vec_temp(i) = group_shares(IDmkt(i),group(i));
end
withinGroup_share = log(share./group_shares_vec_temp);        

% 2sls estimates
XX = [ones(n,1) A price withinGroup_share];
IV = [ones(n,1) A z];
PZ       = IV*inv(IV'*IV)*IV';
Xhat     = PZ*XX;
beta2sls = (Xhat'*Xhat)\(Xhat'*y);
se2sls   = sqrt(diag(mean((y-XX*beta2sls).^2)*((XX'*PZ*XX)\eye(size(XX,2)))));   % same as inv here    

disp('**********************************');
disp('IV 2SLS estimates of nested logit:');
disp('**********************************');
disp(['    Coeff','     ','Std Err']);
disp([beta2sls se2sls]);

% group component not statistically significant
% price coefficient more biased (or less precise) ?

price_coef = beta2sls(5);
sigma_coef = beta2sls(6);

% compute price elasticties matrices
eta_matrices = cell(TM,1);
for m = 1:TM
    total_products_market_m = prods(m);
    elasticities_mat = zeros(total_products_market_m,total_products_market_m);
    for j = 1:total_products_market_m
        for l = 1:total_products_market_m
            j_ID = find(IDmkt==m & IDprod==j);
            l_ID = find(IDmkt==m & IDprod==l);
            if j==l                                    % own price elasticity
                elasticities_mat(j,l) = price_coef*price(j_ID)/(1-sigma_coef)*...
                    (1-sigma_coef*exp(withinGroup_share(j_ID))-(1-sigma_coef)*share(j_ID));
            elseif group(j_ID)==group(l_ID)            % cross price elasticities same group
                elasticities_mat(j,l) = -price_coef*price(l_ID)*share(l_ID);
            else                                       % cross price elasticities different groups
                elasticities_mat(j,l) = -price_coef*price(l_ID)/(1-sigma_coef)*...
                    (sigma_coef*exp(withinGroup_share(l_ID))+(1-sigma_coef)*share(l_ID));
            end
        end
    end
    eta_matrices{m} = elasticities_mat;
end

sum_own_price_elasticities = [];
sum_cross_price_elasticities = [];
for m = 1:TM
    own_price_elasticities = diag(eta_matrices{m});
    cross_price_elasticities = eta_matrices{m}-eye(prods(m)).*eta_matrices{m};
    cross_price_elasticities = reshape(cross_price_elasticities, numel(cross_price_elasticities),1);
    sum_own_price_elasticities = [sum_own_price_elasticities; own_price_elasticities];
    sum_cross_price_elasticities = [sum_cross_price_elasticities; cross_price_elasticities];
end
elast1 = mean(sum_own_price_elasticities);      % mean own price elasticities across markets
elast2 = mean(sum_cross_price_elasticities);    % mean cross price elasticities across markets

disp(['mean own price elasticities across markets: ' num2str(elast1)])
disp(['mean cross price elasticities across markets: ' num2str(elast2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BLP ESTIMATION OF RANDOM COEFFICIENTS LOGIT MODEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% w: use y (not exp(y)) (y is DV of logit) as starting guess of delta
% dim 970*1000

tol_inner = 1.e-14;                         % Tolerance for inner loop (NFXP)
tol_outer = 1.e-6;                          % Tolerance for outer loop (Estimation)
randn('seed',10);                           % Reset normal random number generator
rand('seed',10);                            % Reset uniform random number generator
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
[X,fval_rep,exitflag,output,lambda,grad,hessian] = fmincon(@GMMObjFun_Start_revised,x0,[],[],[],[],x_L,x_U,[],opts);
toc
theta1 = X(1:Kbeta,1);
theta2 = X(Kbeta+1:Kbeta+Ktheta,1);

standard_errors = sqrt(diag(inv(hessian)));         % standard errors incorrect? how to compute?
% ref https://www.mathworks.com/matlabcentral/answers/153414-estimator-standard-errors-using-fmincon-portfolio-optimization-context

disp('**********************************');
disp('IV GMM estimates of RC Logit:');
disp('**********************************');
disp(['      Coeff','   ','Std Err']);
disp([theta1 standard_errors(1:Kbeta,1)]);
disp(['   Deviations','  ','Std Err']);
disp([theta2 standard_errors(Kbeta+1:Kbeta+Ktheta,1)]);

