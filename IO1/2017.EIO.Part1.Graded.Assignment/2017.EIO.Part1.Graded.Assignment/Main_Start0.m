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
[X,fval_rep,exitflag,output,lambda,grad,hessian] = fmincon(@GMMObjFun_Start,x0,[],[],[],[],x_L,x_U,[],opts);
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

