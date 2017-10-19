%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  EMPIRICAL INDUSTRIAL ORGANIZATION (230319) - PART 1 - MATLAB PRACTICE  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code presents 3 types of estimation examples:
% (1) Instrumental variables (IV)
% (2) Seemingly unrelated regression (SUR)
% (3) Maximum likelihood (MLE)

% For each example it does the following:
% - Generates the data and the true parameters
% - Estimates different versions of the models
% - Calculates standard errors


%% (1) INSTRUMENTAL VARIABLES ESTIMATION %%

% This section looks at two types of IV estimation:
% (A) Just-identified case
% (B) Over-identified case

% For each of these, it recovers:
% - First stage OLS, with standard errors and Rsq, AdjRsq, F-test
% - Second stage IV with standard errors
% - Second stage OLS with standard errors


% GENERATE THE DATASET:

% Generate variables:
clc;
clear;
obs    = 10000;
randn('seed',1);
X      = [ones(obs,1) randn(obs,2)];
mu     = [0 0 0 0];
sigma  = [1.0 0.5 0.0 0.0;...
          0.5 1.0 0.5 0.5;...
          0.0 0.5 1.0 0.0;...
          0.0 0.5 0.0 1.0];
W      = mvnrnd(mu,sigma,obs);
eps    = W(:,1);
Y2     = W(:,2);
Z1     = W(:,3);
Z2     = W(:,4);             % set IVs not correlated with epsilon
disp(['Covariance (Y2,eps) = ',num2str(Y2'*eps/obs)]);
disp(['Covariance (Z1,Y2) = ',num2str(Z1'*Y2/obs)]);
disp(['Covariance (Z2,Y2) = ',num2str(Z2'*Y2/obs)]);
    
% Set true parameters:
beta   = [1 1.5 -1 -1.5];
XX     = [X Y2];
Y1     = XX*beta' + eps;
clear W eps mu sigma

% Export the data:
DATA = [X Y1 Y2 Z1 Z2];
csvwrite('PS0_Data.1.csv',DATA);
clear DATA



% (A) JUST-IDENTIFIED CASE (question A.5 in problem set):
IV     = [X Z1];

% First stage IV:
IVfs   = (IV'*IV)\(IV'*Y2);
seIVfs = sqrt(diag(mean((Y2-IV*IVfs).^2)*((IV'*IV)\eye(size(IV,2)))));
str    = [IVfs seIVfs];
RSS    = sum((Y2-IV*IVfs).^2,1);
TSS    = sum((Y2-mean(Y2,1)).^2,1);
Rsq    = 1 - RSS/TSS;
AdjRsq = 1-(1-Rsq)*(obs-1)/(obs-size(IVfs,1));

% F-statistic for instruments:
RIVfs  = (X'*X)\(X'*Y2);
RRSS   = sum((Y2-X*RIVfs).^2,1);
F      = ((RRSS-RSS)/(size(IVfs,1)-size(RIVfs,1)))/(RSS/(obs-size(IVfs,1)));

disp('*************************');
disp('First stage IV estimates:');
disp('*************************');
disp(['    Coeff','     ','Std Err']);
disp(str);
disp(['Rsq:    ',num2str(Rsq)]);
disp(['AdjRsq: ',num2str(AdjRsq)]);
disp(['F-stat: ',num2str(F)]);


% Second stage IV:
beta2sls = (IV'*XX)\(IV'*Y1);
se2sls   = sqrt(diag(mean((Y1-XX*beta2sls).^2)*((IV'*XX)\eye(size(XX,2)))));
str2     = [beta2sls se2sls beta'];


% Second stage OLS:
betaols  = (XX'*XX)\(XX'*Y1);
seols    = sqrt(diag(mean((Y1-XX*betaols).^2)*((XX'*XX)\eye(size(XX,2)))));
str3     = [betaols seols beta'];

disp('*************************');
disp('*************************');
disp('Second stage IV estimates:');
disp('*************************');
disp(['    Coeff','    ','Std Err','    ','True']);
disp(str2);
disp('*************************');
disp('Second stage OLS estimates:');
disp('*************************');
disp(['    Coeff','    ','Std Err','    ','True']);
disp(str3);



% (B) OVER-IDENTIFIED CASE (question B.3 in problem set):
IV     = [X Z1 Z2];

% First stage IV:
IVfs   = (IV'*IV)\(IV'*Y2);
seIVfs = sqrt(diag(mean((Y2-IV*IVfs).^2)*((IV'*IV)\eye(size(IV,2)))));
str    = [IVfs seIVfs];
RSS    = sum((Y2-IV*IVfs).^2,1);
TSS    = sum((Y2-mean(Y2,1)).^2,1);
Rsq    = 1 - RSS/TSS;
AdjRsq = 1-(1-Rsq)*(obs-1)/(obs-size(IVfs,1));

% F-statistic for first stage IV:
RIVfs  = (X'*X)\(X'*Y2);
RRSS   = sum((Y2-X*RIVfs).^2,1);
F      = ((RRSS-RSS)/(size(IVfs,1)-size(RIVfs,1)))/(RSS/(obs-size(IVfs,1)));

disp('*************************');
disp('First stage IV estimates:');
disp('*************************');
disp(['    Coeff','     ','Std Err']);
disp(str);
disp(['Rsq:    ',num2str(Rsq)]);
disp(['AdjRsq: ',num2str(AdjRsq)]);
disp(['F-stat: ',num2str(F)]);


% Secons stage IV:
W        = inv(IV'*IV);
PZ       = IV*W*IV';
Xhat     = PZ*XX;
beta2sls = (Xhat'*Xhat)\(Xhat'*Y1);
se2sls   = sqrt(diag(mean((Y1-XX*beta2sls).^2)*((XX'*PZ*XX)\eye(size(XX,2)))));
str2     = [beta2sls se2sls beta'];

% Formula for robust se (for heteroskedasticity):
%se2slsR = sqrt(diag(inv(Xhat'*Xhat)*(Xhat'*diag((Y1-XX*beta2sls).^2)*Xhat)*inv(Xhat'*Xhat)));


% Second stage OLS:              % w: upward bias in ols endog estimation
betaols  = (XX'*XX)\(XX'*Y1);
seols    = sqrt(diag(mean((Y1-XX*betaols).^2)*((XX'*XX)\eye(size(XX,2)))));
str3     = [betaols seols beta'];

% Formula for robust se (for heteroskedasticity):
%seolsR  = sqrt(diag(inv(XX'*XX)*(XX'*diag((Y1-XX*betaols).^2)*XX)*inv(XX'*XX)));

disp('*************************');
disp('*************************');
disp('Second stage IV estimates:');
disp('*************************');
disp(['    Coeff','    ','Std Err','    ','True']);
disp(str2);
disp('*************************');
disp('Second stage OLS estimates:');
disp('*************************');
disp(['    Coeff','    ','Std Err','    ','True']);
disp(str3);



%% (2) SEEMINGLY UNRELATED REGRESSION %%

% This section estimates the same model with:
% (A) OLS regressions
% (B) Seemingly unrelated regressions


% GENERATE THE DATASET

% Generate variables:
clear;
obs    = 10000;
randn('seed',1);
X1     = [ones(obs,1) randn(obs,2)];
X2     = [ones(obs,1) randn(obs,2)];
mu     = [0 0];
sigma  = [1.0 0.5;...
          0.5 1.0];
W      = mvnrnd(mu,sigma,obs);
eps1   = W(:,1);
eps2   = W(:,2);
disp(['Covariance (eps1,eps2) = ',num2str(eps1'*eps2/obs)]);
    
% Set true parameters:
beta1  = [1 1.5 -1];
beta2  = [-1 1 1.5];
Y1     = X1*beta1' + eps1;
Y2     = X2*beta2' + eps2;
clear W eps1 eps2 mu sigma

% Export the data:
DATA = [X1 X2 Y1 Y2];
csvwrite('PS0_Data.2.csv',DATA);
clear DATA


% Estimate separately by OLS
beta1ols = (X1'*X1)\(X1'*Y1);
se1ols   = sqrt(diag(mean((Y1-X1*beta1ols).^2)*((X1'*X1)\eye(size(X1,2)))));
str1     = [beta1ols se1ols beta1'];
beta2ols = (X2'*X2)\(X2'*Y2);
se2ols   = sqrt(diag(mean((Y2-X2*beta2ols).^2)*((X2'*X2)\eye(size(X2,2)))));
str2     = [beta2ols se2ols beta2'];

disp('*************************');
disp('     OLS estimates:     ');
disp('*************************');
disp(['    Coeff','    ','Std Err','    ','True']);
disp(str1);
disp('*************************');
disp(str2);


% Estimate jointly by SUR
Sigma = mean((Y1 - X1*beta1ols).*(Y2 - X2*beta2ols),1);
W     = kron(inv(Sigma),eye(obs));

beta1sur = (X1'*W*X1)\(X1'*W*Y1);
se1sur   = sqrt(diag(mean((Y1-X1*beta1sur).^2)*((X1'*W*X1)\eye(size(X1,2)))));
str1     = [beta1sur se1sur beta1'];
beta2sur = (X2'*W*X2)\(X2'*W*Y2);
se2sur   = sqrt(diag(mean((Y2-X2*beta2sur).^2)*((X2'*W*X2)\eye(size(X2,2)))));
str2     = [beta2sur se2sur beta2'];
disp('*************************');
disp('     SUR estimates:     ');
disp('*************************');
disp(['    Coeff','    ','Std Err','    ','True']);
disp(str1);
disp('*************************');
disp(str2);



%% (3.a) MAXIMUM LIKELIHOOD ESTIMATION %%

% This section estimates a Probit model by MLE


% GENERATE THE DATASET:

% Generate variables:
clear;
global Y X              % Declare global variables (carried across m-files)
obs  = 10000;
rand('seed',1);
randn('seed',1);
X    = [ones(obs,1) rand(obs,3)];
eps  = randn(obs,1);
%eps = -evrnd(0,1,obs,1);                   % In case you want to do Logit

% Set true parameters:
beta = [1 1.5 -1 -1.5];
Ystr = X*beta' + eps;
Y    = zeros(obs,1);
Y(Ystr>0,1) = 1;
clear eps Ystr

% Export the data:
DATA = [Y X];
csvwrite('PS0_Data.3.csv',DATA);
clear DATA


% Estimate Probit by MLE:

opts    = optimset('display','iter-detailed','Diagnostics','on','TolFun',1e-10,'TolX',1e-10,'GradObj','on','DerivativeCheck','off','Hessian','user-supplied');
start   = (X'*X)\X'*Y;
tic
[betahat,~,~,output,grad,hessian] = fminunc(@PS0_llike,start,opts);
toc
se      = sqrt(diag(hessian).^-1);
str1 = [betahat se beta'];
disp('*************************');
disp(['    Coeff','     ','Std Err','    ','True']);
disp(str1);

