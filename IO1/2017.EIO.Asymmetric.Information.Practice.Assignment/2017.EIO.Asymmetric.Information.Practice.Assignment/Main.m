clc
clear

global DATA1 DATA2 v_IF rep_v_IF NF NI_F NS NJ NObs NI NT 

DATA1       = csvread('Data1.csv');
DATA2       = csvread('Data2.csv');

NF = 5;                    % number of emplyer firms f
NI_F = 500;                % number of individuals working for each firm f
NJ = 5;                    % number of plans to choose
NI = NF*NI_F;              % total number of individuals in data
NObs = size(DATA1,1);      % total number of observations in data
NT = 5;                    % number of years

NS = 200;                  % number of simulation draws

rng(1)
v_IF = randn(NI,NS);       % simulate unobserved risks with NS draws for each individual 

% test
% A = [1 2 3; 4 5 6]
% AA = repmat(A(:),1,2)'
% reshape(AA,[4,3])

v_IF_temp = repmat(v_IF(:),[1 30])';              
rep_v_IF = reshape(v_IF_temp,[75000, NS]);       % replicate v_IF to have the same # of rows as in DATA1

b0 = [1   1   1  1  1];
opts    = optimset('display','iter-detailed','Diagnostics','on','TolFun',1e-10,'TolX',1e-10,'GradObj','off','Hessian','off');
tic
[b, fval,exitflag,output,grad,hessian] = fminunc(@SML_ObjectFun,b0,opts);
toc
se      = sqrt(diag(hessian).^-1);

disp('*************************');
disp('    SML estimates:');
disp('*************************');
disp(['    Coeff','     ','Std Err']);
disp([b' se]);

