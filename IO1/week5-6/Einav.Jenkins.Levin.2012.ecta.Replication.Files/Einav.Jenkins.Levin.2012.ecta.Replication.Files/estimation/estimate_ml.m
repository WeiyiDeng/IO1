%--------------------------------------------------------------------------
% estimate_ml.m
%--------------------------------------------------------------------------

clear all;
diary output.log;

global N X1 X2 X3 X4 K1 K2 K3 K4;
global I less more pK lK dK dK2 tK afK;
global fracpaid none obsdef cenind cenpoint G;
global rho;

load 'full_data';
myoptions_simplex = optimset('MaxIter',1000000,'TolFun',1e-3,'TolX',1e-3,'MaxFunEvals',1000000,'Display','final');
myoptions_initial = optimset('MaxIter',1000000,'TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',1000000,'Display','iter');

%--------------------------------------------------------------------------
% Estimate full demand model
%--------------------------------------------------------------------------

par = load('par_stata.txt');
par = fminunc('likelihood',par,myoptions_initial);
fid = fopen('par_newton.txt','wt');
fprintf(fid,'%10.6f\n',par);
fclose(fid);

par = load('par_newton.txt');
par = fminsearch('likelihood',par,myoptions_simplex);
fid = fopen('par_simplex.txt','wt');
fprintf(fid,'%10.6f\n',par);
fclose(fid);

%--------------------------------------------------------------------------
% end of function
%--------------------------------------------------------------------------
