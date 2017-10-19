%% Static entry assignment
% by Weiyi Deng                            18 Oct 2017

clc
clear

%% Read Data (Code Generated Using MATLAB(R)'s `Import Data`)
filename = 'movieTheaters.txt';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%[^\n\r]';
% fileID = fopen(filename,'r');
%     dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ... 
%                         'TextType', 'string', 'EmptyValue', NaN, ...
%                         'HeaderLines' ,startRow-1, ...
%                         'ReturnOnError', false, 'EndOfLine', '\r\n');
fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ... 
                        'EmptyValue', NaN, ...
                        'HeaderLines' ,startRow-1, ...
                        'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
marketid = dataArray{:, 1};
year = dataArray{:, 2};
estcount = dataArray{:, 3};
population = dataArray{:, 4};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

Cr = population(year==2009)/1000;
Nr = grp2idx(estcount(year==2009))-1;

%% Q3
tau = 1/5*std(log(Cr));
pi = @(n) 0.9^(n-1);
check_n = 8;

log_Cr = mean(log(Cr));
Pr = mean(Nr)/check_n;

Pj = Pr;

% Nr = mean(N);
N = check_n;
 
kappa = exp(log_Cr + log(pi(1+(N-1)*Pr))-log(1+(N-1)*Pr)-tau*norminv(Pj));
disp(kappa)

%% Q4
% plot equilibrium 
p_vec = 0.01:0.01:1;
Prob = zeros(size(p_vec));
for i = 1:length(p_vec)
    Prob(i) = entryProb(log_Cr,p_vec(i),feval(pi,1+(N-1)*Pr),kappa,tau,N);
end  
plot(p_vec,Prob)
hold on
plot(p_vec,p_vec)
hold off

%% Q5
delta = 0.9;
Pr_equil = bayesNash(log_Cr,delta,kappa,tau,N);      % change pr also in pi !!
disp(Pr_equil)
disp(['Equilibrium entry probability Pr is ' num2str(Pr_equil)])

%% Q6
tau_vec = 0.0001:0.0001:tau;
Prob = zeros(size(tau_vec));

for i = 1:length(tau_vec)
%     tic
    Prob(i) = bayesNash(log_Cr,delta,kappa,tau_vec(i),N);
%     toc
end
plot(tau_vec,Prob)
disp(Prob(1))
disp(['When tau is close to zero, equilibrium entry probability Pr is ' num2str(Prob(1))])

%% Q9
rng(1)
check_r = 573;                                          % sample size
simu_log_Cr = mean(log(Cr))+std(log(Cr)).*randn(check_r,1);                % simulate log(Cr)
simu_Cr = exp(simu_log_Cr);

simu_Pr = zeros(check_r,1);                           
simu_Nr_prob = zeros(check_r,check_n+1);                                   % simulate Nr
for i = 1:check_r
    simu_Pr(i) = bayesNash(simu_log_Cr(i),delta,kappa,tau,check_n);
    for j = 1:check_n+1
        simu_Nr_prob(i,j) = factorial(check_n)/(factorial(j-1)*factorial(check_n-(j-1)))...
            *simu_Pr(i)^(j-1)*(1-simu_Pr(i))^(check_n-(j-1));
    end
end
temp_Nr_prob = (cumsum(simu_Nr_prob'))';    
draw_for_Nr=rand(check_r,1);
draw_for_Nr=repmat(draw_for_Nr,1,check_n+1);           
temp_Nr=temp_Nr_prob<draw_for_Nr;
simu_Nr=sum(temp_Nr,2);

%% Q11
x = simu_log_Cr;
y = simu_Nr;
% par0 =  [0.9 21.8 0.1]
par0 =  [0.5 10 0.5]

opts    = optimset('Display','iter');
[estimates,nll,exitflag,output,lambda,grad,hessian]=...
    fmincon(@(par)(LogLikBayesNash(par,x,y)),par0,[],[],[],[],...
                   [zeros(size(par0))],[2 Inf Inf],[],opts);
               
disp(['     estimates','   ','original values']);
disp(['delta   ' num2str(estimates(1)) '    ' num2str(delta)])
disp(['kappa   ' num2str(estimates(2)) '    ' num2str(kappa)])
disp(['tau     ' num2str(estimates(3)) '     ' num2str(tau)])
               
%% Q12
x = log(Cr);
y = Nr;
% par0 = [0.5 10 0.5];
par0 =  [0.9 21.8 0.1]
% par0 = [1 1 1];

opts    = optimset('Display','iter');
[estimates,nll,exitflag,output,lambda,grad,hessian]=...
    fmincon(@(par)(LogLikBayesNash(par,x,y)),par0,[],[],[],[],...
                   [zeros(size(par0))],[2 Inf Inf],[],opts);
               
disp(['     estimates','   ','original values']);
disp(['delta   ' num2str(estimates(1)) '          ' num2str(delta)])
disp(['kappa   ' num2str(estimates(2)) '    ' num2str(kappa)])
disp(['tau     ' num2str(estimates(3)) '     ' num2str(tau)])
disp('w: Estimates seem off, loglikelihood function incorrect ??');

disp(['pi = ' num2str(estimates(1)) '^(n-1)'])
