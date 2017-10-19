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

C = population(year==2009)/1000;
N = grp2idx(estcount(year==2009))-1;

%% Q3
tau = 1/5*std(log(C));
pi = @(n) 0.9^(n-1);
check_n = 8;

log_Cr = mean(log(C));
Pr = mean(N)/check_n;

Pj = Pr;

Nr = mean(N);
 
kappa = exp(log_Cr + log(pi(1+(Nr-1)*Pr))-log(1+(Nr-1)*Pr)-tau*norminv(Pj))

% plot equilibrium 
p_vec = 0.01:0.01:1;
Prob = zeros(size(p_vec));
for i = 1:length(p_vec)
    Prob(i) = entryProb(mean(C),p_vec(i),feval(pi,1+(Nr-1)*Pr),kappa,tau,mean(N));
end  
plot(p_vec,Prob)
hold on
plot(p_vec,p_vec)
hold off

