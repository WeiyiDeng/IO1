%% Estimating an Ordered Probit Entry Model
% Jaap Abbring, 21 July 2017

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

%% MLE of Ordered Probit Model with Log Population (in Thousands)
x = log(population(year==2009)/1000);
y = grp2idx(estcount(year==2009));
parZero = [2.92 0.1*ones(1,max(y)-2) 1];
[estimates,nll,exitflag,output,lambda,grad,hessian]=...
    fmincon(@(par)(-logLikOProbit(par,x,y)),parZero,[],[],[],[],...
                   [-inf zeros(1,max(y)-1)],[]);

standardErrors = sqrt(diag(inv(hessian)));

disp('MLE (s.e.)');
disp([estimates' standardErrors]);

%% Implied Toughness of Competition 
popThresholds = exp(cumsum(estimates(1:end-1)));
popPerFirm = popThresholds./(1:max(y)-1);
disp('entry tresholds in 1000s of population (total/per firm)');
disp([popThresholds' popPerFirm']); 

toughness = popPerFirm(1:end-1)./popPerFirm(2:end)
