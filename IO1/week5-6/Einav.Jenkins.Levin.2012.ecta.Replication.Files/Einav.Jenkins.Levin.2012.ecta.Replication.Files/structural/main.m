%--------------------------------------------------------------------------
% main.m
%--------------------------------------------------------------------------
%
% Solves the borrower's optimal down and repayment problem for a given
% parameter vector, computes moments, and draws scatter plot
%
%--------------------------------------------------------------------------
% variable declarations
%--------------------------------------------------------------------------
clear;

%parameters
global B T G Nsim mindwn ptilda dscale pscale deltaD deltaP APRUSE TERM;
global rho psi cval delta beta betap rra theta;
global saleeps downeps epsilon up ud;
global ystep ymax ymin ystates ynum y;
global lstep lmax lmin lstates lnum l;
global vstep vmax vmin vstates vnum v;

%--------------------------------------------------------------------------
% define discretized state space
%--------------------------------------------------------------------------
ymin=0.050; ystep=0.050; ymax=10.00; ystates=(ymin:ystep:ymax); ynum=size(ystates,2); %income
lmin=0.100; lstep=0.100; lmax=2.000; lstates=(lmin:lstep:lmax); lnum=size(lstates,2); %payment
vmin=-3.00; vstep=0.600; vmax=6.000; vstates=(vmin:vstep:vmax); vnum=size(vstates,2); %value
%--------------------------------------------------------------------------
% define parameters
%--------------------------------------------------------------------------
T = 40; G = T; B = 20; 
dscale=5; pscale=20; deltaD=0.100*dscale; deltaP=0.100*pscale;

%samle offer terms
mindwn=1.000;       %average (rounded)
ptilda=11.000;      %average (rounded)
APRUSE=29.9;        %mode
TERM=42.0;          %mode

%imposed parameter values
delta = 0.880^((TERM/T)/12); 
beta  = 0.750^((TERM/T)/12); 
betap = beta; 
theta = 1.000; 
rra   = 0.999; 
rho   = 1.000; 
cval  = 0.000;
psi   = 0.000;

%calibrated parameter values
vavg = 2.9309;
vvar = 1.5348;
davg = 0.5012;
dvar = 2.0616;
corr =-0.0247;
sig2 = 0.2954;
phi  = 1.0295;

%--------------------------------------------------------------------------
% consumption utilities (YxLxV)
%--------------------------------------------------------------------------
y = repmat(ystates',[1,lnum,vnum]);
l = repmat(lstates ,[ynum,1,vnum]);
v = ones(1,1,vnum); v(1,1,:)=vstates; v=repmat(v,[ynum,lnum,1]);
up = CalcUtility(y,l); ud = CalcUtility(y,0);

%--------------------------------------------------------------------------
% simulate unobservables
%--------------------------------------------------------------------------
Nsim = 1*10^5;
stream = RandStream('mt19937ar','Seed',1); RandStream.setDefaultStream(stream); saleeps = randn(Nsim,1); %outside option
stream = RandStream('mt19937ar','Seed',2); RandStream.setDefaultStream(stream); downeps = randn(Nsim,1); %initial liquidity
stream = RandStream('mt19937ar','Seed',3); RandStream.setDefaultStream(stream); epsilon = randn(Nsim,T); %liquidity shocks

%--------------------------------------------------------------------------
% solve borrower's optimal purchase and repayment problem
%--------------------------------------------------------------------------

%simulate outside option (vout) and liquidity (y0,...,yT) variables
[v0 vi y0 yi yindex yt] = simulate_unobservables(vavg,vvar,davg,dvar,corr,sig2);

%solve repayment problem for all states (note: independent of d,p)
[P EVP EVD D] = optimal_repayment(sig2,delta,beta,psi,phi,theta);

%calculate optimal outcomes for each borrower, given vout and y0
[dopt  lopt  uout ubuy  uopt  simsale  simless  simxtra  simdef  simfrac ] = optimal_down(mindwn,ptilda,v0,vi,y0,yi,yindex,EVP,EVD,D,cval,betap,theta,beta,phi);
[dopt2 lopt2 uout ubuy2 uopt2 simsale2 simless2 simxtra2 simdef2 simfrac2] = optimal_down(mindwn+deltaD,ptilda,v0,vi,y0,yi,yindex,EVP,EVD,D,cval,betap,theta,beta,phi);
[dopt3 lopt3 uout ubuy3 uopt3 simsale3 simless3 simxtra3 simdef3 simfrac3] = optimal_down(mindwn,ptilda+deltaP,v0,vi,y0,yi,yindex,EVP,EVD,D,cval,betap,theta,beta,phi);

%compute moments and print for Table A2
[MOMENTS SIMMOMENTS OBJFN] = CalcMoments(simsale,simless,simxtra,simdef,simfrac,simsale2,simsale3);
[MOMENTS SIMMOMENTS]

%--------------------------------------------------------------------------
% Generate histograms for Figure A1
%--------------------------------------------------------------------------

%output data for down payment histogram
s = simsale;
EDGES = [-Inf 0.00001:0.2:1.80001 Inf];
modhist = histc(simxtra(s),EDGES)/sum(s);
downout = [modhist(1:end-1); 1-mean(s)]

%output data for default time histogram
simpmts = simfrac*T;
EDGES = -0.5:2.0:41.5;
modhist = histc(simpmts,EDGES)/sum(s);
fracout = modhist

%--------------------------------------------------------------------------
% Generate scatter plot
%--------------------------------------------------------------------------
%{
%create curves
yhi=6; use = y0>-2 & y0<yhi & v0<6 & v0>-2; use=y0>-99999;
a = abs(ubuy-uout)<0.2; b = dopt>mindwn & abs(dopt-mindwn)<0.005;
temp = [y0(a) v0(a) use(a)]; temp = sortrows(temp,1);
yforcurve1 = temp(:,1); vthreshold = temp(:,2); use1 = temp(:,3)==1;
temp = [y0(b) v0(b) use(b)]; temp = sortrows(temp,1);
yforcurve2 = temp(:,1); xthreshold = temp(:,2); use2 = temp(:,3)==1;

%scatter plots
hold on;
scatter(y0(use & ~simsale),v0(use & ~simsale),2)
hold on;
scatter(y0(use & simless),v0(use & simless),2)
hold on;
scatter(y0(use & simsale & ~simless),v0(use & simsale & ~simless),2)

%threshold curves
hold on;
line(yforcurve1(use1),vthreshold(use1))
hold on;
line(yforcurve2(use2),xthreshold(use2))

%output data for scatter plot
points = [y0 v0 simsale simless];
fid = fopen('scatter.txt','wt');
for i=1:Nsim,
    fprintf(fid,'%6.3f  ',points(i,:));
    fprintf(fid,'\n');
end;
fclose(fid);
%}
%--------------------------------------------------------------------------
% Figure A3(a): purchase utility and optimal down versus price
%--------------------------------------------------------------------------
%{
numpts = 21;
values = zeros(numpts,5);
simulated_apps = [];
simulated_sales = [];
for i=1:numpts,
    
    %calculate simulated outcomes conditional on sale
    [dopt lopt uout ubuy uopt simsale simless simxtra simdef simfrac] = optimal_down(mindwn,i,v0,vi,y0,yi,yindex,EVP,EVD,D,cval,betap,theta,beta,phi);
    x1 = mean(simsale);
    x2 = mean(simless(simsale));
    x3 = mean(dopt(simsale));
    x4 = mean(simfrac);
    x5 = mean(simdef);
    values(i,:)=[x1 x2 x3 x4 x5];
    %{
    %caluclate simulated outcomes not conditional on sale
    temp_d = -99; %no minimum down restriction
    [dopt lopt uout ubuy uopt simsale simless simxtra simdef simfrac] = optimal_down(temp_d,i,v0,vi,y0,yi,yindex,EVP,EVD,D,cval,betap,theta,beta,phi);
    temp = dopt; temp(temp<0)=0;
    x6 = mean(temp);
    x7 = mean(simfrac);
    x8 = mean(simfrac);
    values(i,:)=[x1 x2 x3 x4 x5 x6 x7 x8];
    %}
    %create data for reduced-form regressions
    new_apps = [y0 i*ones(Nsim,1) mindwn*ones(Nsim,1) simsale simless];
    new_sales = [new_apps(simsale,:) dopt(simsale) simdef simfrac];
    simulated_apps = [simulated_apps; new_apps];
    simulated_sales = [simulated_sales; new_sales];
    
end;
asdf
%--------------------------------------------------------------------------
% Output data
%--------------------------------------------------------------------------
fid = fopen('simoutcomes_vs_price.txt','wt');
for i=1:numpts,
    fprintf(fid,'%10.6f ',values(i,:));
    fprintf(fid,'\n');
end;
fclose(fid);
fid = fopen('simulated_apps_p.txt','wt');
for i=1:size(simulated_apps,1),
    fprintf(fid,'%10.6f,',simulated_apps(i,:));
    fprintf(fid,'\n');
end;
fclose(fid);
fid = fopen('simulated_sales_p.txt','wt');
for i=1:size(simulated_sales,1),
    fprintf(fid,'%10.6f,',simulated_sales(i,:));
    fprintf(fid,'\n');
end;
fclose(fid);
asdf
%}
%--------------------------------------------------------------------------
% Figure A3(b): purchase utility and optimal down versus min down
%--------------------------------------------------------------------------
%{
numpts = 1;
values = zeros(numpts,5);
simulated_apps = [];
simulated_sales = [];
for i=1:numpts,

    %calculate simulated outcomes conditional on sale
    [dopt lopt uout ubuy uopt simsale simless simxtra simdef simfrac] = optimal_down(0.1*(i-1),ptilda,v0,vi,y0,yi,yindex,EVP,EVD,D,cval,betap,theta,beta,phi);
    x1 = mean(simsale);
    x2 = mean(simless(simsale));
    x3 = mean(dopt(simsale));
    x4 = mean(simfrac);
    x5 = mean(simdef);
    values(i,:)=[x1 x2 x3 x4 x5];
    %{
    %caluclate simulated outcomes not conditional on sale
    temp_d = -99; %no minimum down restriction
    [dopt lopt uout ubuy uopt simsale simless simxtra simdef simfrac] = optimal_down(temp_d,ptilda,v0,vi,y0,yi,yindex,EVP,EVD,D,cval,betap,theta,beta,phi);
    temp = dopt; temp(temp<0)=0;
    x6 = mean(temp);
    x7 = mean(simfrac);
    x8 = mean(simfrac);
    values(i,:)=[x1 x2 x3 x4 x5 x6 x7 x8];
    %}
    %create data for reduced-form regressions
    new_apps = [y0 ptilda*ones(Nsim,1) 0.1*(i-1)*ones(Nsim,1) simsale simless];
    new_sales = [new_apps(simsale,:) dopt(simsale) simdef simfrac];
    simulated_apps = [simulated_apps; new_apps];
    simulated_sales = [simulated_sales; new_sales];
    
end;

%--------------------------------------------------------------------------
% Output data
%--------------------------------------------------------------------------
fid = fopen('simoutcomes_vs_mindp.txt','wt');
for i=1:numpts,
    fprintf(fid,'%10.6f ',values(i,:));
    fprintf(fid,'\n');
end;
fclose(fid);
fid = fopen('simulated_apps_d.txt','wt');
for i=1:size(simulated_apps,1),
    fprintf(fid,'%10.6f,',simulated_apps(i,:));
    fprintf(fid,'\n');
end;
fclose(fid);
fid = fopen('simulated_sales_d.txt','wt');
for i=1:size(simulated_sales,1),
    fprintf(fid,'%10.6f,',simulated_sales(i,:));
    fprintf(fid,'\n');
end;
fclose(fid);
%}
%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------
