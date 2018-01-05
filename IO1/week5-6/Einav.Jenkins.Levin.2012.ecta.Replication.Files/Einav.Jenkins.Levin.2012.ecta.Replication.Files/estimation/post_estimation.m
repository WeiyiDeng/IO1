%--------------------------------------------------------------------------
% post_estimation.m
%--------------------------------------------------------------------------

clear all;
diary output.log;

global N X1 X2 X3 X4 X5 X6 K1 K2 K3 K4;
global I list less more none cenind obsdef cenpoint;
global sale default fracpaid appmonth;
global pK lK pstarK dK dK2 tK afK cK APRUSE APRMAX TERM S;
global gnum g appyear censored car_char simeps simnu;
global SLTXRT docfeeK slsfeeK ovrallK tlsfeeK;
global nu xi ep et mxi vxi mep vep met vet;
global dpperiod mdK mperiod lpK; % for price improvement analysis
global G;
global dpar subsample sale_curve prof_curve dmat;

load 'full_data';
par = load('par_final.txt');

%--------------------------------------------------------------------------
% create simulated data set
%--------------------------------------------------------------------------

M = 100; N = N*M; %number of simulations

X1       = repmat(X1,M,1);
X2       = repmat(X2,M,1);
X3       = repmat(X3,M,1);
X4       = repmat(X4,M,1);
lK       = repmat(lK,M,1);
pK       = repmat(pK,M,1);
dK       = repmat(dK,M,1);
tK       = repmat(tK,M,1);
afK      = repmat(afK,M,1);
I        = repmat(I,M,1);
less     = repmat(less,M,1);
more     = repmat(more,M,1);
none     = repmat(none,M,1);
obsdef   = repmat(obsdef,M,1);
cenind   = repmat(cenind,M,1);
fracpaid = repmat(fracpaid,M,1);
cenpoint = repmat(cenpoint,M,1);
censored = repmat(censored,M,1);
default  = repmat(default,M,1);
X5       = repmat(X5,M,1);
X6       = repmat(X6,M,1);
APRUSE   = repmat(APRUSE,M,1);
TERM     = repmat(TERM,M,1);
SLTXRT   = repmat(SLTXRT,M,1);
docfeeK  = repmat(docfeeK,M,1);
slsfeeK  = repmat(slsfeeK,M,1);
tlsfeeK  = repmat(tlsfeeK,M,1);
appmonth = repmat(appmonth,M,1);
gnum     = repmat(gnum,M,1);
cK       = repmat(cK,M,1);
dpperiod = repmat(dpperiod,M,1);
mdK      = repmat(mdK,M,1);

%--------------------------------------------------------------------------
% create tables and figures: DEMAND
%--------------------------------------------------------------------------
S = 10.0;
[nu xi ep et] = simulate_for_model_fit(par);

create_table2(pK,dK,par);            % Model Fit
create_table3(pK,dK,par);            % Marginal Effects
create_figure5(pK,dK,par);           % Model Fit
create_figure8(pK,dK,par,0);         % Down-Default Bubbles

%--------------------------------------------------------------------------
% Estimate supply-side parameter
%--------------------------------------------------------------------------

S = 10.0;
[nu xi ep et] = simulate_unobservables(par);
shadowcost = create_table5(par);

%--------------------------------------------------------------------------
% create tables and figures: SUPPLY
%--------------------------------------------------------------------------
[s down min frac def pvpmt pvrec crev ccost cprof urev uprof] = revenue_expectation(pK,dK,par,shadowcost);

create_table6(pK,dK,par,shadowcost);   % Value of Scoring
create_table7(pK,dK,par,shadowcost);   % Competition
create_figure6(pK,dK,par,shadowcost);  % Counterfactual Downs
create_figure7(pK,dK,par,shadowcost);  % Counterfactual Prices

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------


