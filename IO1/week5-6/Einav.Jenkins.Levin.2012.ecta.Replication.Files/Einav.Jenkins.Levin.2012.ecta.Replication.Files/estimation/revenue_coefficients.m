%--------------------------------------------------------------------------
% revenue_coefficients.m
%--------------------------------------------------------------------------

function [pvlpratio pvrecovery] = revenue_coefficients(fractemp,S,R,T,X5,X6)

global N G Theta1 Theta2;

% coefficients for nonzero recovery probit
Theta1 = [-0.0572433; -0.0194657; 0.0904141; -0.3001009; -0.0283294; 0.0345648; -0.0139271; -0.0013141; 0.0046814; 0.0038578; 0.1283485; 0.0990908; 0.0919281; 0.1064601; 0.098366; 0.0738443; 0.0653585; 0.0490054; 0.030762; 0.0454848; -0.0093926; 0.2774221; 0.420175; 0.4020728; 0.3016064; 0.4410063; 0.3113178; 0.4624002; 0.3882629; 0.2330741; 0.4787957; 0.3404164];

% coefficients for recovery amount regression
Theta2 = [1.551523; -0.0452309; 0.5260421; 0.0792684; -0.0675788; -0.0291352; -0.131353; 0.0478444; 0.0779298; 0.1180724; -0.0501646; -0.1370988; -0.1319799; -0.1114999; -0.0831248; -0.2037505; -0.1784734; -0.1455426; -0.1289853; -0.1363506; -0.1079046; -1.331157; -1.215872; -1.300575; -1.042137; -1.163661; -1.275036; -1.016728; -1.131995; -1.413663; -1.237001; -1.291203];

%--------------------------------------------------------------------------
% accounting definitions
%--------------------------------------------------------------------------
s  = S / 1200;    % firm internal discount rate (monthly)
r  = R / 1200;    % customer's interest rate (monthly)
sg = (T/G).*s;    % firm internal discount rate (per grid point)
rg = (T/G).*r;    % customer's interest rate (per grid point)

%--------------------------------------------------------------------------
% coefficients for expected revenue equation
%--------------------------------------------------------------------------
%discretize fracpaid into grid with G points
g = floor(fractemp*G); paid = fractemp==1; 

%calculate ratio of pv loan payments to loan size
pvfrac = zeros(N,1);
for t=1:G,
    pmt_t  = (1+sg).^(-t).*(rg./(1-(ones(N,1)+rg).^(-G)));
    pvfrac = pvfrac + (t<=g).*pmt_t;
end;
pvlpratio = pvfrac;

%recovery coefficient for default customers
monthspaid = (T/G).*g;
XTheta1 = [ones(N,1) monthspaid X5]*Theta1;
XTheta2 = [ones(N,1) monthspaid X6]*Theta2;
pvrecovery = (1+sg).^(-g).*normcdf(XTheta1,0,1).*XTheta2;
pvrecovery = paid.*0 + ~paid.*pvrecovery;

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------
