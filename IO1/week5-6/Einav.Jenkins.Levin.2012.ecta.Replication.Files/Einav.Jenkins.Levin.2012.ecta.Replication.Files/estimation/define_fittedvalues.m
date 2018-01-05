%--------------------------------------------------------------------------
% define_fittedvalues.m
%--------------------------------------------------------------------------

function [XLtemp XDtemp XBtemp XGtemp ptilda] = define_fittedvalues(p,d,par)

global X1 X2 X3 X4 lK SLTXRT docfeeK slsfeeK tlsfeeK;

% define parameters
[Lambda Delta Beta Gamma] = define_parameters(par);

% calcluate "effective" price used in amount financed and net revenue
ptilda = p.*(1+SLTXRT) + docfeeK + slsfeeK + tlsfeeK;

% compute fitted values
XLtemp = [X1 lK]*Lambda;
XDtemp = [X4 p d]*Delta;
XBtemp = [X2 p]*Beta;
XGtemp = [X3 ptilda]*Gamma;

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------