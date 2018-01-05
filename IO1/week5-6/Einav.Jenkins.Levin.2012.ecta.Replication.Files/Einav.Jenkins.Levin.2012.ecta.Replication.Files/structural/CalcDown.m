%--------------------------------------------------------------------------
% Compute down payment given price and per-period payment size
%--------------------------------------------------------------------------

function [down_payment] = CalcDown(m,ptemp)

global APRUSE TERM;
global G;

T = G; R = APRUSE/1200.*(TERM/G);
down_payment = ptemp - m ./ (R ./ (1 - (1 + R).^-T));

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------
