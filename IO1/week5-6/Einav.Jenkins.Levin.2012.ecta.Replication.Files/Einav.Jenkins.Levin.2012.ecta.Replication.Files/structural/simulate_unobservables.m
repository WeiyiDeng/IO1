%--------------------------------------------------------------------------
% Outputs simulated unobservables v0 and y0 through yT
%--------------------------------------------------------------------------

function [vf vi yf yi yindex yactual] = simulate_unobservables(vavg,vvar,davg,dvar,corr,sig2)

global saleeps downeps epsilon T;
global vmin vmax vstep;
global ymin ymax ystep; 

%simulate car value (normally distributed)
vf = vavg + (vvar^0.5)*saleeps;
vd  = max(vmin,min(vmax,round((vf-vmin)*(1/vstep))*vstep + vmin));
vd(vf<vmin)=vmin; vd(vf>vmax)=vmax;
vi  = round((vd - vmin)/vstep) + 1;

%simulate initial liquidity (normally distributed)
cavg = davg + corr*(dvar/vvar)^0.5*saleeps;
cvar = dvar*(1-corr^2);
yf   = cavg + (cvar^0.5)*downeps;
yd  = max(ymin,min(ymax,round((yf-ymin)*(1/ystep))*ystep + ymin));
yd(yf<ymin)=ymin; yd(yf>ymax)=ymax;
yi  = round((yd - ymin)/ystep) + 1;

%simulate liquidity paths (random walk); first line dimensions: (TxN x NxN)'=(TxN)'=NxT
income_shocks = (epsilon'*diag(sig2^0.5))';
income_shocks = [yf+income_shocks(:,1) income_shocks(:,2:T)];        
yactual = cumsum(income_shocks,2);
yround  = max(ymin,min(ymax,round((yactual-ymin)*(1/ystep))*ystep + ymin));
yindex  = round((yround - ymin)/ystep) + 1;

%--------------------------------------------------------------------------
% end of function
%--------------------------------------------------------------------------
