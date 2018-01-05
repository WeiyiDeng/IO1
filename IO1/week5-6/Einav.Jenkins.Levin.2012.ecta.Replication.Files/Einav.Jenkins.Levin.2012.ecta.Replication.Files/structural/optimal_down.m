%--------------------------------------------------------------------------
% Inputs minimum down, price, y0, and v0; outputs optimal down payments
%--------------------------------------------------------------------------

function [dopt lopt uout ubuy uopt simsale simless simxtra simdef simfrac dsim] = optimal_down(d,p,v,vi,y,yi,yindex,EVP,EVD,D,cval,betap,theta,beta,phi)

global T G APRUSE TERM;
global lstates lnum lmin lstep;
global vmin vmax vstep;
global ymin ymax ystep;

%--------------------------------------------------------------------------
% compute future values from purchase with smoothed income and car value
%--------------------------------------------------------------------------
vf = v; yf = y;

%identify nearest car value gridpoints
vdL = floor((vf-vmin)*(1/vstep))*vstep + vmin;
vdL(vf<vmin)=vmin; vdL(vf>vmax)=vmax-vstep;
viL = round((vdL - vmin)/vstep) + 1;
vdH = ceil((vf-vmin)*(1/vstep))*vstep + vmin;
vdH(vf<vmin)=vmin+vstep; vdH(vf>vmax)=vmax;
viH = round((vdH - vmin)/vstep) + 1;

%identify nearest liquidity gridpoints
ydL = floor((yf-ymin)*(1/ystep))*ystep + ymin;
ydL(yf<ymin)=ymin; ydL(yf>ymax)=ymax-ystep;
yiL = round((ydL - ymin)/ystep) + 1;
ydH = ceil((yf-ymin)*(1/ystep))*ystep + ymin;
ydH(yf<ymin)=ymin+ystep; ydH(yf>ymax)=ymax;
yiH = round((ydH - ymin)/ystep) + 1;

%compute weights for nearest gridpoints
vtemp = max(vmin,min(vf,vmax));
ytemp = max(ymin,min(yf,ymax));
wLL = ((vdH-vtemp)/vstep).*((ydH-ytemp)/ystep);
wHL = ((vtemp-vdL)/vstep).*((ydH-ytemp)/ystep);
wLH = ((vdH-vtemp)/vstep).*((ytemp-ydL)/ystep);
wHH = ((vtemp-vdL)/vstep).*((ytemp-ydL)/ystep);

%reshape weights into matrix form
wLL=repmat(wLL,1,lnum); wHL=repmat(wHL,1,lnum);
wLH=repmat(wLH,1,lnum); wHH=repmat(wHH,1,lnum);

%compute future values at nearest gridpoints
temp = EVP(:,:,:,1);
evpLL = CalcEVP1(viL,yiL,temp);
evpHL = CalcEVP1(viH,yiL,temp);
evpLH = CalcEVP1(viL,yiH,temp);
evpHH = CalcEVP1(viH,yiH,temp);

%compute smoothed future values
%evp1a = CalcEVP1(vi,yi,temp);
evp1b = (wLL.*evpLL + wHL.*evpHL + wLH.*evpLH + wHH.*evpHH);

%--------------------------------------------------------------------------
% compute optimal down payment and repayment decisions
%--------------------------------------------------------------------------
%calculate polynomial EVP1(L) function
x = [ones(lnum,1) lstates' (lstates.^2)']; %Lx3
xx = x'*x;      %3x3
xy = evp1b*x;   %Nx3
coef = xy/xx;   %Nx3
a = coef(:,3); b = coef(:,2); c = coef(:,1);

%calculate first-order condition for down choice
T = G; R = APRUSE/1200.*(TERM/G); F = (1 ./ R) .* (1 - (1 + R).^-T);
d1 = 2*F*betap*a;
d2 = F*betap*b + 2*betap*a.*(y-p);
d3 = theta*F + betap*b.*(y-p);

%check for zeros of first-order condition
disc = d2.^2 - 4*d1.*d3; rflag = disc >= 0; sqrtdisc = disc.^0.5;
if (size(p,1)==1),
    p = p*ones(size(y,1),1);
end;
sqrtdisc(~rflag) = -(2*d1(~rflag).*(p(~rflag)./F)+d2(~rflag));

%compute utility at interior solution
lint = (-d2 - sqrtdisc)./(2*d1);
dint = p - F*lint;

%compute utility with minimum down req.
dopt = max(dint,d);
lopt = (p - dopt)./F;
temp = CalcUtility(y,dopt);
ubuy = v + theta*temp + betap*(a.*lopt.^2 + b.*lopt + c);

%disqualify down payments < liquidity
dopt(y<=d) = -Inf; lopt(y<=d) = -Inf; ubuy(y<=d) = -Inf;

%compute future value of no purchase at nearest gridpoints
temp = EVD(:,:,:,1);
evdLL = CalcEVP1(viL,yiL,temp);
evdHL = CalcEVP1(viH,yiL,temp);
evdLH = CalcEVP1(viL,yiH,temp);
evdHH = CalcEVP1(viH,yiH,temp);
evd1b = (wLL.*evdLL + wHL.*evdHL + wLH.*evdLH + wHH.*evdHH);
evd1b = evd1b + sum(phi*beta.^(1:T));

%disqulify down payments if no sale
temp = CalcUtility(y,0);
uout = cval + theta*temp + betap*evd1b(:,1);
dopt(ubuy<uout) = -Inf; lopt(ubuy<uout) = -Inf;
uopt = max(ubuy,uout);

%calculate simulated purchase and down outcomes
simsale  = dopt >= 0;
simxtra  = dopt - d;
simless  = simxtra<=0 & simsale;

%compute repayment paths using optimal decision rules
simpmti  = round((lopt - lmin)/lstep) + 1;
simpmti(simpmti<1) = 1;
dsim = CalcEVPt(simpmti(simsale),vi(simsale),yindex(simsale,:),D);
dsim = cumprod(dsim,2);

%calculate simulated repayment moments
simad = 1-dsim(:,T);
simdef = simad==1;
simfrac = sum(dsim,2)/G;

%--------------------------------------------------------------------------
% End of function
%--------------------------------------------------------------------------
