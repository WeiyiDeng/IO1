%--------------------------------------------------------------------------
% create_table3.m
%--------------------------------------------------------------------------

function [salemfx downmfx deftmfx] = create_table3(pKtemp,dKtemp,par)

global N X1 X2 X3 X4 K1 K2 K3 K4 lK;
global nu xi ep et mxi vxi mep vep met vet;
global SLTXRT docfeeK slsfeeK tlsfeeK;

%define parameters
[Lambda Delta Beta Gamma] = define_parameters(par);

%--------------------------------------------------------------------------
% Compute marginal effects of explanatory variables on pr(sale/mind/def)
%--------------------------------------------------------------------------
MarginalEffects = zeros(K1,1);
for i=1:K1,
    %create matrices at high and low X values
    if (max(X1(:,i))==1 && min(X1(:,i)==0)) %check for dummy variable
        step = 1;
        XL = X1; XL(:,i)=0;
        XH = X1; XH(:,i)=1;
    else
        step = 0.1;
        XL = X1; XL(:,i)=X1(:,i);
        XH = X1; XH(:,i)=X1(:,i)+step;
    end;

    %calculate negotiated price
    pricL = [XL lK]*Lambda + nu;
    pricH = [XH lK]*Lambda + nu;
    
    %calculate down payment as a function of minimum down
    downL = max([XL pKtemp]*Beta + ep, dKtemp);
    downH = max([XH pKtemp]*Beta + ep, dKtemp);

    %calculate minimum down payment indicator
    mindL = downL <= dKtemp;
    mindH = downH <= dKtemp;

    %calculate sale indicator
    dKtemp2 = dKtemp.*dKtemp;
    saleL = ([XL pricL dKtemp dKtemp2]*Delta + xi > 0);
    saleH = ([XH pricH dKtemp dKtemp2]*Delta + xi > 0);
    
    %calculate amount financed
    amtfL = (pricL.*(1+SLTXRT) + docfeeK + slsfeeK + tlsfeeK)-downL;
    amtfH = (pricH.*(1+SLTXRT) + docfeeK + slsfeeK + tlsfeeK)-downH;

    %calculate fraction of payments made
    fracL = exp([XL amtfL]*Gamma + et); deftL = fracL<1;
    fracH = exp([XH amtfH]*Gamma + et); deftH = fracH<1;
    
    %calculate marginal effects
    MarginalEffects(i,1) = (mean(saleH) - mean(saleL))/step;
    MarginalEffects(i,2) = (mean(mindH(saleH)) - mean(mindL(saleL)))/step;
    MarginalEffects(i,3) = (mean(deftH(saleH)) - mean(deftL(saleL)))/step;

end;


%--------------------------------------------------------------------------
% output data
%--------------------------------------------------------------------------

fid = fopen('table3.txt','wt');
for i=1:K1,
    fprintf(fid,'%6.4f  ',MarginalEffects(i,:));
    fprintf(fid,'\n');
end
fclose(fid);

%--------------------------------------------------------------------------
% create output variables for bootstrap
%--------------------------------------------------------------------------

salemfx = MarginalEffects(:,1);
downmfx = MarginalEffects(:,2);
deftmfx = MarginalEffects(:,3);

%--------------------------------------------------------------------------
% end of program
%--------------------------------------------------------------------------
