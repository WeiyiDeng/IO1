% computed from conditional sample frequency
function DeltaU = estimateDeltaU(choices,iX,nSuppX,nFirms,nPeriods)

count_At0 = zeros(nSuppX,2);
count_At1 = zeros(nSuppX,2);
count_Xt_AtPrev = zeros(nSuppX,2);

for i = 1:nFirms
    for t = 2:nPeriods
        if choices(t,i)==0
            count_At0(iX(t,i),choices(t-1,i)+1) = count_At0(iX(t,i),choices(t-1,i)+1)+1;
        else
            count_At1(iX(t,i),choices(t-1,i)+1) = count_At1(iX(t,i),choices(t-1,i)+1)+1;
        end
        count_Xt_AtPrev(iX(t,i),choices(t-1,i)+1) = count_Xt_AtPrev(iX(t,i),choices(t-1,i)+1)+1;
    end
end
   
ConditionalPr_At0 = count_At0./count_Xt_AtPrev;
ConditionalPr_At1 = count_At1./count_Xt_AtPrev;

DeltaU = log(ConditionalPr_At1)-log(ConditionalPr_At0);