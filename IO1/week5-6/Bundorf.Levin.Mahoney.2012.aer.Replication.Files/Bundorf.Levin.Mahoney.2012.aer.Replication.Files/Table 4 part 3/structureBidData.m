clear all; clc;

cd '/Users/nmahoney/Documents/Jon/Joint/Data';
load bidIndex;

I = length(data);
J = max(data(:,9));

index = ones(J,1);

for i=1:I
    
    bInd = data(i,9);
    ind = index(bInd,1);
    
    bidData(bInd).group(ind)=data(i,1); %#ok<AGROW>
    bidData(bInd).year(ind)=data(i,2); %#ok<AGROW>
    bidData(bInd).carrier(ind)=data(i,3); %#ok<AGROW>
    bidData(bInd).proPredict(ind)=data(i,4); %#ok<AGROW>
    bidData(bInd).co(ind)=data(i,5); %#ok<AGROW>
    bidData(bInd).ded(ind)=data(i,6); %#ok<AGROW>
    bidData(bInd).plan(ind)=data(i,7); %#ok<AGROW>
    bidData(bInd).bid(ind)=data(i,8); %#ok<AGROW>

    index(bInd,1) = index(bInd,1)+1;
        
end

save ../structureBidData bidData
