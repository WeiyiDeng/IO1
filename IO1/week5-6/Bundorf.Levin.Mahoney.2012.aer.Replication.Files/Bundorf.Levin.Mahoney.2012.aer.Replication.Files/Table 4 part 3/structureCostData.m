clear all; clc;

cd '/Users/nmahoney/Documents/Jon/Joint/Data';
load costIndex;

I = length(data);
J = max(data(:,29));

index = ones(J,1);

for i=1:I
    
    cInd = data(i,29);
    ind = index(cInd,1);
    
    costData(cInd).group(ind)=data(i,1); %#ok<AGROW>
    costData(cInd).year(ind)=data(i,2); %#ok<AGROW>
    costData(cInd).cost(ind)=data(i,3); %#ok<AGROW>
    costData(cInd).carrier(ind)=data(i,4); %#ok<AGROW>
    
    for k=1:20
       costData(cInd).mu(ind,k)=data(i,4+k);  %#ok<AGROW>
    end
    
    costData(cInd).pro(ind)=data(i,25); %#ok<AGROW>
    costData(cInd).co(ind)=data(i,26); %#ok<AGROW>
    costData(cInd).ded(ind)=data(i,27); %#ok<AGROW>
    costData(cInd).plan(ind)=data(i,28);  %#ok<AGROW>

    index(cInd,1) = index(cInd,1)+1;
        
end

save ../structureCostData costData
                
                
                
                
    
    
