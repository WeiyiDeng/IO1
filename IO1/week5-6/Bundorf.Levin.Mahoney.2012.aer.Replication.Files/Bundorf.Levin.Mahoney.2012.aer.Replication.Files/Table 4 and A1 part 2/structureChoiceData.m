clear all; clc;

cd '/Users/nmahoney/Documents/Jon/Joint/IV';
load choiceIndex;

I = length(data);
J = max(data(:,19));

index = ones(J,1);

count = 0;

for i=1:I
    
    cInd = data(i,19);
    ind = index(cInd,1);
    
    choiceData(cInd).group(ind)=data(i,1); %#ok<AGROW>
    choiceData(cInd).year(ind)=data(i,2); %#ok<AGROW>
    choiceData(cInd).tier(ind)=data(i,3); %#ok<AGROW>
    choiceData(cInd).fam(ind)=data(i,4); %#ok<AGROW>
    choiceData(cInd).ph(ind)=data(i,5); %#ok<AGROW>
    choiceData(cInd).bid(ind)=data(i,6)/100; %#ok<AGROW>
    choiceData(cInd).premium(ind)=data(i,7)/100; %#ok<AGROW>
    choiceData(cInd).bidHatDem(ind)=data(i,8)/100; %#ok<AGROW>
    choiceData(cInd).bidHatPro(ind)=data(i,9)/100; %#ok<AGROW>
    choiceData(cInd).dchoice(ind)=data(i,10); %#ok<AGROW>
    choiceData(cInd).co(ind)=data(i,11)/100; %#ok<AGROW>
    choiceData(cInd).ded(ind)=data(i,12)/1000; %#ok<AGROW>
    choiceData(cInd).premiumAT(ind)=data(i,13)/100; %#ok<AGROW> 
    choiceData(cInd).premiumIV(ind)=data(i,14)/100; %#ok<AGROW>
    choiceData(cInd).bidAT(ind)=data(i,15)/100; %#ok<AGROW>
    choiceData(cInd).bidHatDemAT(ind)=data(i,16)/100; %#ok<AGROW>
    choiceData(cInd).bidHatProAT(ind)=data(i,17)/100; %#ok<AGROW>
    choiceData(cInd).plan(ind)=data(i,18); %#ok<AGROW>

    if choiceData(cInd).bidHatDem(ind)==99.99
        
        count = count+1;
        choiceData(cInd).premiumIV(ind)=99.99; %#ok<AGROW>
        
    end
    
    
    for k=1:4

        choiceData(cInd).planX(ind,k)=data(i,19+k); %#ok<AGROW>
        choiceData(cInd).proX(ind,k)=data(i,23+k); %#ok<AGROW>
        choiceData(cInd).pro90X(ind,k)=data(i,27+k); %#ok<AGROW>
        choiceData(cInd).ageX(ind,k)=data(i,31+k)/100; %#ok<AGROW>
        choiceData(cInd).maleX(ind,k)=data(i,35+k); %#ok<AGROW>
        choiceData(cInd).incX(ind,k)=data(i,39+k)/100000; %#ok<AGROW>

%         for l=1:20
% 
%             choiceData(cInd).mu(ind,k,l)=data(i,31+4*l+k);  %#ok<AGROW>
%             
%         end
    
    end

    index(cInd,1) = index(cInd,1)+1;
        
end

save structureChoiceDataIV_3 choiceData; 
                
                
                
                
    
    
