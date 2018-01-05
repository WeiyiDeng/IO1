clear all; clc;

cd '/Users/nmahoney/Documents/Jon/Joint/Data';
load choiceIndex;

I = length(data);
J = max(data(:,41));

index = ones(J,1);

for i=1:I

    cInd = data(i,41);
    ind = index(cInd,1);

    choiceData(cInd).group(ind)=data(i,1); %#ok<AGROW>
    choiceData(cInd).year(ind)=data(i,2); %#ok<AGROW>
    choiceData(cInd).tier(ind)=data(i,4); %#ok<AGROW>
    choiceData(cInd).maxPro(ind)=data(i,5); %#ok<AGROW>

    for l = 1:20

        choiceData(cInd).maxMu(ind,l)=data(i,5+l);  %#ok<AGROW>

    end

    choiceData(cInd).fam(ind)=data(i,26); %#ok<AGROW>
    choiceData(cInd).ph(ind)=data(i,27); %#ok<AGROW>
    choiceData(cInd).bid(ind)=data(i,28)/100; %#ok<AGROW>
    choiceData(cInd).premium(ind)=data(i,29)/100; %#ok<AGROW>
    choiceData(cInd).bidHatDem(ind)=data(i,30)/100; %#ok<AGROW>
    choiceData(cInd).bidHatPro(ind)=data(i,31)/100; %#ok<AGROW>
    choiceData(cInd).dchoice(ind)=data(i,32); %#ok<AGROW>
    choiceData(cInd).co(ind)=data(i,33)/100; %#ok<AGROW>
    choiceData(cInd).ded(ind)=data(i,34)/1000; %#ok<AGROW>
    choiceData(cInd).premiumAT(ind)=data(i,35)/100; %#ok<AGROW>
    choiceData(cInd).premiumIV(ind)=data(i,36)/100; %#ok<AGROW>
    choiceData(cInd).bidAT(ind)=data(i,37)/100; %#ok<AGROW>
    choiceData(cInd).bidHatDemAT(ind)=data(i,38)/100; %#ok<AGROW>
    choiceData(cInd).bidHatProAT(ind)=data(i,39)/100; %#ok<AGROW>
    choiceData(cInd).plan(ind)=data(i,40); %#ok<AGROW>


    for k=1:4

        choiceData(cInd).planX(ind,k)=data(i,41+k); %#ok<AGROW>
        choiceData(cInd).proX(ind,k)=data(i,45+k); %#ok<AGROW>
        choiceData(cInd).pro90X(ind,k)=data(i,49+k); %#ok<AGROW>
        choiceData(cInd).ageX(ind,k)=data(i,53+k)/100; %#ok<AGROW>
        choiceData(cInd).maleX(ind,k)=data(i,57+k); %#ok<AGROW>
        choiceData(cInd).incX(ind,k)=data(i,61+k)/100000; %#ok<AGROW>

        for l=1:20

            choiceData(cInd).mu(ind,k,l)=data(i,61+4*l+k);  %#ok<AGROW>

        end

    end

    index(cInd,1) = index(cInd,1)+1;

end

save ../structureChoiceData choiceData






