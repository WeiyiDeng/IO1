% square all premiums in SE calculations
clear all; clc;

load jointOut;

% adjust coefficients

premium = [-beta(19) se(19)]/100;

sigmaEpsilon = pi/sqrt(6);

% se by the delta-method g(theta) = sigmaEpsilon/premium(1)
% se(g(theta)) = premium(2) sigmaEpsilon/premium(1)^2

choice(23).string = '$\sigma_{\epsilon}$';
choice(23).beta = sigmaEpsilon/premium(1);
choice(23).se = (sigmaEpsilon*premium(2))/(premium(1)^2);

choice(1).string = 'premium';
choice(1).beta = -1;
choice(1).se = premium(2)/premium(1);

choice(2).string = 'Coinsurance';
choice(2).beta = beta(20)/(premium(1)*100);
choice(2).se = se(20)/(premium(1)*100);

choice(3).string = 'Deductible';
choice(3).beta = beta(21)/(premium(1)*1000);
choice(3).se = se(21)/(premium(1)*1000);

choice(4).string = 'Ph nonstandard';
choice(4).beta = beta(22)/(premium(1));
choice(4).se = se(22)/(premium(1));

choice(5).string = 'AHMO*Risk';
choice(5).beta = beta(4)/premium(1);
choice(5).se = se(4)/premium(1);

choice(6).string = 'AHMO*Age';
choice(6).beta = beta(10)/(premium(1)*100);
choice(6).se = se(10)/(premium(1)*100);

choice(7).string = 'AHMO*Male';
choice(7).beta = beta(13)/premium(1);
choice(7).se = se(13)/premium(1);

choice(8).string = 'AHMO*Income';
choice(8).beta = beta(16)/(premium(1)*100000);
choice(8).se = se(16)/(premium(1)*100000);

choice(9).string = 'AHMO*Risk90';
choice(9).beta = beta(7)/premium(1);
choice(9).se = se(7)/premium(1);

choice(10).string = 'AHMO*Constant';
choice(10).beta = beta(1)/premium(1);
choice(10).se = se(1)/premium(1);

choice(11).string = 'APPO*Risk';
choice(11).beta = beta(5)/premium(1);
choice(11).se = se(5)/premium(1);

choice(12).string = 'APPO*Age';
choice(12).beta = beta(11)/(premium(1)*100);
choice(12).se = se(11)/(premium(1)*100);

choice(13).string = 'APPO*Male';
choice(13).beta = beta(14)/premium(1);
choice(13).se = se(14)/premium(1);

choice(14).string = 'APPO*Income';
choice(14).beta = beta(17)/(premium(1)*100000);
choice(14).se = se(17)/(premium(1)*100000);

choice(15).string = 'APPO*Risk90';
choice(15).beta = beta(8)/premium(1);
choice(15).se = se(8)/premium(1);

choice(16).string = 'APPO*Constant';
choice(16).beta = beta(2)/premium(1);
choice(16).se = se(2)/premium(1);

choice(17).string = 'BPOS*Risk';
choice(17).beta = beta(6)/premium(1);
choice(17).se = se(6)/premium(1);

choice(18).string = 'BPOS*Age';
choice(18).beta = beta(12)/(premium(1)*100);
choice(18).se = se(12)/(premium(1)*100);

choice(19).string = 'BPOS*Male';
choice(19).beta = beta(15)/premium(1);
choice(19).se = se(15)/premium(1);

choice(20).string = 'BPOS*Income';
choice(20).beta = beta(18)/(premium(1)*100000);
choice(20).se = se(18)/(premium(1)*100000);

choice(21).string = 'BPOS*Risk90';
choice(21).beta = beta(9)/premium(1);
choice(21).se = se(9)/premium(1);

choice(22).string = 'BPOS*Constant';
choice(22).beta = beta(3)/premium(1);
choice(22).se = se(3)/premium(1);

fid = fopen('out.txt','w');

fprintf(fid,'\n\n%s\n\n','Choice');

for i = 1:length(choice)

    fprintf(fid,'%s\t',choice(i).string);
    fprintf(fid,'%6.4f\n',choice(i).beta);
    fprintf(fid,'\t\t%6.4f\n',choice(i).se);

end

fclose(fid);


