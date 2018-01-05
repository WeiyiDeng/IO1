clear all; clc;

load jointOutMuHighRisk

% adjust coefficients

premium = [-beta(19) se(19)]/100;

sigmaEpsilon = pi/sqrt(6);

%se by the delta-method g(theta) = sigmaEpsilon/premium(1)
%se(g(theta)) = premium(2) sigmaEpsilon/premium(1)^2

choice(1).string = 'Contribution';
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

choice(5).string = 'NHMO*Risk';
choice(5).beta = beta(4)/premium(1);
choice(5).se = se(4)/premium(1);

choice(6).string = 'NHMO*Age';
choice(6).beta = beta(10)/(premium(1)*100);
choice(6).se = se(10)/(premium(1)*100);

choice(7).string = 'NHMO*Male';
choice(7).beta = beta(13)/premium(1);
choice(7).se = se(13)/premium(1);

choice(8).string = 'NHMO*Income';
choice(8).beta = beta(16)/(premium(1)*100000);
choice(8).se = se(16)/(premium(1)*100000);

choice(9).string = 'NHMO*Risk90';
choice(9).beta = beta(7)/premium(1);
choice(9).se = se(7)/premium(1);

choice(10).string = 'NHMO*Constant';
choice(10).beta = beta(1)/premium(1);
choice(10).se = se(1)/premium(1);

choice(11).string = 'NPPO*Risk';
choice(11).beta = beta(5)/premium(1);
choice(11).se = se(5)/premium(1);

choice(12).string = 'NPPO*Age';
choice(12).beta = beta(11)/(premium(1)*100);
choice(12).se = se(11)/(premium(1)*100);

choice(13).string = 'NPPO*Male';
choice(13).beta = beta(14)/premium(1);
choice(13).se = se(14)/premium(1);

choice(14).string = 'NPPO*Constant';
choice(14).beta = beta(2)/premium(1);
choice(14).se = se(2)/premium(1);

choice(15).string = 'NPPO*Income';
choice(15).beta = beta(17)/(premium(1)*100000);
choice(15).se = se(17)/(premium(1)*100000);

choice(16).string = 'NPPO*Risk90';
choice(16).beta = beta(8)/premium(1);
choice(16).se = se(8)/premium(1);

choice(17).string = 'IPOS*Risk';
choice(17).beta = beta(6)/premium(1);
choice(17).se = se(6)/premium(1);

choice(18).string = 'IPOS*Age';
choice(18).beta = beta(12)/(premium(1)*100);
choice(18).se = se(12)/(premium(1)*100);

choice(19).string = 'IPOS*Male';
choice(19).beta = beta(15)/premium(1);
choice(19).se = se(15)/premium(1);

choice(20).string = 'IPOS*Income';
choice(20).beta = beta(18)/(premium(1)*100000);
choice(20).se = se(18)/(premium(1)*100000);

choice(21).string = 'IPOS*Risk90';
choice(21).beta = beta(9)/premium(1);
choice(21).se = se(9)/premium(1);

choice(22).string = 'NPPO*Constant';
choice(22).beta = beta(3)/premium(1);
choice(22).se = se(3)/premium(1);

choice(23).string = '$\sigma_{\epsilon}$';
choice(23).beta = sigmaEpsilon/premium(1);
choice(23).se = (sigmaEpsilon*premium(2))/(premium(1)^2);

choice(24).string = '$\sigma_{\mu}$';
choice(24).beta = beta(33);
choice(24).se = se(33)/20;

costs(1).string = '$\delta_A$';
costs(1).beta = beta(23);
costs(1).se = se(23);

costs(2).string = '$\delta_B$';
costs(2).beta = beta(24);
costs(2).se = se(24);

costs(3).string = '$\alpha_{NHMO}$';
costs(3).beta = beta(25);
costs(3).se = se(25);

costs(4).string = '$\alpha_{NPPO}$';
costs(4).beta = beta(26);
costs(4).se = se(26);

costs(5).string = '$\alpha_{IHMO}$';
costs(5).beta = beta(27);
costs(5).se = se(27);

costs(6).string = '$\alpha_{IPOS}$';
costs(6).beta = beta(28);
costs(6).se = se(28);

costs(7).string = '$\beta_{NHMO}$';
costs(7).beta = beta(29);
costs(7).se = se(29);

costs(8).string = '$\beta_{NPPO}$';
costs(8).beta = beta(30);
costs(8).se = se(30);

costs(9).string = '$\beta_{IHMO}$';
costs(9).beta = beta(31);
costs(9).se = se(31);

costs(10).string = '$\beta_{IPOS}$';
costs(10).beta = beta(32);
costs(10).se = se(32);

fid = fopen('outMuHighRisk.txt','w');

fprintf(fid,'\n\n%s\n\n','Choice');

for i = 1:length(choice)

    fprintf(fid,'%s\t',choice(i).string);
    fprintf(fid,'%6.4f\n',choice(i).beta);
    fprintf(fid,'\t%6.4f\n',choice(i).se);

end

fprintf(fid,'\n\n%s\n\n','Costs');

for i = 1:length(costs)

    fprintf(fid,'%s\t',costs(i).string);
    fprintf(fid,'%6.4f\n',costs(i).beta);
    fprintf(fid,'\t\t%6.4f\n',costs(i).se);

end

fclose(fid);


